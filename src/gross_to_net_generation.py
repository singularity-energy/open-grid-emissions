import numpy as np
import os
import pandas as pd
import statsmodels.formula.api as smf
import sqlalchemy as sa
import warnings

# import pudl packages
import pudl.analysis.allocate_net_gen as allocate_gen_fuel
import pudl.output.pudltabl

# import other modules
import load_data
import data_cleaning
import validation
from column_checks import get_dtypes
from filepaths import outputs_folder, downloads_folder
from logging_util import get_logger

logger = get_logger(__name__)


def convert_gross_to_net_generation(cems, eia923_allocated, plant_attributes, year):
    """
    Converts hourly gross generation in CEMS to hourly net generation by calculating a gross to net generation ratio
    Inputs:

    Returns:
        cems df with an added column for net_generation_mwh and a column indicated the method used to calculate net generation
    """

    gtn_conversions = calculate_gross_to_net_conversion_factors(
        cems, eia923_allocated, plant_attributes, year
    )

    factors_to_use = filter_gtn_conversion_factors(gtn_conversions)

    # merge the conversion factors we want to use into the cems data
    cems = cems.merge(
        factors_to_use[
            [
                "plant_id_eia",
                "subplant_id",
                "report_date",
                "data_source",
                "annual_subplant_shift_mw",
                "annual_plant_shift_mw",
                "annual_subplant_ratio",
                "annual_plant_ratio",
                "annual_fuel_ratio",
                "default_gtn_ratio",
            ]
        ],
        how="left",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="m:1",
    )

    cems["gtn_method"] = "1_annual_subplant_ratio"
    cems["net_generation_mwh"] = (
        cems["gross_generation_mwh"] * cems["annual_subplant_ratio"]
    )

    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "2_annual_plant_ratio"
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["annual_plant_ratio"]
    )

    cems.loc[
        cems["net_generation_mwh"].isna(), "gtn_method"
    ] = "3_annual_subplant_shift_factor"
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] + cems["annual_subplant_shift_mw"]
    )

    cems.loc[
        cems["net_generation_mwh"].isna(), "gtn_method"
    ] = "4_annual_plant_shift_factor"
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] + cems["annual_plant_shift_mw"]
    )

    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "5_annual_fuel_ratio"
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["annual_fuel_ratio"]
    )

    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "6_default_eia_ratio"
    # warn if there are any missing default gtn ratios for plants that would use them.
    missing_defaults = cems.loc[
        (cems["gtn_method"] == "6_default_eia_ratio")
        & (cems["default_gtn_ratio"].isna())
    ]
    if len(missing_defaults) > 0:
        logger.warning(
            "The following subplants are missing default GTN ratios. Using a default value of 0.97"
        )
        logger.warning(
            "\n"
            + missing_defaults[["plant_id_eia", "subplant_id"]]
            .drop_duplicates()
            .to_string()
        )
    # if there is a missing default gtn ratio, fill with 0.97
    cems["default_gtn_ratio"] = cems["default_gtn_ratio"].fillna(0.97)
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["default_gtn_ratio"]
    )

    # drop intermediate columns
    cems = cems.drop(
        columns=[
            "data_source",
            "annual_subplant_shift_mw",
            "annual_plant_shift_mw",
            "annual_subplant_ratio",
            "annual_plant_ratio",
            "annual_fuel_ratio",
            "default_gtn_ratio",
        ]
    )

    validation.validate_gross_to_net_conversion(cems, eia923_allocated)

    return cems, gtn_conversions


def calculate_gross_to_net_conversion_factors(
    cems, eia923_allocated, plant_attributes, year
):
    """
    Calculates gross to net ratios and shift factors
    """
    # aggregate the hourly cems data by subplant
    gross_gen_data = cems[
        [
            "plant_id_eia",
            "subplant_id",
            "report_date",
            "datetime_utc",
            "gross_generation_mwh",
        ]
    ].copy()
    # identify the 2nd percentile lowest hourly gross generation value in a month
    min_gross = (
        gross_gen_data.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .agg({"gross_generation_mwh": lambda x: x.quantile(0.02)})
        .reset_index()
        .rename(columns={"gross_generation_mwh": "minimum_gross_generation_mwh"})
    )
    # identify the 98th percentile highest hourly gross generation value in a month
    max_gross = (
        gross_gen_data.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .agg({"gross_generation_mwh": lambda x: x.quantile(0.98)})
        .reset_index()
        .rename(columns={"gross_generation_mwh": "maximum_gross_generation_mwh"})
    )
    subplant_capacity = calculate_subplant_nameplate_capacity(year)
    # aggregate the cems data to the monthly level
    gross_gen_data = (
        gross_gen_data.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .agg({"datetime_utc": "count", "gross_generation_mwh": "sum"})
        .reset_index()
        .rename(columns={"datetime_utc": "hours_in_month"})
    )
    net_gen_data = (
        eia923_allocated.dropna(subset=["net_generation_mwh"])
        .groupby(["plant_id_eia", "subplant_id", "report_date"], dropna=False)[
            "net_generation_mwh"
        ]
        .sum()
        .reset_index()
    )

    # combine monthly gross and net generation data where we have data for both
    combined_gen_data = (
        gross_gen_data.merge(
            min_gross,
            how="left",
            on=["plant_id_eia", "subplant_id", "report_date"],
            validate="1:1",
        )
        .merge(
            max_gross,
            how="left",
            on=["plant_id_eia", "subplant_id", "report_date"],
            validate="1:1",
        )
        .merge(
            subplant_capacity,
            how="left",
            on=["plant_id_eia", "subplant_id"],
            validate="m:1",
        )
        .merge(
            net_gen_data,
            how="outer",
            on=["plant_id_eia", "subplant_id", "report_date"],
            indicator="data_source",
            validate="1:1",
        )
    )
    combined_gen_data["data_source"] = combined_gen_data["data_source"].replace(
        {"left_only": "cems_only", "right_only": "eia_only"}
    )

    # validate the data
    validation.check_missing_or_zero_generation_matches(combined_gen_data)

    # calculate other groupings at the plant and annual levels
    annual_subplant_ratio = (
        combined_gen_data.dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia", "subplant_id"], dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
        ]
        .sum()
        .reset_index()
    )
    monthly_plant_ratio = (
        combined_gen_data.dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia", "report_date"], dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh"]
        ]
        .sum()
        .reset_index()
    )
    annual_plant_ratio = (
        combined_gen_data.dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia"], dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
        ]
        .sum()
        .reset_index()
    )

    # calculate the ratios at each aggregation level
    # fill missing values (due to divide by zero) with zero
    # replace infinite values with missing
    combined_gen_data["monthly_subplant_ratio"] = (
        combined_gen_data["net_generation_mwh"]
        / combined_gen_data["gross_generation_mwh"]
    ).replace([np.inf, -np.inf], np.nan)
    annual_subplant_ratio["annual_subplant_ratio"] = (
        annual_subplant_ratio["net_generation_mwh"]
        / annual_subplant_ratio["gross_generation_mwh"]
    ).replace([np.inf, -np.inf], np.nan)
    monthly_plant_ratio["monthly_plant_ratio"] = (
        monthly_plant_ratio["net_generation_mwh"]
        / monthly_plant_ratio["gross_generation_mwh"]
    ).replace([np.inf, -np.inf], np.nan)
    annual_plant_ratio["annual_plant_ratio"] = (
        annual_plant_ratio["net_generation_mwh"]
        / annual_plant_ratio["gross_generation_mwh"]
    ).replace([np.inf, -np.inf], np.nan)

    # flag anomalous plant ratios
    validation.identify_anomalous_annual_plant_gtn_ratios(annual_plant_ratio)

    # calculate a monthly and annual shift factor
    combined_gen_data["monthly_subplant_shift_mw"] = (
        combined_gen_data["net_generation_mwh"]
        - combined_gen_data["gross_generation_mwh"]
    ) / (combined_gen_data["hours_in_month"])
    annual_subplant_ratio["annual_subplant_shift_mw"] = (
        annual_subplant_ratio["net_generation_mwh"]
        - annual_subplant_ratio["gross_generation_mwh"]
    ) / (annual_subplant_ratio["hours_in_month"])
    annual_plant_ratio["annual_plant_shift_mw"] = (
        annual_plant_ratio["net_generation_mwh"]
        - annual_plant_ratio["gross_generation_mwh"]
    ) / (annual_plant_ratio["hours_in_month"])

    # drop the gross and net generation data from the dataframes at teh other aggregation levels
    annual_subplant_ratio = annual_subplant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
    )
    monthly_plant_ratio = monthly_plant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh"]
    )
    annual_plant_ratio = annual_plant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh"]
    )

    # merge the various ratios back into a single dataframe
    gtn_conversions = combined_gen_data.merge(
        annual_subplant_ratio,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )
    gtn_conversions = gtn_conversions.merge(
        monthly_plant_ratio,
        how="left",
        on=["plant_id_eia", "report_date"],
        validate="m:1",
    )
    gtn_conversions = gtn_conversions.merge(
        annual_plant_ratio,
        how="left",
        on=["plant_id_eia"],
        suffixes=("_subplant", "_plant"),
        validate="m:1",
    )

    # where gross or net generation data was missing in a month, change the monthly ratios to missing
    gtn_conversions.loc[
        gtn_conversions[["gross_generation_mwh", "net_generation_mwh"]]
        .isna()
        .any(axis=1),
        ["monthly_subplant_ratio", "monthly_plant_ratio"],
    ] = np.NaN

    # calculate the mean ratio for all plants of a single fuel type
    annual_fuel_ratio = (
        annual_plant_ratio.merge(
            plant_attributes[["plant_id_eia", "plant_primary_fuel"]],
            how="left",
            on="plant_id_eia",
            validate="m:1",
        )
        .groupby("plant_primary_fuel")["annual_plant_ratio"]
        .mean()
        .reset_index()
        .rename(columns={"annual_plant_ratio": "annual_fuel_ratio"})
    )

    # merge the plant primary fuel and the fuel ratios into the conversion table
    gtn_conversions = gtn_conversions.merge(
        plant_attributes[["plant_id_eia", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )
    gtn_conversions = gtn_conversions.merge(
        annual_fuel_ratio, how="left", on="plant_primary_fuel", validate="m:1"
    )

    # add regression values
    gtn_regression_subplant = gross_to_net_regression(combined_gen_data, "subplant")
    gtn_regression_subplant = gtn_regression_subplant[
        gtn_regression_subplant["rsquared_adj"] > 0.9
    ]
    gtn_regression_subplant = gtn_regression_subplant.rename(
        columns={
            "slope": "subplant_regression_ratio",
            "intercept": "subplant_regression_shift_mw",
            "rsquared_adj": "subplant_regression_rsq_adj",
        }
    )
    gtn_regression_subplant = gtn_regression_subplant.drop(
        columns=["rsquared", "observations"]
    )

    gtn_regression_plant = gross_to_net_regression(combined_gen_data, "plant")
    gtn_regression_plant = gtn_regression_plant[
        gtn_regression_plant["rsquared_adj"] > 0.9
    ]
    gtn_regression_plant = gtn_regression_plant.rename(
        columns={
            "slope": "plant_regression_ratio",
            "intercept": "plant_regression_shift_mw",
            "rsquared_adj": "plant_regression_rsq_adj",
        }
    )
    gtn_regression_plant = gtn_regression_plant.drop(
        columns=["rsquared", "observations"]
    )

    gtn_conversions = gtn_conversions.merge(
        gtn_regression_subplant,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )
    gtn_conversions = gtn_conversions.merge(
        gtn_regression_plant, how="left", on=["plant_id_eia"], validate="m:1"
    )

    # add default GTN ratios from EIA
    gtn_conversions = gtn_conversions.merge(
        load_data.load_default_gtn_ratios(),
        how="left",
        on="prime_mover_code",
        validate="m:1",
    )

    # drop the prime mover code column
    gtn_conversions = gtn_conversions.drop(columns=["prime_mover_code"])

    return gtn_conversions


def calculate_subplant_nameplate_capacity(year):
    """Calculates the total nameplate capacity and primary prime mover for each CEMS subplant."""
    # load generator data
    pudl_out = load_data.initialize_pudl_out(year)
    gen_capacity = pudl_out.gens_eia860()[
        ["plant_id_eia", "generator_id", "prime_mover_code", "capacity_mw"]
    ]

    # add subplant ids to the generator data
    subplant_crosswalk = pd.read_csv(
        outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
        dtype=get_dtypes(),
    )[["plant_id_eia", "generator_id", "subplant_id"]].drop_duplicates()
    gen_capacity = gen_capacity.merge(
        subplant_crosswalk,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="1:1",
    )
    validation.test_for_missing_subplant_id(gen_capacity)
    subplant_capacity = (
        gen_capacity.groupby(["plant_id_eia", "subplant_id"])["capacity_mw"]
        .sum()
        .reset_index()
    )

    # identify the primary prime mover for each subplant based on capacity
    subplant_prime_mover = gen_capacity[
        gen_capacity.groupby(["plant_id_eia", "subplant_id"], dropna=False)[
            "capacity_mw"
        ].transform(max)
        == gen_capacity["capacity_mw"]
    ][["plant_id_eia", "subplant_id", "prime_mover_code"]].drop_duplicates(
        subset=["plant_id_eia", "subplant_id"], keep="first"
    )

    # add the prime mover information
    subplant_capacity = subplant_capacity.merge(
        subplant_prime_mover,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="1:1",
    )

    return subplant_capacity


def filter_gtn_conversion_factors(gtn_conversions):
    factors_to_use = gtn_conversions[
        [
            "plant_id_eia",
            "subplant_id",
            "report_date",
            "data_source",
            "gross_generation_mwh",
            "net_generation_mwh",
            "minimum_gross_generation_mwh",
            "maximum_gross_generation_mwh",
            "capacity_mw",
            "annual_subplant_shift_mw",
            "annual_plant_shift_mw",
            "annual_subplant_ratio",
            "annual_plant_ratio",
            "annual_fuel_ratio",
            "default_gtn_ratio",
        ]
    ]

    for shift_factor in ["annual_subplant_shift_mw", "annual_plant_shift_mw"]:
        # remove any shift factors that would lead net generation in any hour to be less than -50 MW
        factors_to_use.loc[
            factors_to_use["minimum_gross_generation_mwh"]
            + factors_to_use[shift_factor]
            <= -50,
            shift_factor,
        ] = np.NaN
        # remove any shift factors that would lead net generation in any hour to be greater than 150% of nameplate capacity
        factors_to_use.loc[
            factors_to_use["maximum_gross_generation_mwh"]
            + factors_to_use[shift_factor]
            > (factors_to_use["capacity_mw"] * 1.50),
            shift_factor,
        ] = np.NaN

    for scaling_factor in [
        "annual_subplant_ratio",
        "annual_plant_ratio",
        "annual_fuel_ratio",
    ]:
        # remove any factors that would scale net generation to less than 75% of gross generation
        # In general, the IQR of GTN ratios is between 0.75 and 1, with an upper bound around 1.25
        # remove any ratios that are negative to avoid flipping the shape of the profile
        factors_to_use.loc[
            factors_to_use[scaling_factor] < 0.75, scaling_factor
        ] = np.NaN
        # remove any factors that would cause the generation in any hour to exceed 125% of nameplate capacity
        factors_to_use.loc[
            (
                factors_to_use[scaling_factor]
                * factors_to_use["maximum_gross_generation_mwh"]
            )
            > (factors_to_use["capacity_mw"] * 1.25),
            scaling_factor,
        ] = np.NaN

    # All subplant-months at each plant should use the same method
    # if any annual_subplant_ratio are missing at a plant, revert to using annual_plant_ratio for the entire plant
    # if any annual_plant_ratio are missing for a plant, revert to using annual_subplant_shift
    # if any annual_subplant_shift are missing for a plant, revert to using annual_plant_shift
    method_hierarchy = [
        "annual_subplant_ratio",
        "annual_plant_ratio",
        "annual_subplant_shift_mw",
        "annual_plant_shift_mw",
    ]

    for method in method_hierarchy:
        # get a count of the number of non-na factor values and non-na net generation values for each plant
        incomplete_factors = (
            factors_to_use.groupby(["plant_id_eia", "data_source"], dropna=False)
            .count()[[method, "net_generation_mwh"]]
            .reset_index()
        )
        # only keep the plants where there are any missing factors for that plant
        incomplete_factors = incomplete_factors[
            (incomplete_factors[method] < incomplete_factors["net_generation_mwh"])
        ]

        # merge this list into factors_to_use, and set the factor to na for any plants that exist in the right df
        factors_to_use = factors_to_use.merge(
            incomplete_factors[["plant_id_eia", "data_source"]],
            how="outer",
            on=["plant_id_eia", "data_source"],
            indicator="incomplete_flag",
            validate="m:1",
        )
        factors_to_use.loc[factors_to_use["incomplete_flag"] == "both", method] = np.NaN
        factors_to_use = factors_to_use.drop(columns=["incomplete_flag"])

    return factors_to_use


def gross_to_net_regression(combined_gen_data, agg_level):
    """
    Regresses net generation on gross generation at either the plant or subplant level.
    """

    (
        gen_data_for_regression,
        plant_aggregation_columns,
    ) = aggregate_combined_gen_data_for_regression(combined_gen_data, agg_level)

    # calculate the hourly average generation values
    # gen_data["hours_in_month"] = gen_data["report_date"].dt.daysinmonth * 24
    gen_data_for_regression["gross_generation_mw"] = (
        gen_data_for_regression["gross_generation_mwh"]
        / gen_data_for_regression["hours_in_month"]
    )
    gen_data_for_regression["net_generation_mw"] = (
        gen_data_for_regression["net_generation_mwh"]
        / gen_data_for_regression["hours_in_month"]
    )

    # calculate the ratio for each plant and create a dataframe
    gtn_regression = (
        gen_data_for_regression.dropna()
        .groupby(plant_aggregation_columns, dropna=False)
        .apply(model_gross_to_net)
    )
    gtn_regression = gtn_regression.dropna()
    try:
        gtn_regression = pd.DataFrame(
            gtn_regression.tolist(),
            index=gtn_regression.index,
            columns=["slope", "intercept", "rsquared", "rsquared_adj", "observations"],
        ).reset_index()
    # if the dataframe is empty
    except AttributeError:
        gtn_regression = pd.DataFrame(
            [],
            index=gtn_regression.index,
            columns=["slope", "intercept", "rsquared", "rsquared_adj", "observations"],
        ).reset_index()
    """if not os.path.exists(outputs_folder(f"gross_to_net"):
        os.mkdir(outputs_folder(f"gross_to_net")

    gtn_regression.to_csv(
        outputs_folder(f"gross_to_net/{agg_level}_gross_to_net_regression.csv",
        index=False,
    )"""

    return gtn_regression


def aggregate_combined_gen_data_for_regression(combined_gen_data, agg_level):
    if agg_level == "plant":
        plant_aggregation_columns = ["plant_id_eia"]
    elif agg_level == "subplant":
        plant_aggregation_columns = ["plant_id_eia", "subplant_id"]
    else:
        raise ValueError("agg_level must be either 'plant' or 'subplant'")

    groupby_columns = plant_aggregation_columns + ["report_date"]

    # drop columns with missing data
    gen_data_for_regression = combined_gen_data.dropna(
        subset=["gross_generation_mwh", "net_generation_mwh"]
    )

    gen_data_for_regression = (
        gen_data_for_regression.groupby(groupby_columns, dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
        ]
        .sum()
        .reset_index()
    )

    return gen_data_for_regression, plant_aggregation_columns


def model_gross_to_net(df):
    """
    Performs a linear regression model of monthly gross to net generation.

    Performs recursive outlier removal up to two times if the absolute value of
    the studentizes residual > 3

    Args:
        df: dataframe containing all values of gross and net generation that should be regressed
    Returns:
        various model parameters
    """

    def model_data(df):
        # get a linear model for the data points
        model = smf.ols("net_generation_mw ~ gross_generation_mw", data=df).fit()
        outliers = model.outlier_test()

        return model, outliers

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # if a model is not able to be created with the first iteration, skip this data
        try:
            # get a linear model for the data points
            model, outliers = model_data(df)

            # find and remove outliers recursively up to two times
            # the first time removes any obvious outliers
            # the second time removes any outliers that may have been masked by the first outliers
            if abs(outliers["student_resid"]).max() > 3:
                # remove any outlier values
                df = df[
                    ~df.index.isin(outliers[abs(outliers["student_resid"]) > 3].index)
                ]

                # get a linear model of the corrected data
                try:
                    model, outliers = model_data(df)
                except ValueError:
                    pass

                # perform this removal one more time in case any outliers were masked by the first outlier(s)
                if abs(outliers["student_resid"]).max() > 3:
                    # remove any outlier values
                    df = df[
                        ~df.index.isin(
                            outliers[abs(outliers["student_resid"]) > 3].index
                        )
                    ]
                    # get a linear model for the data points
                    try:
                        model, outliers = model_data(df)
                    except ValueError:
                        pass

            # get outputs of final adjusted model
            slope = model.params[1]
            intercept = model.params[0]
            rsquared = model.rsquared
            rsquared_adj = model.rsquared_adj
            number_observations = model.nobs

            # if the intercept is positive, it may result in calculated net generation being larger than gross generation
            # thus, in this case, we will re-run the regression, forcing the intercept through zero
            """if intercept > 0:
                # get a linear model for the data points, forcing the intercept through zero
                model = smf.ols(
                    "net_generation_mw ~ gross_generation_mw - 1", data=df
                ).fit()
                # get outputs of final adjusted model
                slope = model.params[0]
                intercept = 0
                rsquared = model.rsquared
                rsquared_adj = model.rsquared_adj
                number_observations = model.nobs"""

            return slope, intercept, rsquared, rsquared_adj, number_observations

        except ValueError:
            pass


# Currently unused code for exploring gross to net conversions over multiple years
########################################################################################


def calculate_multiyear_gtn_factors(year, number_of_years):
    """This is the coordinating function for loading and calculating subplant IDs, GTN regressions, and GTN ratios."""
    start_year = year - (number_of_years - 1)
    end_year = year

    # load 5 years of monthly data from CEMS and EIA-923
    cems_monthly, gen_fuel_allocated = load_monthly_gross_and_net_generation(
        start_year, end_year
    )

    # add subplant ids to the data
    logger.info("Creating subplant IDs")
    cems_monthly, gen_fuel_allocated = data_cleaning.generate_subplant_ids(
        start_year, end_year, cems_monthly, gen_fuel_allocated
    )

    logger.info("Calculating Gross to Net regressions and ratios")
    # perform regression at subplant level
    gross_to_net_regression(
        gross_gen_data=cems_monthly,
        net_gen_data=gen_fuel_allocated,
        agg_level="subplant",
    )

    # perform regression at plant level
    gross_to_net_regression(
        gross_gen_data=cems_monthly, net_gen_data=gen_fuel_allocated, agg_level="plant"
    )

    # calculate monthly ratios at subplant level
    gross_to_net_ratio(
        gross_gen_data=cems_monthly,
        net_gen_data=gen_fuel_allocated,
        agg_level="subplant",
        year=year,
    )

    # calculate monthly ratios at plant level
    gross_to_net_ratio(
        gross_gen_data=cems_monthly,
        net_gen_data=gen_fuel_allocated,
        agg_level="plant",
        year=year,
    )


def load_monthly_gross_and_net_generation(start_year, end_year):
    # load cems data
    cems_monthly = load_data.load_cems_gross_generation(start_year, end_year)

    # load and clean EIA data
    # create pudl_out
    pudl_db = "sqlite:///" + downloads_folder("pudl/pudl.sqlite") 
    pudl_engine = sa.create_engine(pudl_db)
    pudl_out = pudl.output.pudltabl.PudlTabl(
        pudl_engine,
        freq="MS",
        start_date=f"{start_year}-01-01",
        end_date=f"{end_year}-12-31",
    )

    # allocate net generation and heat input to each generator-fuel grouping
    logger.info("Allocating EIA-923 generation data")
    gen_fuel_allocated = allocate_gen_fuel.allocate_gen_fuel_by_generator_energy_source(
        pudl_out, drop_interim_cols=True
    )
    # aggregate the allocated data to the generator level
    gen_fuel_allocated = allocate_gen_fuel.agg_by_generator(
        gen_fuel_allocated, sum_cols=["net_generation_mwh"]
    )

    return cems_monthly, gen_fuel_allocated


def gross_to_net_ratio(gross_gen_data, net_gen_data, agg_level, year):

    if agg_level == "plant":
        plant_aggregation_columns = ["plant_id_eia"]
    elif agg_level == "subplant":
        plant_aggregation_columns = ["plant_id_eia", "subplant_id"]
    else:
        raise ValueError("agg_level must be either 'plant' or 'subplant'")

    groupby_columns = plant_aggregation_columns + ["report_date"]

    gen_data = gross_gen_data.merge(
        net_gen_data,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="m:m",
    )

    # identify any rows where gross or net generation are missing
    incomplete_data = gen_data[
        gen_data[["gross_generation_mwh", "net_generation_mwh"]].isnull().any(axis=1)
    ]

    # load the activation and retirement dates into the data
    subplant_crosswalk = pd.read_csv(
        outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
        dtype=get_dtypes(),
    ).dropna(subset="emissions_unit_id_epa")
    incomplete_data = incomplete_data.merge(
        subplant_crosswalk,
        how="left",
        on=(["plant_id_eia", "subplant_id", "emissions_unit_id_epa", "generator_id"]),
        validate="m:1",
    ).drop(columns="plant_id_epa")

    # drop any of these rows where the retirement date is before the report date (only applies if net generation missing)
    incomplete_data = incomplete_data[
        ~(incomplete_data["retirement_date"] < incomplete_data["report_date"])
    ]

    # drop any of these rows where the report date is before the planned operating date
    incomplete_data = incomplete_data[
        ~(
            incomplete_data["report_date"]
            < incomplete_data["current_planned_operating_date"]
        )
    ]

    # get a list of unique subplant ids and report dates - this identifies where we have missing data we shouldn't calculate a ratio for
    incomplete_data = incomplete_data.drop_duplicates(subset=groupby_columns)[
        groupby_columns
    ]

    # only keep gen data that is not incomplete
    gtn_ratio = gen_data.merge(
        incomplete_data,
        how="outer",
        on=groupby_columns,
        indicator="source",
        validate="1:1",
    )
    gtn_ratio = gtn_ratio[gtn_ratio["source"] == "left_only"].drop(columns="source")

    # group data by aggregation columns
    gtn_ratio = gtn_ratio.groupby(groupby_columns, dropna=False).sum().reset_index()

    # calculate gross to net ratios for the remaining data
    gtn_ratio["gtn_ratio"] = (
        gtn_ratio["net_generation_mwh"] / gtn_ratio["gross_generation_mwh"]
    )

    gtn_ratio = gtn_ratio[groupby_columns + ["gtn_ratio"]]

    os.makedirs(outputs_folder("gross_to_net"), exist_ok=True)

    gtn_ratio.to_csv(
        outputs_folder(f"gross_to_net/{agg_level}_gross_to_net_ratio.csv"),
        index=False,
    )
