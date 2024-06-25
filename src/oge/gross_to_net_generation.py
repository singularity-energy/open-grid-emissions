import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import warnings

# import other modules
import oge.load_data as load_data
import oge.helpers as helpers
import oge.validation as validation

from oge.data_cleaning import assign_fuel_type_to_cems
from oge.helpers import create_plant_ba_table, add_subplant_ids_to_df
from oge.logging_util import get_logger

logger = get_logger(__name__)


def convert_gross_to_net_generation(
    cems: pd.DataFrame,
    eia923_allocated: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Converts hourly gross generation in CEMS to hourly net generation.

    This function first calculates various types of gross to net conversion factors,
    filters them for quality, and them applies them according to a method hierarchy:
        1. Annual GTN ratio at the subplant level
        2. Annual GTN ratio at the plant level
        3. Annual GTN ratio at the national fleet (fuel category-prime mover) level
        4. Default GTN ratio published by EIA based on prime mover
        5. Assumed GTN ratio of 0.97

    Args:
        cems (pd.DataFrame): hourly CEMS data aggregated to the subplant level
        eia923_allocated (pd.DataFrame): Cleaned EIA-923 data at the subplant level
        primary_fuel_table (pd.DataFrame): Table indicating subplant primary fuel
        year (int): data year

    Returns:
        pd.DataFrame: cems df with an added column for net_generation_mwh and a column
            indicating the method used to calculate net generation
    """

    gtn_conversions = calculate_gross_to_net_conversion_factors(
        cems, eia923_allocated, primary_fuel_table, year
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
                "annual_subplant_ratio",
                "annual_plant_ratio",
                "annual_fleet_ratio",
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

    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "3_annual_fleet_ratio"
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["annual_fleet_ratio"]
    )

    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "4_default_eia_ratio"
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["default_gtn_ratio"]
    )

    # warn if there are any missing default gtn ratios for plants that would use them.
    missing_defaults = cems.loc[
        (cems["gtn_method"] == "4_default_eia_ratio")
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
            .merge(
                create_plant_ba_table(year)[["plant_id_eia", "ba_code"]],
                how="left",
                on="plant_id_eia",
                validate="m:1",
            )
            .to_string()
        )
    # if there is a missing default gtn ratio, fill with 0.97
    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "5_assumed_gtn_ratio"
    cems["default_gtn_ratio"] = cems["default_gtn_ratio"].fillna(0.97)
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["default_gtn_ratio"]
    )

    # drop intermediate columns
    cems = cems.drop(
        columns=[
            "data_source",
            "annual_subplant_ratio",
            "annual_plant_ratio",
            "annual_fleet_ratio",
            "default_gtn_ratio",
        ]
    )

    validation.validate_gross_to_net_conversion(cems, eia923_allocated, year)

    return cems, gtn_conversions


def calculate_gross_to_net_conversion_factors(
    cems: pd.DataFrame,
    eia923_allocated: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Calculates gross to net ratios based on gross generation data reported in CEMS
    and net generation data reported in EIA-923.

    Calculates ratios for specific subplants, plants, and fleets (fuel-PM).
    When calculating ratios, we ensure to only keep data where there is data both for
    CEMS and EIA.

    Args:
        cems (pd.DataFrame): hourly CEMS data aggregated to the subplant level
        eia923_allocated (pd.DataFrame): Cleaned EIA-923 data at the subplant level
        primary_fuel_table (pd.DataFrame): Table indicating subplant primary fuel
        year (int): data year

    Returns:
        pd.DataFrame: table containing monthly and annual subplant and plant ratios,
            and annual fleet ratios for each subplant.
    """

    # aggregate the hourly cems data by subplant
    gross_gen_data = cems[
        [
            "plant_id_eia",
            "subplant_id",
            "report_date",
            "datetime_utc",
            "gross_generation_mwh",
            "fuel_consumed_mmbtu",
        ]
    ].copy()
    # add energy source codes to teh data
    gross_gen_data = assign_fuel_type_to_cems(gross_gen_data, year, primary_fuel_table)
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
        .agg(
            {
                "datetime_utc": "count",
                "gross_generation_mwh": "sum",
                "energy_source_code": "first",
            }
        )
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
            net_gen_data,
            how="outer",
            on=["plant_id_eia", "subplant_id", "report_date"],
            indicator="data_source",
            validate="1:1",
        )
        .merge(
            subplant_capacity,
            how="left",
            on=["plant_id_eia", "subplant_id"],
            validate="m:1",
        )
    )
    combined_gen_data["data_source"] = combined_gen_data[
        "data_source"
    ].cat.rename_categories(
        {
            "left_only": "cems_only",
            "right_only": "eia_only",
        }
    )

    # validate the data
    validation.check_missing_or_zero_generation_matches(combined_gen_data, year)

    # calculate other groupings at the plant and annual levels
    annual_subplant_ratio = (
        combined_gen_data[combined_gen_data["data_source"] == "both"]
        .dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia", "subplant_id"], dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
        ]
        .sum()
        .reset_index()
    )
    monthly_plant_ratio = (
        combined_gen_data[combined_gen_data["data_source"] == "both"]
        .dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia", "report_date"], dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh"]
        ]
        .sum()
        .reset_index()
    )
    annual_plant_ratio = (
        combined_gen_data[combined_gen_data["data_source"] == "both"]
        .dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia"], dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
        ]
        .sum()
        .reset_index()
    )
    # calculate a fleet-level ratio
    # first assign fuel categories to the data
    annual_fleet_ratio = helpers.assign_fuel_category_to_esc(
        combined_gen_data[combined_gen_data["data_source"] == "both"].dropna(
            subset=["gross_generation_mwh", "net_generation_mwh"]
        )
    )
    # first group the data by subplant and remove any data that is anomalous
    annual_fleet_ratio = (
        annual_fleet_ratio.groupby(
            ["plant_id_eia", "subplant_id", "fuel_category", "prime_mover_code"],
            dropna=False,
        )[["gross_generation_mwh", "net_generation_mwh"]]
        .sum()
        .reset_index()
    )
    annual_fleet_ratio["annual_fleet_ratio"] = (
        annual_fleet_ratio["net_generation_mwh"]
        / annual_fleet_ratio["gross_generation_mwh"]
    ).replace([np.inf, -np.inf], np.nan)
    annual_fleet_ratio = annual_fleet_ratio[
        (annual_fleet_ratio["annual_fleet_ratio"] <= 1.25)
        & (annual_fleet_ratio["annual_fleet_ratio"] >= 0.75)
    ]

    # sum data by fuel category
    annual_fleet_ratio = (
        annual_fleet_ratio.groupby(["fuel_category", "prime_mover_code"], dropna=False)[
            ["gross_generation_mwh", "net_generation_mwh"]
        ]
        .sum()
        .reset_index()
    )

    annual_fleet_ratio["annual_fleet_ratio"] = (
        annual_fleet_ratio["net_generation_mwh"]
        / annual_fleet_ratio["gross_generation_mwh"]
    ).replace([np.inf, -np.inf], np.nan)

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
    validation.identify_anomalous_annual_plant_gtn_ratios(annual_plant_ratio, year)

    # drop the gross and net generation data from the dataframes at the other
    # aggregation levels
    annual_subplant_ratio = annual_subplant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
    )
    monthly_plant_ratio = monthly_plant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh"]
    )
    annual_plant_ratio = annual_plant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh"]
    )
    annual_fleet_ratio = annual_fleet_ratio.drop(
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
    # merge the plant primary fuel and the fuel ratios into the conversion table
    gtn_conversions = gtn_conversions.merge(
        primary_fuel_table[
            ["plant_id_eia", "subplant_id", "subplant_primary_fuel"]
        ].drop_duplicates(),
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )
    gtn_conversions = helpers.assign_fuel_category_to_esc(gtn_conversions)
    gtn_conversions = gtn_conversions.merge(
        annual_fleet_ratio,
        how="left",
        on=["fuel_category", "prime_mover_code"],
        validate="m:1",
    )

    # where gross or net generation data was missing in a month, change the monthly
    # ratios to missing
    gtn_conversions.loc[
        gtn_conversions[["gross_generation_mwh", "net_generation_mwh"]]
        .isna()
        .any(axis=1),
        ["monthly_subplant_ratio", "monthly_plant_ratio"],
    ] = np.NaN

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

    # drop intermediate columns
    gtn_conversions = gtn_conversions.drop(
        columns=[
            "prime_mover_code",
            "subplant_primary_fuel",
            "fuel_category",
            "fuel_category_eia930",
        ]
    )

    return gtn_conversions


def calculate_subplant_nameplate_capacity(year):
    """Calculates the total nameplate capacity and primary prime mover for each CEMS subplant."""
    # load generator data
    gen_capacity = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year,
        columns=[
            "plant_id_eia",
            "generator_id",
            "prime_mover_code",
            "capacity_mw",
            "operational_status_code",
        ],
    )

    # add subplant ids to the generator data
    gen_capacity = add_subplant_ids_to_df(
        gen_capacity,
        year,
        plant_part_to_map="generator_id",
        how_merge="inner",
        validate_merge="1:1",
    )
    subplant_capacity = (
        gen_capacity.groupby(["plant_id_eia", "subplant_id"])["capacity_mw"]
        .sum()
        .reset_index()
    )

    # identify the primary prime mover for each subplant based on capacity
    subplant_prime_mover = gen_capacity[
        gen_capacity.groupby(["plant_id_eia", "subplant_id"], dropna=False)[
            "capacity_mw"
        ].transform("max")
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


def filter_gtn_conversion_factors(gtn_conversions: pd.DataFrame) -> pd.DataFrame:
    """Filters the calculated GTN ratios to remove anomalous or incomplete factors.

    First, we remove any ratios that are less than 0.75 or greater than 1.25.
    We also want to ensure that at each plant, we either use all annual_subplant_ratio
    or all annual_plant_ratio so that the annual plant total net generation matches. We
    remove any annual_subplant_ratios if they are not available for all subplants at a
    plant.

    Args:
        gtn_conversions (pd.DataFrame): output of
            calculate_gross_to_net_conversion_factors()

    Returns:
        pd.DataFrame: filtered gtn_conversions
    """
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
            "annual_subplant_ratio",
            "annual_plant_ratio",
            "annual_fleet_ratio",
            "default_gtn_ratio",
        ]
    ]

    for scaling_factor in [
        "annual_subplant_ratio",
        "annual_plant_ratio",
        "annual_fleet_ratio",
    ]:
        # remove any factors that would scale net generation to less than 75% of gross
        # generation. In general, the IQR of GTN ratios is between 0.75 and 1, with an
        # upper bound around 1.25. Remove any ratios that are negative to avoid flipping
        # the shape of the profile
        factors_to_use.loc[factors_to_use[scaling_factor] < 0.75, scaling_factor] = (
            np.NaN
        )
        # remove any factors that would cause the generation in any hour to exceed 125%
        # of nameplate capacity
        factors_to_use.loc[
            (
                factors_to_use[scaling_factor]
                * factors_to_use["maximum_gross_generation_mwh"]
            )
            > (factors_to_use["capacity_mw"] * 1.25),
            scaling_factor,
        ] = np.NaN

    # All subplant-months at each plant should use the same method
    # if any annual_subplant_ratio are missing at a plant, revert to using
    # annual_plant_ratio for the entire plant
    method_hierarchy = [
        "annual_subplant_ratio",
        "annual_plant_ratio",
    ]

    for method in method_hierarchy:
        # get a count of the number of non-na factor values and non-na net generation
        # values for each plant
        incomplete_factors = (
            factors_to_use.groupby(
                ["plant_id_eia", "data_source"], dropna=False, observed=False
            )
            .count()[[method, "net_generation_mwh"]]
            .reset_index()
        )
        # only keep the plants where there are any missing factors for that plant
        incomplete_factors = incomplete_factors[
            (incomplete_factors[method] < incomplete_factors["net_generation_mwh"])
        ]

        # merge this list into factors_to_use, and set the factor to na for any plants
        # that exist in the right df
        factors_to_use = factors_to_use.merge(
            incomplete_factors[["plant_id_eia", "data_source"]],
            how="outer",
            on=["plant_id_eia", "data_source"],
            indicator="incomplete_flag",
            validate="m:1",
        )
        factors_to_use.loc[factors_to_use["incomplete_flag"] == "both", method] = np.NaN
        factors_to_use = factors_to_use.drop(
            columns=["incomplete_flag", "energy_source_code"]
        )

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
