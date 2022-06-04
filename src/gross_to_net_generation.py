from audioop import cross
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import os
from pathlib import Path
import src.data_cleaning as data_cleaning
import sqlalchemy as sa
import warnings

import pudl.analysis.epa_crosswalk as epa_crosswalk
import pudl.analysis.allocate_net_gen as allocate_gen_fuel
import pudl.output.pudltabl


def identify_subplants_and_gtn_conversions(year, number_of_years):
    """This is the coordinating function for loading and calculating subplant IDs, GTN regressions, and GTN ratios."""
    start_year = year - (number_of_years - 1)
    end_year = year

    # TODO: move the following code to a separate function so that it does not hold these dataframes in memory after calculation

    # load 5 years of monthly data from CEMS and EIA-923
    cems_monthly, gen_fuel_allocated = load_monthly_gross_and_net_generation(
        start_year, end_year
    )

    # add subplant ids to the data
    print("Creating subplant IDs")
    cems_monthly, gen_fuel_allocated = generate_subplant_ids(
        start_year, end_year, cems_monthly, gen_fuel_allocated
    )

    print("Calculating Gross to Net regressions and ratios")
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
    )

    # calculate monthly ratios at plant level
    gross_to_net_ratio(
        gross_gen_data=cems_monthly, net_gen_data=gen_fuel_allocated, agg_level="plant"
    )


def load_monthly_gross_and_net_generation(start_year, end_year):
    # load cems data
    cems_monthly = load_cems_gross_generation(start_year, end_year)

    # load and clean EIA data
    # create pudl_out
    pudl_db = "sqlite:///../data/downloads/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)
    pudl_out = pudl.output.pudltabl.PudlTabl(
        pudl_engine,
        freq="MS",
        start_date=f"{start_year}-01-01",
        end_date=f"{end_year}-12-31",
    )

    # allocate net generation and heat input to each generator-fuel grouping
    print("Allocating EIA-923 generation data")
    gen_fuel_allocated = allocate_gen_fuel.allocate_gen_fuel_by_generator_energy_source(
        pudl_out, drop_interim_cols=True
    )
    # aggregate the allocated data to the generator level
    gen_fuel_allocated = allocate_gen_fuel.agg_by_generator(
        gen_fuel_allocated, sum_cols=["net_generation_mwh"]
    )

    return cems_monthly, gen_fuel_allocated


def load_cems_gross_generation(start_year, end_year):
    """Loads hourly CEMS gross generation data for multiple years."""
    cems_all = []

    for year in range(start_year, end_year + 1):
        print(f"loading {year} CEMS data")
        # specify the path to the CEMS data
        cems_path = f"../data/downloads/pudl/pudl_data/parquet/epacems/year={year}"

        # specify the columns to use from the CEMS database
        cems_columns = [
            "plant_id_eia",
            "unitid",
            "unit_id_epa",
            "datetime_utc",
            "operating_time_hours",
            "gross_load_mw",
        ]

        # load the CEMS data
        cems = pd.read_parquet(cems_path, columns=cems_columns)

        # only keep values when the plant was operating
        # this will help speed up calculations and allow us to add this data back later
        cems = cems[(cems["gross_load_mw"] > 0) | (cems["operating_time_hours"] > 0)]

        # rename cems plant_id_eia to plant_id_epa (PUDL simply renames the ORISPL_CODE column from the raw CEMS data as 'plant_id_eia' without actually crosswalking to the EIA id)
        # rename the heat content column to use the convention used in the EIA data
        cems = cems.rename(columns={"plant_id_eia": "plant_id_epa",})

        # if the unitid has any leading zeros, remove them
        cems["unitid"] = cems["unitid"].str.lstrip("0")

        # crosswalk the plant IDs and add a plant_id_eia column
        cems = data_cleaning.crosswalk_epa_eia_plant_ids(cems, year)

        # fill any missing values for operating time or steam load with zero
        cems["operating_time_hours"] = cems["operating_time_hours"].fillna(0)

        # calculate gross generation by multiplying gross_load_mw by operating_time_hours
        cems["gross_generation_mwh"] = (
            cems["gross_load_mw"] * cems["operating_time_hours"]
        )

        # add a report date
        cems = data_cleaning.add_report_date(cems)

        cems = cems[
            [
                "plant_id_eia",
                "unitid",
                "unit_id_epa",
                "report_date",
                "gross_generation_mwh",
            ]
        ]

        # group data by plant, unit, month
        cems = cems.groupby(
            ["plant_id_eia", "unitid", "unit_id_epa", "report_date"]
        ).sum()

        cems_all.append(cems)

    cems = pd.concat(cems_all, axis=0).reset_index()

    return cems


def manual_crosswalk_updates(crosswalk):
    # load manual matches
    crosswalk_manual = pd.read_csv("../data/manual/epa_eia_crosswalk_manual.csv").drop(
        columns=["notes"]
    )
    crosswalk_manual = crosswalk_manual.rename(
        columns={
            "plant_id_epa": "CAMD_PLANT_ID",
            "unitid": "CAMD_UNIT_ID",
            "plant_id_eia": "EIA_PLANT_ID",
            "generator_id": "EIA_GENERATOR_ID",
        }
    )

    # The EPA's crosswalk document incorrectly maps plant_id_epa 55248 to plant_id_eia 55248
    # the correct plant_id_eia is 2847
    crosswalk.loc[crosswalk["CAMD_PLANT_ID"] == 55248, "EIA_PLANT_ID"] = 2847

    # move missing crosswalk matches to a different dataframe
    unmatched = crosswalk.copy()[crosswalk["EIA_GENERATOR_ID"].isna()]
    crosswalk = crosswalk[~crosswalk["EIA_GENERATOR_ID"].isna()]

    # append the manual data to the crosswalk
    crosswalk = pd.concat([crosswalk, crosswalk_manual], axis=0)

    # filter the list of unmatched generators to those that were not in our manual list
    unmatched = unmatched.merge(
        crosswalk_manual,
        how="outer",
        on=["CAMD_PLANT_ID", "CAMD_UNIT_ID"],
        indicator="source",
        suffixes=(None, "_manual"),
    )
    unmatched = unmatched[unmatched["source"] == "left_only"].drop(
        columns=["EIA_PLANT_ID_manual", "EIA_GENERATOR_ID_manual", "source"]
    )

    # add these back to the crosswalk
    crosswalk = pd.concat([crosswalk, unmatched], axis=0)

    return crosswalk


def generate_subplant_ids(start_year, end_year, cems_monthly, gen_fuel_allocated):
    """
    Groups units and generators into unique subplant groups.

    This function consists of three primary parts:
    1. Identify a list of all unique plant-units that exist in the CEMS data
        for the years in question. This will be used to filter the crosswalk.
    2. Load the EPA-EIA crosswalk and filter it based on the units that exist
        in the CEMS data for the years in question
    3. Use graph analysis to identify distinct groupings of EPA units and EIA
        generators based on 1:1, 1:m, m:1, or m:m relationships.

    Returns:
        exports the subplant crosswalk to a csv file
        cems_monthly and gen_fuel_allocated with subplant_id added
    
    """

    ids = cems_monthly[["plant_id_eia", "unitid", "unit_id_epa"]].drop_duplicates()

    # load the crosswalk and filter it by the data that actually exists in cems
    crosswalk = pudl.output.epacems.epa_crosswalk()

    # update the crosswalk with manual matches
    crosswalk = manual_crosswalk_updates(crosswalk)

    # filter the crosswalk to drop any units that don't exist in CEMS
    filtered_crosswalk = epa_crosswalk.filter_crosswalk(crosswalk, ids)[
        [
            "plant_id_eia",
            "unitid",
            "CAMD_PLANT_ID",
            "CAMD_UNIT_ID",
            "CAMD_GENERATOR_ID",
            "EIA_PLANT_ID",
            "EIA_GENERATOR_ID",
        ]
    ]

    # change the plant id to an int
    filtered_crosswalk["EIA_PLANT_ID"] = filtered_crosswalk["EIA_PLANT_ID"].astype(int)

    # filter to generators that exist in the EIA data
    # get a list of unique generators in the EIA-923 data
    unique_eia_ids = gen_fuel_allocated[
        ["plant_id_eia", "generator_id"]
    ].drop_duplicates()
    filtered_crosswalk = unique_eia_ids.merge(
        filtered_crosswalk,
        left_on=["plant_id_eia", "generator_id"],
        right_on=["EIA_PLANT_ID", "EIA_GENERATOR_ID"],
        how="inner",
        suffixes=("_actual", None),
    ).drop(columns=["plant_id_eia_actual", "generator_id"])

    crosswalk_with_subplant_ids = epa_crosswalk.make_subplant_ids(filtered_crosswalk)
    # fix the column names
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.drop(
        columns=["plant_id_eia", "unitid", "CAMD_GENERATOR_ID"]
    )
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.rename(
        columns={
            "CAMD_PLANT_ID": "plant_id_epa",
            "EIA_PLANT_ID": "plant_id_eia",
            "CAMD_UNIT_ID": "unitid",
            "EIA_GENERATOR_ID": "generator_id",
        }
    )
    # change the eia plant id to int
    crosswalk_with_subplant_ids["plant_id_eia"] = crosswalk_with_subplant_ids[
        "plant_id_eia"
    ].astype(int)

    # change the order of the columns
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids[
        ["plant_id_epa", "unitid", "plant_id_eia", "generator_id", "subplant_id"]
    ]

    # add proposed operating dates and retirements to the subplant id crosswalk
    pudl_db = "sqlite:///../data/downloads/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)
    # get values starting with the year prior to teh start year so that we can get proposed operating dates for the start year (which are reported in year -1)
    pudl_out_status = pudl.output.pudltabl.PudlTabl(
        pudl_engine,
        freq="MS",
        start_date=f"{start_year-1}-01-01",
        end_date=f"{end_year}-12-31",
    )
    generator_status = pudl_out_status.gens_eia860().loc[
        :,
        [
            "plant_id_eia",
            "generator_id",
            "report_date",
            "operational_status",
            "current_planned_operating_date",
            "retirement_date",
        ],
    ]
    # only keep values that have a planned operating date or retirement date
    generator_status = generator_status[
        (~generator_status["current_planned_operating_date"].isna())
        | (~generator_status["retirement_date"].isna())
    ]
    # drop any duplicate entries
    generator_status = generator_status.sort_values(
        by=["plant_id_eia", "generator_id", "report_date"]
    ).drop_duplicates(
        subset=[
            "plant_id_eia",
            "generator_id",
            "current_planned_operating_date",
            "retirement_date",
        ],
        keep="last",
    )
    # for any generators that have different retirement or planned dates reported in different years, keep the most recent value
    generator_status = generator_status.sort_values(
        by=["plant_id_eia", "generator_id", "report_date"]
    ).drop_duplicates(subset=["plant_id_eia", "generator_id"], keep="last")

    # merge the dates into the crosswalk
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.merge(
        generator_status[
            [
                "plant_id_eia",
                "generator_id",
                "current_planned_operating_date",
                "retirement_date",
            ]
        ],
        how="left",
        on=["plant_id_eia", "generator_id"],
    )

    if not os.path.exists("../data/outputs/subplant_crosswalk"):
        os.mkdir("../data/outputs/subplant_crosswalk")

    # export the crosswalk to csv
    crosswalk_with_subplant_ids.to_csv(
        "../data/outputs/subplant_crosswalk.csv", index=False
    )

    # merge the subplant ids into each dataframe
    gen_fuel_allocated = gen_fuel_allocated.merge(
        crosswalk_with_subplant_ids[["plant_id_eia", "generator_id", "subplant_id"]],
        how="left",
        on=["plant_id_eia", "generator_id"],
    )
    cems_monthly = cems_monthly.merge(
        crosswalk_with_subplant_ids[["plant_id_eia", "unitid", "subplant_id"]],
        how="left",
        on=["plant_id_eia", "unitid"],
    )

    return cems_monthly, gen_fuel_allocated


def combine_gross_and_net_generation_data(gross_gen_data, net_gen_data, agg_level):
    if agg_level == "plant":
        plant_aggregation_columns = ["plant_id_eia"]
    elif agg_level == "subplant":
        plant_aggregation_columns = ["plant_id_eia", "subplant_id"]
    else:
        raise ValueError("agg_level must be either 'plant' or 'subplant'")

    groupby_columns = plant_aggregation_columns + ["report_date"]

    # aggregate the data and merge it together
    net_gen = (
        net_gen_data.groupby(groupby_columns)
        .sum(min_count=1)["net_generation_mwh"]
        .reset_index()
    )
    gross_gen = (
        gross_gen_data.groupby(groupby_columns)
        .sum()["gross_generation_mwh"]
        .reset_index()
    )
    gen_data = gross_gen.merge(net_gen, how="outer", on=groupby_columns)

    return gen_data, plant_aggregation_columns


def gross_to_net_regression(gross_gen_data, net_gen_data, agg_level):
    """
    Regresses net generation on gross generation at either the plant or subplant level.
    """

    gen_data, plant_aggregation_columns = combine_gross_and_net_generation_data(
        gross_gen_data, net_gen_data, agg_level
    )

    # calculate the hourly average generation values
    gen_data["hours_in_month"] = gen_data["report_date"].dt.daysinmonth * 24
    gen_data["gross_generation_mw"] = (
        gen_data["gross_generation_mwh"] / gen_data["hours_in_month"]
    )
    gen_data["net_generation_mw"] = (
        gen_data["net_generation_mwh"] / gen_data["hours_in_month"]
    )

    # calculate the ratio for each plant and create a dataframe
    gtn_regression = (
        gen_data.dropna().groupby(plant_aggregation_columns).apply(model_gross_to_net)
    )
    gtn_regression = pd.DataFrame(
        gtn_regression.tolist(),
        index=gtn_regression.index,
        columns=["slope", "intercept", "rsquared", "rsquared_adj", "observations"],
    ).reset_index()

    if not os.path.exists("../data/outputs/gross_to_net"):
        os.mkdir("../data/outputs/gross_to_net")

    gtn_regression.to_csv(
        f"../data/outputs/gross_to_net/{agg_level}_gross_to_net_regression.csv",
        index=False,
    )


def gross_to_net_ratio(gross_gen_data, net_gen_data, agg_level):

    if agg_level == "plant":
        plant_aggregation_columns = ["plant_id_eia"]
    elif agg_level == "subplant":
        plant_aggregation_columns = ["plant_id_eia", "subplant_id"]
    else:
        raise ValueError("agg_level must be either 'plant' or 'subplant'")

    groupby_columns = plant_aggregation_columns + ["report_date"]

    gen_data = gross_gen_data.merge(
        net_gen_data, how="outer", on=["plant_id_eia", "subplant_id", "report_date"]
    )

    # identify any rows where gross or net generation are missing
    incomplete_data = gen_data[
        gen_data[["gross_generation_mwh", "net_generation_mwh"]].isnull().any(axis=1)
    ]

    # load the activation and retirement dates into the data
    subplant_crosswalk = pd.read_csv(
        f"../data/outputs/subplant_crosswalk.csv"
    )
    incomplete_data = incomplete_data.merge(
        subplant_crosswalk,
        how="left",
        on=(["plant_id_eia", "subplant_id", "unitid", "generator_id"]),
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
        incomplete_data, how="outer", on=groupby_columns, indicator="source"
    )
    gtn_ratio = gtn_ratio[gtn_ratio["source"] == "left_only"].drop(columns="source")

    # group data by aggregation columns
    gtn_ratio = gtn_ratio.groupby(groupby_columns).sum().reset_index()

    # calculate gross to net ratios for the remaining data
    gtn_ratio["gtn_ratio"] = (
        gtn_ratio["net_generation_mwh"] / gtn_ratio["gross_generation_mwh"]
    )

    gtn_ratio = gtn_ratio[groupby_columns + ["gtn_ratio"]]

    if not os.path.exists("../data/outputs/gross_to_net"):
        os.mkdir("../data/outputs/gross_to_net")

    gtn_ratio.to_csv(
        f"../data/outputs/gross_to_net/{agg_level}_gross_to_net_ratio.csv", index=False
    )


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
            if intercept > 0:
                # get a linear model for the data points, forcing the intercept through zero
                model = smf.ols(
                    "net_generation_mw ~ gross_generation_mw - 1", data=df
                ).fit()
                # get outputs of final adjusted model
                slope = model.params[0]
                intercept = 0
                rsquared = model.rsquared
                rsquared_adj = model.rsquared_adj
                number_observations = model.nobs

            return slope, intercept, rsquared, rsquared_adj, number_observations

        except ValueError:
            pass

