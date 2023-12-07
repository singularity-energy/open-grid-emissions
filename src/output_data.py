import math
import pandas as pd
import numpy as np
import shutil
import os
import load_data
import column_checks
import validation
from filepaths import outputs_folder, results_folder, data_folder
from logging_util import get_logger

logger = get_logger(__name__)


GENERATED_EMISSION_RATE_COLS = [
    "generated_co2_rate_lb_per_mwh_for_electricity",
    "generated_ch4_rate_lb_per_mwh_for_electricity",
    "generated_n2o_rate_lb_per_mwh_for_electricity",
    "generated_co2e_rate_lb_per_mwh_for_electricity",
    "generated_nox_rate_lb_per_mwh_for_electricity",
    "generated_so2_rate_lb_per_mwh_for_electricity",
    "generated_co2_rate_lb_per_mwh_for_electricity_adjusted",
    "generated_ch4_rate_lb_per_mwh_for_electricity_adjusted",
    "generated_n2o_rate_lb_per_mwh_for_electricity_adjusted",
    "generated_co2e_rate_lb_per_mwh_for_electricity_adjusted",
    "generated_nox_rate_lb_per_mwh_for_electricity_adjusted",
    "generated_so2_rate_lb_per_mwh_for_electricity_adjusted",
]

CONSUMED_EMISSION_RATE_COLS = [
    "consumed_co2_rate_lb_per_mwh_for_electricity",
    "consumed_ch4_rate_lb_per_mwh_for_electricity",
    "consumed_n2o_rate_lb_per_mwh_for_electricity",
    "consumed_co2e_rate_lb_per_mwh_for_electricity",
    "consumed_nox_rate_lb_per_mwh_for_electricity",
    "consumed_so2_rate_lb_per_mwh_for_electricity",
    "consumed_co2_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_ch4_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_n2o_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_co2e_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_nox_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_so2_rate_lb_per_mwh_for_electricity_adjusted",
]

UNIT_CONVERSIONS = {"lb": ("kg", 0.453592), "mmbtu": ("GJ", 1.055056)}

TIME_RESOLUTIONS = {"hourly": "H", "monthly": "M", "annual": "A"}


def prepare_files_for_upload(years):
    """
    Zips files in preparation for upload to cloud storage and Zenodo.

    This should only be run when releasing a new minor or major version of the repo.
    """

    for year in years:
        zip_results_for_s3(year)
        zip_data_for_zenodo(year)


def zip_results_for_s3(year):
    """
    Zips results directories that contain more than a single file for hosting on an Amazon S3 bucket.
    """
    os.makedirs(data_folder("s3_upload"), exist_ok=True)
    for data_type in ["power_sector_data", "carbon_accounting", "plant_data"]:
        for aggregation in ["hourly", "monthly", "annual"]:
            for unit in ["metric_units", "us_units"]:
                if (
                    (data_type == "plant_data")
                    & (aggregation == "hourly")
                    & (unit == "metric_units")
                ):
                    # skip the metric hourly plant data since we do not create those outputs
                    pass
                else:
                    logger.info(
                        f"zipping {year}_{data_type}_{aggregation}_{unit} for s3"
                    )
                    folder = (
                        f"{results_folder()}/{year}/{data_type}/{aggregation}/{unit}"
                    )
                    shutil.make_archive(
                        f"{data_folder()}/s3_upload/{year}_{data_type}_{aggregation}_{unit}",
                        "zip",
                        root_dir=folder,
                        # base_dir="",
                    )
    # move and rename the plant attributes files
    shutil.copy(
        f"{results_folder()}/{year}/plant_data/plant_static_attributes.csv",
        f"{data_folder()}/s3_upload/plant_static_attributes_{year}.csv",
    )
    # archive the data quality metrics
    shutil.make_archive(
        f"{data_folder()}/s3_upload/{year}_data_quality_metrics",
        "zip",
        root_dir=f"{results_folder()}/{year}/data_quality_metrics",
        # base_dir="",
    )


def zip_data_for_zenodo(year):
    """
    Zips each of the four data directories for archiving on Zenodo.
    """
    os.makedirs(data_folder("zenodo"), exist_ok=True)
    for directory in ["outputs", "results"]:
        logger.info(f"zipping {directory}_{year} for zenodo")
        shutil.make_archive(
            data_folder(f"zenodo/{directory}_{year}"),
            "zip",
            root_dir=data_folder(f"{directory}/{year}"),
            # base_dir="",
        )


def output_intermediate_data(df, file_name, path_prefix, year, skip_outputs):
    column_checks.check_columns(df, file_name)
    if not skip_outputs:
        logger.info(f"Exporting {file_name} to data/outputs")
        df.to_csv(outputs_folder(f"{path_prefix}{file_name}_{year}.csv"), index=False)


def output_to_results(
    df, file_name, subfolder, path_prefix, skip_outputs, include_metric=True
):
    # Always check columns that should not be negative.
    small = "small" in path_prefix
    logger.info(f"Exporting {file_name} to data/results/{path_prefix}{subfolder}")

    if include_metric:
        metric = convert_results(df)
        metric = round_table(metric)
    df = round_table(df)

    # Check for negatives after rounding
    validation.test_for_negative_values(df, small)
    # check that there are no missing values
    validation.test_for_missing_values(df, small)

    if not skip_outputs:
        df.to_csv(
            results_folder(f"{path_prefix}{subfolder}us_units/{file_name}.csv"),
            index=False,
        )
        if include_metric:
            metric.to_csv(
                results_folder(f"{path_prefix}{subfolder}metric_units/{file_name}.csv"),
                index=False,
            )


def output_data_quality_metrics(df, file_name, path_prefix, skip_outputs):
    if not skip_outputs:
        logger.info(
            f"Exporting {file_name} to data/results/{path_prefix}data_quality_metrics"
        )

        # TODO: Add column checks

        df.to_csv(
            results_folder(f"{path_prefix}data_quality_metrics/{file_name}.csv"),
            index=False,
        )


def output_plant_data(df, path_prefix, resolution, skip_outputs, plant_attributes):
    """
    Helper function for plant-level output.
    Output for each time granularity, and output separately for real and shaped plants
    `df` contains all plant-level data, both CEMS and synthetic.

    Note: plant-level does not include rates, so all aggregation is summation.
    """
    if not skip_outputs:
        if resolution == "hourly":
            # output hourly data
            validation.validate_unique_datetimes(
                df, "individual_plant_data", ["plant_id_eia"]
            )
            validation.check_for_complete_timeseries(
                df, "individual_plant_data", ["plant_id_eia"], "year"
            )
            # Separately save real and aggregate plants
            output_to_results(
                df[df.plant_id_eia > 900000],
                "shaped_fleet_data",
                "plant_data/hourly/",
                path_prefix,
                skip_outputs,
            )
            output_to_results(
                df[df.plant_id_eia < 900000],
                "individual_plant_data",
                "plant_data/hourly/",
                path_prefix,
                skip_outputs,
            )

        elif resolution == "monthly":
            # output monthly data
            output_to_results(
                df,
                "plant_data",
                "plant_data/monthly/",
                path_prefix,
                skip_outputs,
            )
        elif resolution == "annual":
            # output annual data
            df = (
                df.groupby(["plant_id_eia"], dropna=False)
                .sum(numeric_only=True)
                .reset_index()
            )
            # check for anomalous looking co2 rates
            validation.check_for_anomalous_co2_factors(
                df, plant_attributes, min_threshold=10, max_threshold=15000
            )
            # Separately save real and aggregate plants
            output_to_results(
                df,
                "plant_data",
                "plant_data/annual/",
                path_prefix,
                skip_outputs,
            )


def convert_results(df):
    """
    Take df in US units (used throughout pipeline).
    Return a df with metric units.

    ASSUMPTIONS:
        * Columns to convert have names of form
            `co2_lb_per_mwh_produced` (mass),
            `co2_lb_per_mwh_produced_for_electricity` (rate),
            `fuel_consumed_mmbtu` (mass)
          meaning that unit to convert is ALWAYS in numerator
    """
    converted = df.copy(deep=True)
    for column in converted.columns:
        unit = ""
        for u in UNIT_CONVERSIONS.keys():  # What to convert?
            if u in column.split("_"):
                unit = u
                break
        if unit == "":
            continue  # nothing to convert, next column
        new_unit, factor = UNIT_CONVERSIONS[unit]
        new_col = column.replace(unit, new_unit)
        converted.rename(columns={column: new_col}, inplace=True)
        converted.loc[:, new_col] = converted.loc[:, new_col] * factor
    return converted


def write_generated_averages(ba_fuel_data, year, path_prefix, skip_outputs):
    if not skip_outputs:
        avg_fuel_type_production = (
            ba_fuel_data.groupby(["fuel_category"]).sum(numeric_only=True).reset_index()
        )
        # Add row for total before taking rates
        total = avg_fuel_type_production.mean(numeric_only=True).to_frame().T
        total.loc[0, "fuel_category"] = "total"
        avg_fuel_type_production = pd.concat([avg_fuel_type_production, total], axis=0)

        # Find rates
        for emission_type in ["_for_electricity", "_for_electricity_adjusted"]:
            for emission in ["co2", "ch4", "n2o", "co2e", "nox", "so2"]:
                avg_fuel_type_production[
                    f"generated_{emission}_rate_lb_per_mwh{emission_type}"
                ] = (
                    (
                        avg_fuel_type_production[f"{emission}_mass_lb{emission_type}"]
                        / avg_fuel_type_production["net_generation_mwh"]
                    )
                    .replace(np.inf, np.NaN)
                    .replace(-np.inf, np.NaN)
                    .fillna(0)
                )
        output_intermediate_data(
            avg_fuel_type_production,
            "annual_generation_averages_by_fuel",
            path_prefix,
            year,
            skip_outputs,
        )


def write_plant_metadata(
    plant_static_attributes,
    eia923_allocated,
    cems,
    partial_cems_subplant,
    partial_cems_plant,
    shaped_eia_data,
    path_prefix,
    skip_outputs,
):
    """
    Outputs metadata for each subplant-month.

    Include rows for subplants aggregated to a synthetic plant,
    so users can see when a plant's subplants are split across plant-level and synthetic hourly data files
    """

    KEY_COLUMNS = [
        "plant_id_eia",
        "subplant_id",
        "report_date",
    ]

    METADATA_COLUMNS = [
        "data_source",
        "hourly_profile_source",
        "net_generation_method",
    ]

    if not skip_outputs:
        # From monthly EIA data, we want only the EIA-only subplants -- these are the ones that got shaped
        eia_only_subplants = eia923_allocated[
            (eia923_allocated["hourly_data_source"] == "eia")
            & ~(eia923_allocated["fuel_consumed_mmbtu"].isna())
        ].copy()

        # identify the source
        cems["data_source"] = "CEMS"
        partial_cems_subplant["data_source"] = "EIA"
        partial_cems_plant["data_source"] = "EIA"
        shaped_eia_data["data_source"] = "EIA"
        eia_only_subplants["data_source"] = "EIA"

        # identify net generation method
        cems = cems.rename(columns={"gtn_method": "net_generation_method"})
        shaped_eia_data["net_generation_method"] = shaped_eia_data["profile_method"]
        eia_only_subplants["net_generation_method"] = "<See shaped plant ID>"
        partial_cems_subplant["net_generation_method"] = "scaled_partial_cems_subplant"
        partial_cems_plant["net_generation_method"] = "shaped_from_partial_cems_plant"

        # identify hourly profile method
        cems["hourly_profile_source"] = "CEMS"
        partial_cems_subplant["hourly_profile_source"] = "partial CEMS subplant"
        partial_cems_plant["hourly_profile_source"] = "partial CEMS plant"
        shaped_eia_data = shaped_eia_data.rename(
            columns={"profile_method": "hourly_profile_source"}
        )
        eia_only_subplants["hourly_profile_source"] = "<See shaped plant ID>"

        # only keep one metadata row per plant/subplant-month
        cems_meta = cems.copy()[KEY_COLUMNS + METADATA_COLUMNS].drop_duplicates(
            subset=KEY_COLUMNS
        )
        partial_cems_subplant_meta = partial_cems_subplant.copy()[
            KEY_COLUMNS + METADATA_COLUMNS
        ].drop_duplicates(subset=KEY_COLUMNS)
        partial_cems_plant_meta = partial_cems_plant.copy()[
            KEY_COLUMNS + METADATA_COLUMNS
        ].drop_duplicates(subset=KEY_COLUMNS)
        shaped_eia_data_meta = shaped_eia_data.copy()[
            ["plant_id_eia", "report_date"] + METADATA_COLUMNS
        ].drop_duplicates(subset=["plant_id_eia", "report_date"])
        monthly_eia_meta = eia_only_subplants.copy()[
            ["plant_id_eia", "report_date"] + METADATA_COLUMNS
        ].drop_duplicates(subset=["plant_id_eia", "report_date"])

        # For monthly only: specify which synthetic plant we were aggregated to
        monthly_eia_meta = monthly_eia_meta.merge(
            plant_static_attributes[["plant_id_eia", "shaped_plant_id"]],
            how="left",
            on="plant_id_eia",
            validate="m:1",  # There can be multiple subplants in monthly EIA for each plant in static attributes
        )

        # concat the metadata into a one file and export
        metadata = pd.concat(
            [
                cems_meta,
                partial_cems_subplant_meta,
                partial_cems_plant_meta,
                shaped_eia_data_meta,
                monthly_eia_meta,
            ],
            axis=0,
        )

        column_checks.check_columns(metadata, "plant_metadata")

        metadata.to_csv(
            results_folder(f"{path_prefix}plant_data/plant_metadata.csv"), index=False
        )

        # drop the metadata columns from each dataframe
        cems = cems.drop(columns=METADATA_COLUMNS)
        partial_cems_subplant = partial_cems_subplant.drop(columns=METADATA_COLUMNS)
        partial_cems_plant = partial_cems_plant.drop(columns=METADATA_COLUMNS)
        shaped_eia_data = shaped_eia_data.drop(columns=METADATA_COLUMNS)


def round_table(table):
    """
    Round each numeric column.
    All values in a column have the same rounding.
    Rounding for each col is based on the median non-zero value: if < 1, sigfigs = 3, else 2 decimal places
    """
    decimals = {}
    # Iterate through numeric columns
    for c in table.select_dtypes(include=np.number).columns:
        # Non-zero minimum
        val = table.loc[table[c] > 0, c].median()
        if pd.isna(val):  # if val is NaN, then this col has only NaN or only 0 values
            decimals[c] = 4
        # >1 gets 2 decimals
        elif val > 1:
            decimals[c] = 2
        else:
            try:
                decimals[c] = abs(math.floor(math.log10(val))) + 2
            # Always 3 sigfigs (for median)
            except ValueError:
                logger.error(val)
                raise Exception
    return table.round(decimals)


def write_power_sector_results(ba_fuel_data, path_prefix, skip_outputs):
    """
    Helper function to write combined data by BA
    """

    data_columns = [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "co2e_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "co2e_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "co2_mass_lb_for_electricity_adjusted",
        "ch4_mass_lb_for_electricity_adjusted",
        "n2o_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
    ]

    if not skip_outputs:
        for ba in list(ba_fuel_data.ba_code.unique()):
            if type(ba) is not str:
                logger.warning(
                    f"not aggregating {sum(ba_fuel_data.ba_code.isna())} plants with numeric BA {ba}"
                )
                continue

            # filter the data for a single BA
            ba_table = ba_fuel_data[ba_fuel_data["ba_code"] == ba].drop(
                columns="ba_code"
            )

            # convert the datetime_utc column back to a datetime
            ba_table["datetime_utc"] = pd.to_datetime(
                ba_table["datetime_utc"], utc=True
            ).astype("datetime64[s]")

            # calculate a total for the BA
            # grouping by datetime_utc and report_date will create some duplicate datetime_utc
            # values for certain bas where there are plants located in multiple timezones
            # the report date column is necessary for monthly aggregation, but we will have to
            # remove it and group values by datetime_utc for the hourly calculations
            ba_total = (
                ba_table.groupby(["datetime_utc", "report_date"], dropna=False)[
                    data_columns
                ]
                .sum()
                .reset_index()
            )
            ba_total["fuel_category"] = "total"

            # concat the totals to the fuel-specific totals
            ba_table = pd.concat([ba_table, ba_total], axis=0, ignore_index=True)

            # create a dataframe for the hourly values that groups duplicate datetime_utc values
            ba_table_hourly = ba_table.copy().drop(columns=["report_date"])
            ba_table_hourly = (
                ba_table_hourly.groupby(["fuel_category", "datetime_utc"])
                .sum()
                .reset_index()
            )

            def add_generated_emission_rate_columns(df):
                for emission_type in ["_for_electricity", "_for_electricity_adjusted"]:
                    for emission in ["co2", "ch4", "n2o", "co2e", "nox", "so2"]:
                        col_name = (
                            f"generated_{emission}_rate_lb_per_mwh{emission_type}"
                        )
                        df[col_name] = (
                            (
                                df[f"{emission}_mass_lb{emission_type}"]
                                / df["net_generation_mwh"]
                            )
                            .replace(np.inf, np.NaN)
                            .replace(-np.inf, np.NaN)
                        )
                        # where the rate is missing because of a divide by zero (i.e.
                        # net_generation_mwh is zero), replace the emission rate with
                        # zero. We want to keep all other NAs so that they get flagged
                        # by our validation checks since this indicates an unexpected
                        # issue
                        df.loc[df["net_generation_mwh"] == 0, col_name] = df.loc[
                            df["net_generation_mwh"] == 0, col_name
                        ].fillna(0)
                        # Set negative rates to zero, following eGRID methodology
                        df.loc[df[col_name] < 0, col_name] = 0
                return df

            # output the hourly data
            ba_table_hourly = add_generated_emission_rate_columns(ba_table_hourly)

            # create a local datetime column
            try:
                local_tz = load_data.ba_timezone(ba, "local")
                ba_table_hourly["datetime_local"] = ba_table_hourly[
                    "datetime_utc"
                ].dt.tz_convert(local_tz)
            except ValueError:
                ba_table_hourly["datetime_local"] = pd.NaT

            # re-order columns
            ba_table_hourly = ba_table_hourly[
                ["fuel_category", "datetime_local", "datetime_utc"]
                + data_columns
                + GENERATED_EMISSION_RATE_COLS
            ]

            validation.validate_unique_datetimes(
                df=ba_table_hourly,
                df_name="power sector hourly ba table",
                keys=["fuel_category"],
            )
            validation.check_for_complete_timeseries(
                ba_table_hourly,
                "power sector hourly ba table",
                ["fuel_category"],
                "year",
            )

            # export to a csv
            output_to_results(
                ba_table_hourly,
                ba,
                "power_sector_data/hourly/",
                path_prefix,
                skip_outputs,
            )

            # aggregate data to monthly
            ba_table_monthly = (
                ba_table.groupby(["fuel_category", "report_date"], dropna=False)
                .sum(numeric_only=True)
                .reset_index()
            )
            ba_table_monthly = add_generated_emission_rate_columns(ba_table_monthly)
            # re-order columns
            ba_table_monthly = ba_table_monthly[
                ["fuel_category", "report_date"]
                + data_columns
                + GENERATED_EMISSION_RATE_COLS
            ]
            output_to_results(
                ba_table_monthly,
                ba,
                "power_sector_data/monthly/",
                path_prefix,
                skip_outputs,
            )

            # aggregate data to annual
            ba_table_annual = (
                ba_table.groupby(["fuel_category"], dropna=False)
                .sum(numeric_only=True)
                .reset_index()
            )
            ba_table_annual = add_generated_emission_rate_columns(ba_table_annual)
            # re-order columns
            ba_table_annual = ba_table_annual[
                ["fuel_category"] + data_columns + GENERATED_EMISSION_RATE_COLS
            ]
            output_to_results(
                ba_table_annual,
                ba,
                "power_sector_data/annual/",
                path_prefix,
                skip_outputs,
            )
