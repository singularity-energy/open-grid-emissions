import math
import pandas as pd
import numpy as np
import shutil
import os

import oge.load_data as load_data
import oge.column_checks as column_checks
import oge.validation as validation
from oge.filepaths import outputs_folder, results_folder, data_folder
from oge.logging_util import get_logger
from oge.constants import ConversionFactors

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

UNIT_CONVERSIONS = {
    "lb": ("kg", ConversionFactors.lb_to_kg),
    "mmbtu": ("GJ", ConversionFactors.mmbtu_to_GJ),
}


def prepare_files_for_upload(years: list):
    """Zips files in preparation for upload to cloud storage and Zenodo.
    This should only be run when releasing a new minor or major version of the repo.

    Args:
        years (list): list of four-digit year indicating when the data were taken.
    """

    for year in years:
        zip_results_for_s3(year)
        zip_data_for_zenodo(year)


def zip_results_for_s3(year: int):
    """Zips results directories that contain more than a single file for hosting on an
    Amazon S3 bucket.


    Args:
        year (int): a four-digit year indicating when the data were taken.
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


def zip_data_for_zenodo(year: int):
    """Zips each of the three data directories for archiving on Zenodo.

    Args:
        year (int): a four-digit year indicating when the data were taken.
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


def output_intermediate_data(
    df: pd.DataFrame, file_name: str, path_prefix: str, year: int, skip_outputs: bool
):
    """Save data frame as ZIP into the outputs directory.

    Args:
        df (pd.DataFrame): data frame that will be saved.
        file_name (str): name of file without file extension.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        year (int): a four-digit year indicating when the data were taken.
        skip_outputs (bool): whether to save data or not.
    """
    column_checks.check_columns(df, file_name)
    if not skip_outputs:
        logger.info(f"Exporting {file_name} to data/outputs")
        df.to_csv(
            outputs_folder(f"{path_prefix}{file_name}_{year}.csv.zip"),
            index=False,
            compression="zip",
        )


def output_to_results(
    df: pd.DataFrame,
    year: int,
    file_name: str,
    subfolder: str,
    path_prefix: str,
    skip_outputs: bool,
    include_metric=True,
):
    """Save data franme as CSV into the results directory.

    Args:
        df (pd.DataFrame): data frame that will be saved.
        year (int):a four-digit year indicating when the data were taken.
        file_name (str): name of file without file extension.
        subfolder (str): name of directory following `path_prefix`.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        skip_outputs (bool): whether to save the data or not.
        include_metric (bool, optional): whether to create a file with values in metric
            units. Defaults to True.
    """
    # Always check columns that should not be negative.
    small = "small" in path_prefix
    logger.info(f"Exporting {file_name} to data/results/{path_prefix}{subfolder}")

    if include_metric:
        metric = convert_results(df)
        metric = round_table(metric)
    df = round_table(df)

    # Check for negatives after rounding
    validation.test_for_negative_values(df, year, small)
    # Check that there are no missing values
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


def output_data_quality_metrics(
    df: pd.DataFrame, file_name: str, path_prefix: str, skip_outputs: bool
):
    """Output data quality metrics.

    Args:
        df (pd.DataFrame): data frame that will be saved.
        file_name (str): name of file without file extension.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        skip_outputs (bool): whether to save data or not.
    """
    if not skip_outputs:
        logger.info(
            f"Exporting {file_name} to data/results/{path_prefix}data_quality_metrics"
        )

        # TODO: Add column checks

        df.to_csv(
            results_folder(f"{path_prefix}data_quality_metrics/{file_name}.csv"),
            index=False,
        )


def output_plant_data(
    df: pd.DataFrame,
    year: int,
    path_prefix: str,
    resolution: str,
    skip_outputs: bool,
    plant_attributes: pd.DataFrame,
):
    """Helper function for plant-level output. Different output will be produced for
    real and shaped plants.

    Note:
        plant-level does not include rates, so all aggregation is summation.

    Args:
        df (pd.DataFrame): plant-level data both CEMS and synthetic.
        year (int): a four-digit year indicating when data were taken.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        resolution (str): temporal resolution. Wither 'hourly', 'monthly' or 'annual'.
        skip_outputs (bool): whether to save data or not.
        plant_attributes (pd.DataFrame): the plant static attributes table.
    """
    if not skip_outputs:
        if resolution == "hourly":
            # output hourly data
            validation.validate_unique_datetimes(
                year, df, "individual_plant_data", ["plant_id_eia"]
            )
            validation.check_for_complete_hourly_timeseries(
                df, "individual_plant_data", ["plant_id_eia"], "year"
            )
            # Separately save real and aggregate plants
            output_to_results(
                df[df.plant_id_eia > 900000],
                year,
                "shaped_fleet_data",
                "plant_data/hourly/",
                path_prefix,
                skip_outputs,
            )
            output_to_results(
                df[df.plant_id_eia < 900000],
                year,
                "individual_plant_data",
                "plant_data/hourly/",
                path_prefix,
                skip_outputs,
            )
        elif resolution == "monthly":
            # output monthly data
            output_to_results(
                df,
                year,
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
                df, plant_attributes, year, min_threshold=10, max_threshold=15000
            )
            # Separately save real and aggregate plants
            output_to_results(
                df,
                year,
                "plant_data",
                "plant_data/annual/",
                path_prefix,
                skip_outputs,
            )


def convert_results(df: pd.DataFrame) -> pd.DataFrame:
    """Convert values in data frame from US units to metric units.

    Note:
        Columns to convert have names of form: 'co2_lb_per_mwh_produced' (mass),
        'co2_lb_per_mwh_produced_for_electricity' (rate) and 'uel_consumed_mmbtu'
        (mass) meaning that nit to convert is always in numerator.

    Args:
        df (pd.DataFrame): data frame in US units.

    Returns:
        pd.DataFrame: data frame in metric units.
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


def write_generated_averages(
    ba_fuel_data: pd.DataFrame, year: int, path_prefix: str, skip_outputs: bool
):
    """Outputs generated averaged emission.

    Args:
        ba_fuel_data (pd.DataFrame): plant data aggregated by BA and fuel type.
        year (int): a four-digit year indicating when the data were taken.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        skip_outputs (bool): whether to save data or not.
    """
    if not skip_outputs:
        fuel_type_production = (
            ba_fuel_data.groupby(["fuel_category"]).sum(numeric_only=True).reset_index()
        )

        # Add row for total
        total_production = fuel_type_production.sum(numeric_only=True).to_frame().T
        total_production.loc[0, "fuel_category"] = "total"
        production = pd.concat([fuel_type_production, total_production], axis=0)

        # Calculate rates
        for emission_type in ["_for_electricity", "_for_electricity_adjusted"]:
            for emission in ["co2", "ch4", "n2o", "co2e", "nox", "so2"]:
                production[f"generated_{emission}_rate_lb_per_mwh{emission_type}"] = (
                    (
                        production[f"{emission}_mass_lb{emission_type}"]
                        / production["net_generation_mwh"]
                    )
                    .replace(np.inf, np.NaN)
                    .replace(-np.inf, np.NaN)
                    .fillna(0)
                )

        output_intermediate_data(
            production,
            "annual_generation_averages_by_fuel",
            path_prefix,
            year,
            skip_outputs,
        )


def write_plant_metadata(
    plant_static_attributes: pd.DataFrame,
    eia923_allocated: pd.DataFrame,
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
    shaped_eia_data: pd.DataFrame,
    path_prefix: str,
    skip_outputs: bool,
    year: int,
):
    """Outputs plant metadata or each subplant-month. Includes rows for subplants
    aggregated to a synthetic plant, so users can see when a plant's subplants are
    split across plant-level and synthetic hourly data files

    Args:
        plant_static_attributes (pd.DataFrame): plant static attributes table.
        eia923_allocated (pd.DataFrame): allocated EIA-923 generation and fuel
            consumption data.
        cems (pd.DataFrame): CEMS data.
        partial_cems_subplant (pd.DataFrame): subplant data for which there is partial
            reporting in CEMS.
        partial_cems_plant (pd.DataFrame): subplant data for which there is partial
            plant data reporting in CEMS.
        shaped_eia_data (pd.DataFrame): hourly generation profile derived from
            monthly-level EIA data.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        skip_outputs (bool): whether to save data or not.
        year (int): a four-digit year indicating when data were taken.
    """

    KEY_COLUMNS = [
        "plant_id_eia",
        "subplant_id",
        "report_date",
    ]

    # create CEMS metadata
    # only keep one metadata row per plant/subplant-month
    cems_meta = cems.copy()[KEY_COLUMNS + ["gtn_method"]].drop_duplicates(
        subset=KEY_COLUMNS
    )
    cems_meta["data_source"] = "CEMS"
    cems_meta = cems_meta.rename(columns={"gtn_method": "net_generation_method"})
    cems_meta["hourly_profile_source"] = "CEMS"

    # create partial cems subplant metadata
    partial_cems_subplant_meta = partial_cems_subplant.copy()[
        KEY_COLUMNS
    ].drop_duplicates(subset=KEY_COLUMNS)
    partial_cems_subplant_meta["data_source"] = "EIA"
    partial_cems_subplant_meta["net_generation_method"] = "scaled_partial_cems_subplant"
    partial_cems_subplant_meta["hourly_profile_source"] = "partial CEMS subplant"

    # create partial cems plant metadata
    partial_cems_plant_meta = partial_cems_plant.copy()[KEY_COLUMNS].drop_duplicates(
        subset=KEY_COLUMNS
    )
    partial_cems_plant_meta["data_source"] = "EIA"
    partial_cems_plant_meta["net_generation_method"] = "shaped_from_partial_cems_plant"
    partial_cems_plant_meta["hourly_profile_source"] = "partial CEMS plant"

    # create EIA-only metadata
    # From monthly EIA data, we want only the EIA-only subplants
    # these are the ones that got shaped
    monthly_eia_meta = (
        eia923_allocated.copy()
        .loc[
            (eia923_allocated["hourly_data_source"] == "eia")
            & ~(eia923_allocated["fuel_consumed_mmbtu"].isna()),
            ["plant_id_eia", "subplant_id", "report_date"],
        ]
        .drop_duplicates(subset=["plant_id_eia", "subplant_id", "report_date"])
    )

    # For monthly only: specify which synthetic plant we were aggregated to
    monthly_eia_meta = monthly_eia_meta.merge(
        plant_static_attributes[["plant_id_eia", "shaped_plant_id"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",  # There can be multiple subplants for each plant
    )

    monthly_eia_meta["data_source"] = "EIA"

    # merge in the net generation method and hourly profile from the shaped metadata
    shaped_eia_data_meta = shaped_eia_data.copy()[
        ["plant_id_eia", "report_date", "profile_method"]
    ].drop_duplicates(subset=["plant_id_eia", "report_date"])
    shaped_eia_data_meta = shaped_eia_data_meta.rename(
        columns={"plant_id_eia": "shaped_plant_id"}
    )
    shaped_eia_data_meta["net_generation_method"] = shaped_eia_data_meta[
        "profile_method"
    ]
    shaped_eia_data_meta = shaped_eia_data_meta.rename(
        columns={"profile_method": "hourly_profile_source"}
    )
    monthly_eia_meta = monthly_eia_meta.merge(
        shaped_eia_data_meta,
        how="left",
        on=["shaped_plant_id", "report_date"],
        validate="m:1",  # There can be multiple subplants for each plant
    )

    monthly_eia_meta = monthly_eia_meta.drop(columns=["shaped_plant_id"])

    # concat the metadata into a one file and export
    metadata = pd.concat(
        [
            cems_meta,
            partial_cems_subplant_meta,
            partial_cems_plant_meta,
            monthly_eia_meta,
        ],
        axis=0,
    ).sort_values(by=["plant_id_eia", "subplant_id", "report_date"], ascending=True)

    validation.check_for_complete_monthly_timeseries(
        df=metadata,
        df_name="plant_metadata",
        keys=["plant_id_eia", "subplant_id"],
        columns_to_check=["hourly_profile_source"],
        year=year,
    )

    column_checks.check_columns(metadata, "plant_metadata")

    if not skip_outputs:
        metadata.to_csv(
            results_folder(f"{path_prefix}plant_data/plant_metadata.csv"), index=False
        )


def round_table(table: pd.DataFrame) -> pd.DataFrame:
    """Round each numeric columns. All values in a column have the same rounding.
    Rounding for each column is based on the median non-zero value.

    Args:
        table (pd.DataFrame): table whose numeric columns will be rounded.

    Raises:
        ValueError: if a column cannot be rouned.

    Returns:
        pd.DataFrame: data frame with rounded values.
    """
    decimals = {}
    # Iterate through numeric columns
    for c in table.select_dtypes(include=np.number).columns:
        # Non-zero minimum
        val = table.loc[table[c] > 0, c].median()
        # if val is NaN, then this col has only NaN or only 0 values
        if pd.isna(val):
            decimals[c] = 4
        elif val > 1:
            decimals[c] = 2
        else:
            try:
                decimals[c] = abs(math.floor(math.log10(val))) + 2
            except ValueError:
                logger.error(val)
                raise ValueError
    return table.round(decimals)


def write_power_sector_results(
    ba_fuel_data: pd.DataFrame,
    year: int,
    path_prefix: str,
    skip_outputs: bool,
    include_hourly: bool,
):
    """Helper function to write combined data by BA

    Args:
        ba_fuel_data (pd.DataFrame): plant data aggregated by BA and fuel type.
        year (int): a four-digit year indicating when data were taken.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        skip_outputs (bool): whether to save power sector results or not.
        include_hourly (bool): whether to include hourly results in addition to
            monthly and yearly results
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
            if not isinstance(ba, str):
                logger.warning(
                    f"not aggregating {sum(ba_fuel_data.ba_code.isna())} plants "
                    f"with numeric BA {ba}"
                )
                continue

            # filter the data for a single BA
            ba_table = ba_fuel_data[ba_fuel_data["ba_code"] == ba].drop(
                columns="ba_code"
            )

            if include_hourly:
                # convert the datetime_utc column back to a datetime
                ba_table["datetime_utc"] = (
                    pd.to_datetime(ba_table["datetime_utc"])
                    .dt.tz_localize(None)
                    .astype("datetime64[s]")
                    .dt.tz_localize("UTC")
                )

                # calculate a total for the BA
                # grouping by datetime_utc and report_date will create some duplicate
                # datetime_utc values for certain bas where there are plants located
                # in multiple timezones the report date column is necessary for monthly
                # aggregation, but we will have to remove it and group values by
                # datetime_utc for the hourly calculations
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

                # create a dataframe for the hourly values that groups duplicate
                # datetime_utc values
                ba_table_hourly = ba_table.copy().drop(columns=["report_date"])
                ba_table_hourly = (
                    ba_table_hourly.groupby(["fuel_category", "datetime_utc"])
                    .sum()
                    .reset_index()
                )

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
                    year,
                    df=ba_table_hourly,
                    df_name="power sector hourly ba table",
                    keys=["fuel_category"],
                )
                validation.check_for_complete_hourly_timeseries(
                    ba_table_hourly,
                    "power sector hourly ba table",
                    ["fuel_category"],
                    "year",
                )

                # export to a csv
                output_to_results(
                    ba_table_hourly,
                    year,
                    ba,
                    "power_sector_data/hourly/",
                    path_prefix,
                    skip_outputs,
                )
            else:
                ba_total = (
                    ba_table.groupby(["report_date"], dropna=False)[data_columns]
                    .sum()
                    .reset_index()
                )
                ba_total["fuel_category"] = "total"

                # concat the totals to the fuel-specific totals
                ba_table = pd.concat([ba_table, ba_total], axis=0, ignore_index=True)

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
                year,
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
                year,
                ba,
                "power_sector_data/annual/",
                path_prefix,
                skip_outputs,
            )


def add_generated_emission_rate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add generated emission rate columns to the input data frame.

    Args:
        df (pd.DataFrame): data frame with emission data.

    Returns:
        pd.DataFrame: data frame with the additional generated emission rate columns.
    """
    for emission_type in ["_for_electricity", "_for_electricity_adjusted"]:
        for emission in ["co2", "ch4", "n2o", "co2e", "nox", "so2"]:
            col_name = f"generated_{emission}_rate_lb_per_mwh{emission_type}"
            df[col_name] = (
                (df[f"{emission}_mass_lb{emission_type}"] / df["net_generation_mwh"])
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
