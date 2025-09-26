import math
import pandas as pd
import numpy as np
import shutil
import os

import oge.load_data as load_data
from oge.column_checks import check_columns, DATA_COLUMNS
import oge.validation as validation
from oge.filepaths import outputs_folder, results_folder, data_folder
from oge.helpers import assign_fleet_to_subplant_data
from oge.logging_util import get_logger
from oge.constants import (
    ConversionFactors,
    CLEAN_FUELS,
    earliest_validated_year,
    earliest_hourly_data_year,
    latest_validated_year,
    current_early_release_year,
)

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
        zip_results_for_s3()
        zip_data_for_zenodo(year)


def zip_results_for_s3():
    """Zips results directories that contain more than a single file for hosting on an
    Amazon S3 bucket.

    """
    os.makedirs(data_folder("s3_upload"), exist_ok=True)
    historical_years = list(range(earliest_validated_year, earliest_hourly_data_year))
    year_range = f"{earliest_validated_year}-{earliest_hourly_data_year - 1}"
    for data_type in ["power_sector_data", "plant_data"]:
        for aggregation in ["monthly", "annual"]:
            for unit in ["metric_units", "us_units"]:
                logger.info(
                    f"zipping {year_range}_{data_type}_{aggregation}_{unit} for s3"
                )
                for year in historical_years:
                    # copy the annual file to a combined folder
                    shutil.copytree(
                        (results_folder(f"{year}/{data_type}/{aggregation}/{unit}")),
                        data_folder(
                            f"s3_upload/{year_range}_{data_type}_{aggregation}_{unit}/{year}"
                        ),
                    )
                # now create an archive
                shutil.make_archive(
                    data_folder(
                        f"s3_upload/{year_range}_{data_type}_{aggregation}_{unit}"
                    ),
                    "zip",
                    root_dir=data_folder(
                        f"s3_upload/{year_range}_{data_type}_{aggregation}_{unit}"
                    ),
                )
                shutil.rmtree(
                    data_folder(
                        f"s3_upload/{year_range}_{data_type}_{aggregation}_{unit}"
                    )
                )
    # move and rename the plant attributes files
    os.makedirs(
        data_folder(f"s3_upload/{year_range}_plant_attributes"),
        exist_ok=True,
    )
    for year in historical_years:
        # data quality metrics
        shutil.copytree(
            (results_folder(f"{year}/data_quality_metrics")),
            data_folder(f"s3_upload/{year_range}_data_quality_metrics/{year}"),
        )

        shutil.copy(
            results_folder(f"{year}/plant_data/plant_static_attributes.csv"),
            data_folder(
                f"s3_upload/{year_range}_plant_attributes/plant_static_attributes_{year}.csv"
            ),
        )
    shutil.make_archive(
        data_folder(f"s3_upload/{year_range}_data_quality_metrics"),
        "zip",
        root_dir=data_folder(f"s3_upload/{year_range}_data_quality_metrics"),
    )
    shutil.rmtree(data_folder(f"s3_upload/{year_range}_data_quality_metrics"))
    shutil.make_archive(
        data_folder(f"s3_upload/{year_range}_plant_attributes"),
        "zip",
        root_dir=data_folder(f"s3_upload/{year_range}_plant_attributes"),
    )
    shutil.rmtree(data_folder(f"s3_upload/{year_range}_plant_attributes"))
    for year in range(2019, max(latest_validated_year, current_early_release_year) + 1):
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
                        folder = f"{results_folder()}/{year}/{data_type}/{aggregation}/{unit}"
                        shutil.make_archive(
                            data_folder(
                                f"s3_upload/{year}_{data_type}_{aggregation}_{unit}"
                            ),
                            "zip",
                            root_dir=folder,
                            # base_dir="",
                        )
        # move and rename the plant attributes files
        shutil.copy(
            f"{results_folder()}/{year}/plant_data/plant_static_attributes.csv",
            data_folder(f"s3_upload/plant_static_attributes_{year}.csv"),
        )
        # archive the data quality metrics
        shutil.make_archive(
            data_folder(f"s3_upload/{year}_data_quality_metrics"),
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
    check_columns(df, file_name)
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


def write_plant_data_to_results(
    monthly_subplant_data: pd.DataFrame,
    year: int,
    path_prefix: str,
    plant_part: str,
    resolution: str,
    skip_outputs: bool,
    plant_attributes: pd.DataFrame,
):
    """
    Helper function to write subplant or plant data at the monthly or annual resolution.

    Args:
        monthly_subplant_data (pd.DataFrame): Combined monthly data for all subplants
        year (int): a four-digit year indicating when data were taken.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        plant_part (str): either "plant" or "subplant"
        resolution (str): temporal resolution. Wither 'hourly', 'monthly' or 'annual'.
        skip_outputs (bool): whether to save data or not.
        plant_attributes (pd.DataFrame): the plant static attributes table.
    """
    if not skip_outputs:
        groupby_cols = ["plant_id_eia"]
        if plant_part == "subplant":
            groupby_cols += ["subplant_id"]
        if resolution == "monthly":
            groupby_cols += ["report_date"]

        # group data to specified groups
        df = (
            monthly_subplant_data.groupby(groupby_cols, dropna=False)
            .sum(numeric_only=True)
            .reset_index()
        )

        # add some basic plant data to the annual plant output table
        if resolution == "annual" and plant_part == "plant":
            df = df.merge(
                plant_attributes[
                    [
                        "plant_id_eia",
                        "plant_name_eia",
                        "fuel_category",
                        "capacity_mw",
                        "ba_code",
                        "city",
                        "county",
                        "state",
                    ]
                ],
                how="left",
                on="plant_id_eia",
                validate="m:1",
            )

            # rearrange columns
            df = df[
                [
                    "plant_id_eia",
                    "plant_name_eia",
                    "fuel_category",
                    "capacity_mw",
                    "ba_code",
                    "city",
                    "county",
                    "state",
                ]
                + DATA_COLUMNS
            ]

        # calculate emission rates
        df = add_generated_emission_rate_columns(df)

        # check for anomalous looking co2 rates
        validation.check_for_anomalous_co2_factors(
            df, plant_attributes, year, min_threshold=10, max_threshold=15000
        )

        # output monthly data
        output_to_results(
            df,
            year,
            file_name=f"{plant_part}_data",
            subfolder=f"plant_data/{resolution}/",
            path_prefix=path_prefix,
            skip_outputs=skip_outputs,
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


def write_national_fleet_averages(
    ba_fuel_data: pd.DataFrame, year: int, path_prefix: str, skip_outputs: bool
):
    """Outputs annual, national-average fleet-level emissions data

    Args:
        ba_fuel_data (pd.DataFrame): plant data aggregated by BA and fuel type.
        year (int): a four-digit year indicating when the data were taken.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        skip_outputs (bool): whether to save data or not.
    """
    if not skip_outputs:
        # sum all of the columns by fuel before calculating emission rates
        national_avg = (
            ba_fuel_data.groupby(["fuel_category"])[DATA_COLUMNS]
            .sum(numeric_only=True)
            .reset_index()
        )

        # Add row for total
        national_total = pd.DataFrame(national_avg[DATA_COLUMNS].sum()).T
        national_total["fuel_category"] = "total"

        # concat the totals to the fuel-specific totals
        national_avg = pd.concat(
            [national_avg, national_total], axis=0, ignore_index=True
        )

        national_avg = add_generated_emission_rate_columns(national_avg)

        output_to_results(
            national_avg,
            year,
            "US",
            "power_sector_data/annual/",
            path_prefix,
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

    monthly_eia_meta["data_source"] = "EIA"

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

    check_columns(metadata, "plant_metadata")

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
    fleet_data: pd.DataFrame,
    year: int,
    path_prefix: str,
    skip_outputs: bool,
    include_hourly: bool,
    include_monthly: bool,
    include_annual: bool,
):
    """Helper function to write combined data by BA

    Args:
        fleet_data (pd.DataFrame): subplant data aggregated by BA and fuel type.
        year (int): a four-digit year indicating when data were taken.
        path_prefix (str): name of base directory prefixing directory where data will
            be saved.
        skip_outputs (bool): whether to save power sector results or not.
        include_hourly (bool): whether to include hourly results.
        include_monthly (bool): whether to include monthly results.
        include_annual (bool): whether to include annual results.
    """

    if not skip_outputs:
        for ba in list(fleet_data.ba_code.unique()):
            if not isinstance(ba, str):
                logger.warning(
                    f"not aggregating {sum(fleet_data.ba_code.isna())} plants "
                    f"with numeric BA {ba}"
                )
                continue

            # filter the data for a single BA
            ba_table = fleet_data[fleet_data["ba_code"] == ba].drop(columns="ba_code")

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
                        DATA_COLUMNS
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
                    + DATA_COLUMNS
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
            elif include_monthly or include_annual:
                ba_total = (
                    ba_table.groupby(["report_date"], dropna=False)[DATA_COLUMNS]
                    .sum()
                    .reset_index()
                )
                ba_total["fuel_category"] = "total"

                # concat the totals to the fuel-specific totals
                ba_table = pd.concat([ba_table, ba_total], axis=0, ignore_index=True)

                agg = {}
                if include_monthly:
                    agg["monthly"] = ["fuel_category", "report_date"]

                if include_annual:
                    agg["annual"] = ["fuel_category"]

                for agg_level, groupby_cols in agg.items():
                    # aggregate data
                    ba_table = (
                        ba_table.groupby(groupby_cols, dropna=False)
                        .sum(numeric_only=True)
                        .reset_index()
                    )
                    ba_table = add_generated_emission_rate_columns(ba_table)
                    # re-order columns
                    ba_table = ba_table[
                        groupby_cols + DATA_COLUMNS + GENERATED_EMISSION_RATE_COLS
                    ]
                    output_to_results(
                        ba_table,
                        year,
                        ba,
                        f"power_sector_data/{agg_level}/",
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


def identify_percent_of_data_by_input_source(
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
    eia_only_data: pd.DataFrame,
    year: int,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
) -> pd.DataFrame:
    """For each BA, identifies the percentage of input data (at the subplant level) from
    one of 4 sources, for each data column:
        - cems_hourly
        - eia_annual
        - eia_monthly
        - eia_multiple

    Args:
        cems (pd.DataFrame): used to identify the data source
        partial_cems_subplant (pd.DataFrame): used to identify the data source
        partial_cems_plant (pd.DataFrame): used to identify the data source
        eia_only_data (pd.DataFrame): used to identify the data source
        year (int): the data year
        plant_attributes (pd.DataFrame): used to assign fleet keys
        primary_fuel_table (pd.DataFrame): used to assign fleet keys

    Returns:
        pd.DataFrame: Table of data source percentages
    """

    columns_to_use = [
        "net_generation_mwh",
        "emitting_net_generation_mwh",
        "co2_mass_lb",
        "co2_mass_lb_for_electricity",
        "co2e_mass_lb",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb",
        "so2_mass_lb_for_electricity",
    ]

    # add data resolution column to data that is based on EIA
    eia_only_data = identify_reporting_frequency(eia_only_data, year)
    partial_cems_subplant = identify_reporting_frequency(partial_cems_subplant, year)
    partial_cems_plant = identify_reporting_frequency(partial_cems_plant, year)

    # add ba codes and plant primary fuel to all of the data
    eia_only_data = assign_fleet_to_subplant_data(
        eia_only_data, plant_attributes, primary_fuel_table, year
    )
    cems = assign_fleet_to_subplant_data(
        cems, plant_attributes, primary_fuel_table, year, drop_primary_fuel_col=False
    )
    partial_cems_subplant = assign_fleet_to_subplant_data(
        partial_cems_subplant,
        plant_attributes,
        primary_fuel_table,
        year,
        drop_primary_fuel_col=False,
    )
    partial_cems_plant = assign_fleet_to_subplant_data(
        partial_cems_plant,
        plant_attributes,
        primary_fuel_table,
        year,
        drop_primary_fuel_col=False,
    )

    # add a column for fossil-based generation
    # this copies the net generation data if the associated fuel is not clean or
    # geothermal, and otherwise adds a zero use the generator-specific energy source
    # code for the eia data, otherwise use the plant primary fuel
    eia_only_data = eia_only_data.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.energy_source_code.isin(CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
        )
    )
    cems = cems.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.subplant_primary_fuel.isin(CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
        )
    )
    partial_cems_subplant = partial_cems_subplant.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.subplant_primary_fuel.isin(CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
        )
    )
    partial_cems_plant = partial_cems_plant.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.subplant_primary_fuel.isin(CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
        )
    )

    # associate each dataframe with a data source label
    data_sources = {
        "cems": cems,
        "partial_cems_subplant": partial_cems_subplant,
        "partial_cems_plant": partial_cems_plant,
        "eia": eia_only_data,
    }
    # get a count of the number of observations (subplant-hours) from each source
    source_of_input_data = []
    for name, df in data_sources.items():
        if len(df) == 0:  # Empty df. May occur when running `small`
            logger.warning(f"data source {name} has zero entries")
            continue
        if name == "eia":
            subplant_data = df.groupby(
                ["ba_code", "plant_id_eia", "subplant_id", "eia_data_resolution"],
                dropna=False,
            )[columns_to_use].sum()
            # because EIA data is not hourly, we have to multiply the number of subplants by the number of hours in a year
            if year % 4 == 0:
                hours_in_year = 8784
            else:
                hours_in_year = 8760
            subplant_data["subplant_hours"] = hours_in_year
            # group the data by resolution
            subplant_data = (
                subplant_data.reset_index()
                .groupby(["ba_code", "eia_data_resolution"], dropna=False)[
                    ["subplant_hours"] + columns_to_use
                ]
                .sum()
                .reset_index()
            )
            subplant_data = subplant_data.rename(
                columns={"eia_data_resolution": "source"}
            )
            subplant_data["source"] = subplant_data["source"].replace(
                {
                    "annual": "eia_annual",
                    "monthly": "eia_monthly",
                    "multiple": "eia_multiple",
                }
            )
            source_of_input_data.append(subplant_data)
        # for the partial cems data
        elif (name == "partial_cems_subplant") | (name == "partial_cems_plant"):
            subplant_data = df.groupby(
                [
                    "ba_code",
                    "plant_id_eia",
                    "subplant_id",
                    "datetime_utc",
                    "eia_data_resolution",
                ],
                dropna=False,
            )[columns_to_use].sum()
            subplant_data["subplant_hours"] = 1
            # group the data by resolution
            subplant_data = (
                subplant_data.reset_index()
                .groupby(["ba_code", "eia_data_resolution"], dropna=False)[
                    ["subplant_hours"] + columns_to_use
                ]
                .sum()
                .reset_index()
            )
            subplant_data = subplant_data.rename(
                columns={"eia_data_resolution": "source"}
            )
            subplant_data["source"] = subplant_data["source"].replace(
                {
                    "annual": "eia_annual",
                    "monthly": "eia_monthly",
                    "multiple": "eia_multiple",
                }
            )
            source_of_input_data.append(subplant_data)
        # for the cems data
        else:
            subplant_data = df.groupby(
                ["ba_code", "plant_id_eia", "subplant_id", "datetime_utc"], dropna=False
            )[columns_to_use].sum()
            subplant_data["subplant_hours"] = 1
            subplant_data["source"] = "cems_hourly"
            # group the data by resolution
            subplant_data = (
                subplant_data.reset_index()
                .groupby(["ba_code", "source"], dropna=False)[
                    ["subplant_hours"] + columns_to_use
                ]
                .sum()
                .reset_index()
            )
            source_of_input_data.append(subplant_data)

    # concat the dataframes together
    source_of_input_data = pd.concat(source_of_input_data, axis=0)

    # groupby and calculate percentages for the entire country
    national_source = source_of_input_data.groupby("source").sum(numeric_only=True)
    national_source = (national_source / national_source.sum(axis=0)).reset_index()
    national_source["ba_code"] = "US Total"

    # calculate percentages by ba
    source_of_input_data = (
        source_of_input_data.groupby(["ba_code", "source"]).sum(numeric_only=True)
        / source_of_input_data.groupby(["ba_code"]).sum(numeric_only=True)
    ).reset_index()
    # concat the national data to the ba data
    source_of_input_data = pd.concat([source_of_input_data, national_source], axis=0)

    return source_of_input_data


def identify_reporting_frequency(eia923_allocated, year):
    """Identifies if EIA data was reported as an annual total or monthly totals.
    Returns input dataframe with `eia_data_resolution` column added"""

    # load data about the respondent frequency for each plant and merge into the EIA-923 data
    plant_frequency = load_data.load_pudl_table(
        "out_eia__yearly_plants",
        year,
        columns=["plant_id_eia", "reporting_frequency_code"],
    )
    plant_frequency["reporting_frequency_code"] = plant_frequency[
        "reporting_frequency_code"
    ].fillna("multiple")
    # rename the column and recode the values
    plant_frequency = plant_frequency.rename(
        columns={"reporting_frequency_code": "eia_data_resolution"}
    )
    plant_frequency["eia_data_resolution"] = plant_frequency[
        "eia_data_resolution"
    ].replace({"A": "annual", "AM": "monthly", "M": "monthly"})
    # merge the data resolution column into the EIA data
    eia_data = eia923_allocated.merge(
        plant_frequency, how="left", on="plant_id_eia", validate="m:1"
    )
    return eia_data


def summarize_annually_reported_eia_data(eia923_allocated, year):
    """Creates table summarizing the percent of final data from annually-reported EIA data."""

    columns_to_summarize = [
        "fuel_consumed_mmbtu",
        "net_generation_mwh",
        "co2_mass_lb",
        "co2_mass_lb_for_electricity",
        "nox_mass_lb",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb",
        "so2_mass_lb_for_electricity",
    ]

    eia_data = identify_reporting_frequency(eia923_allocated, year)

    data_from_annual = (
        eia_data.groupby(["eia_data_resolution"], dropna=False)[
            columns_to_summarize
        ].sum()
        / eia_data[columns_to_summarize].sum()
        * 100
    ).reset_index()

    annual_eia_used = (
        eia_data[eia_data["hourly_data_source"] != "cems"]
        .groupby(["eia_data_resolution"], dropna=False)[columns_to_summarize]
        .sum()
        / eia_data[columns_to_summarize].sum()
        * 100
    ).reset_index()

    multi_source_subplants = (
        eia_data[["plant_id_eia", "subplant_id", "hourly_data_source"]]
        .drop_duplicates()
        .drop(columns="hourly_data_source")
    )
    multi_source_subplants = multi_source_subplants[
        multi_source_subplants.duplicated(subset=["plant_id_eia", "subplant_id"])
    ]
    multi_source_subplants = eia_data.merge(
        multi_source_subplants, how="inner", on=["plant_id_eia", "subplant_id"]
    )
    multi_source_summary = (
        multi_source_subplants.groupby(["eia_data_resolution"], dropna=False)[
            columns_to_summarize
        ].sum()
        / eia_data[columns_to_summarize].sum()
        * 100
    ).reset_index()

    annual_data_summary = pd.concat(
        [
            pd.DataFrame(
                data_from_annual.loc[
                    data_from_annual["eia_data_resolution"] == "annual", :
                ]
                .set_index("eia_data_resolution")
                .rename(
                    index={
                        "annual": "% of EIA-923 input data from EIA annual reporters"
                    }
                )
                .round(2)
            ),
            pd.DataFrame(
                annual_eia_used.loc[
                    annual_eia_used["eia_data_resolution"] == "annual", :
                ]
                .set_index("eia_data_resolution")
                .rename(index={"annual": "% of output data from EIA annual reporters"})
                .round(2)
            ),
            pd.DataFrame(
                multi_source_summary.loc[
                    multi_source_summary["eia_data_resolution"] == "annual", :
                ]
                .set_index("eia_data_resolution")
                .rename(
                    index={
                        "annual": "% of output data mixing CEMS and annually-reported EIA data"
                    }
                )
                .round(2)
            ),
        ],
        axis=0,
    )

    annual_data_summary.rename(columns={"eia_data_resolution": "category"})

    annual_data_summary = annual_data_summary.reset_index()

    return annual_data_summary
