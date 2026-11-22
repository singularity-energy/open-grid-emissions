import calendar
import math
import pandas as pd
import numpy as np
import shutil
import os

import oge.load_data as load_data
from oge.column_checks import (
    check_columns,
    DATA_COLUMNS,
    HG_DATA_COLUMNS,
    STORAGE_DATA_COLUMNS,
)
import oge.validation as validation
from oge.filepaths import outputs_folder, results_folder, data_folder
from oge.helpers import (
    assign_fleet_to_subplant_data,
    flag_missing_hg,
    remove_incomplete_hg,
)
from oge.logging_util import get_logger
from oge.constants import (
    ConversionFactors,
    CLEAN_FUELS,
    RESIDUAL_MIX_FUEL_CATEGORIES,
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

# plant factor columns of the annual plant results. They are blank where they don't
# apply (e.g. the nonbaseload factor of a wind plant), and are rounded to a fixed
# number of decimals because they are multiplied by generation and emissions
PLANT_FACTOR_COLUMNS = ["capacity_factor", "nonbaseload_factor"]


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
    logger.info(f"Exporting {file_name} to data/results/{path_prefix}{subfolder}")

    if include_metric:
        metric = convert_results(df)
        metric = round_table(metric)
    df = round_table(df)

    # Check for negatives after rounding
    validation.test_for_negative_values(df, year)
    # Check that there are no missing values. Energy storage data is only reported for
    # energy storage resources, and Hg data is only reported for some subplants, so
    # they are expected to be missing for other resources
    validation.test_for_missing_values(
        df, skip_cols=STORAGE_DATA_COLUMNS + HG_DATA_COLUMNS + PLANT_FACTOR_COLUMNS
    )

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
    plant_factors: pd.DataFrame,
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
        plant_factors (pd.DataFrame): each plant's annual capacity factor and
            nonbaseload factor, from `calculate_plant_factors`, which are added to the
            annual plant data.
    """
    if not skip_outputs:
        groupby_cols = ["plant_id_eia"]
        if plant_part == "subplant":
            groupby_cols += ["subplant_id"]
        if resolution == "monthly":
            groupby_cols += ["report_date"]

        # group data to specified groups. Only sum the data columns so that numeric
        # identifiers (e.g. subplant_id in plant-level data) are not summed
        # Hg is only summed where it is complete for every subplant-month in the group
        data_columns = (
            DATA_COLUMNS
            + HG_DATA_COLUMNS
            + ["missing_hg"]
            + [
                col
                for col in STORAGE_DATA_COLUMNS
                if col in monthly_subplant_data.columns
            ]
        )
        df = (
            monthly_subplant_data.assign(
                missing_hg=flag_missing_hg(monthly_subplant_data)
            )
            .groupby(groupby_cols, dropna=False)[data_columns]
            .sum(min_count=1)
            .reset_index()
        )
        df = remove_incomplete_hg(df)

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
            ).merge(plant_factors, how="left", on="plant_id_eia", validate="1:1")

            # rearrange columns
            df = df[
                [
                    "plant_id_eia",
                    "plant_name_eia",
                    "fuel_category",
                    "capacity_mw",
                    "capacity_factor",
                    "nonbaseload_factor",
                    "ba_code",
                    "city",
                    "county",
                    "state",
                ]
                + DATA_COLUMNS
                + HG_DATA_COLUMNS
                + [column for column in STORAGE_DATA_COLUMNS if column in df.columns]
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


def calculate_plant_factors(
    monthly_subplant_data: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    eia923_allocated: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Calculates each plant's annual capacity factor and nonbaseload factor.

    This is the coordinating function for the plant factors of the annual plant data.
    The capacity factor only counts the months when each generator was operating (see
    `calculate_operating_capacity_mwh`), and is zero where net generation is negative,
    as in eGRID. The nonbaseload factor follows eGRID's method (see
    `calculate_nonbaseload_factor`).

    Args:
        monthly_subplant_data (pd.DataFrame): monthly data for all subplants, whose net
            generation is summed to each plant's annual net generation.
        plant_attributes (pd.DataFrame): the plant static attributes table, with each
            plant's 'capacity_mw' and 'plant_primary_fuel'.
        eia923_allocated (pd.DataFrame): monthly EIA-923 data allocated to generators.
        year (int): data year.

    Returns:
        pd.DataFrame: 'plant_id_eia', 'capacity_factor' and 'nonbaseload_factor' of
            each plant in `monthly_subplant_data`.
    """
    # sum each plant's monthly subplant net generation to its annual net generation,
    # in the same way as the annual plant data
    plants = (
        monthly_subplant_data.groupby("plant_id_eia", dropna=False)[
            "net_generation_mwh"
        ]
        .sum(min_count=1)
        .reset_index()
    )

    # add the nameplate capacity, which is used where a plant has no generator data,
    # and the primary fuel, which decides whether a plant has a nonbaseload factor.
    # The name and fuel category are only used to describe plants in warnings
    plants = plants.merge(
        plant_attributes[
            [
                "plant_id_eia",
                "plant_name_eia",
                "fuel_category",
                "plant_primary_fuel",
                "capacity_mw",
            ]
        ],
        how="left",
        on="plant_id_eia",
        validate="1:1",
    )

    # add the generation at nameplate capacity in the months when each of the plant's
    # generators was operating or generating, which is the denominator of the
    # capacity factor
    operating_capacity = calculate_operating_capacity_mwh(eia923_allocated, year)
    plants = plants.merge(
        operating_capacity, how="left", on="plant_id_eia", validate="1:1"
    )

    # capacity factors are zero where net generation is negative, as in eGRID
    plants["capacity_factor"] = calculate_capacity_factor(
        plants["net_generation_mwh"],
        plants["capacity_mw"],
        year,
        operating_capacity_mwh=plants["operating_capacity_mwh"],
    ).clip(lower=0)

    # the nonbaseload factor is the share of the plant's generation and emissions that
    # is considered nonbaseload, based on its capacity factor and primary fuel
    plants["nonbaseload_factor"] = calculate_nonbaseload_factor(
        plants["capacity_factor"], plants["plant_primary_fuel"]
    )

    # capacity factors above 1 point to missing capacity or net generation that is too
    # high, which are not corrected but are listed in a warning
    validation.check_for_capacity_factors_above_one(plants)

    return plants[["plant_id_eia", "capacity_factor", "nonbaseload_factor"]]


def calculate_capacity_factor(
    net_generation: pd.Series,
    capacity: pd.Series,
    year: int,
    operating_capacity_mwh: pd.Series | None = None,
) -> pd.Series:
    """Calculates capacity factors as net generation divided by the generation at
    nameplate capacity in the hours when the capacity was operating.

    Where `operating_capacity_mwh` is given and positive, it is the generation at
    nameplate capacity in the months when each generator operated (see
    `calculate_operating_capacity_mwh`), so that plants that only operated part of the
    year are not counted as idle for the rest of it. Otherwise, the capacity is assumed
    to operate in every hour of the year, i.e. 8,784 hours in a leap year.

    Capacity factors are negative where net generation is negative. eGRID keeps these
    negative capacity factors for generators, but sets them to zero for plants.

    Args:
        net_generation (pd.Series): annual net generation in MWh.
        capacity (pd.Series): nameplate capacity in MW.
        year (int): data year.
        operating_capacity_mwh (pd.Series, optional): generation at nameplate capacity
            in the months when the capacity was operating, in MWh. Defaults to None.

    Returns:
        pd.Series: capacity factors, blank where the capacity is zero or blank.
    """
    hours = 24 * (366 if calendar.isleap(year) else 365)
    full_year_capacity_mwh = capacity.where(capacity > 0) * hours
    if operating_capacity_mwh is not None:
        full_year_capacity_mwh = operating_capacity_mwh.where(
            operating_capacity_mwh > 0, full_year_capacity_mwh
        )
    return net_generation / full_year_capacity_mwh


def calculate_operating_capacity_mwh(
    eia923_allocated: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Calculates each plant's generation at nameplate capacity in the months when each
    of its generators was operating, which is the denominator of its capacity factor.

    A generator is counted as operating in a month if either:
    - its EIA-860 dates cover the month: it is existing or retired, its operating date
      is blank or on or before the end of the month, and its retirement date is blank
      or on or after the start of the month. A generator that starts or retires
      during a month is counted for the whole month, since EIA-860 dates are by month.
    - it produced electricity in the month (positive net generation in the EIA-923
      data allocated to generators), even if it has not reached its commercial
      operation date. For example, generators that are still testing (EIA-860 status
      TS) are reported as proposed, but can generate for months before they are
      considered operating.

    Args:
        eia923_allocated (pd.DataFrame): monthly EIA-923 data allocated to generators,
            with 'plant_id_eia', 'generator_id', 'report_date' and
            'net_generation_mwh' columns.
        year (int): data year.

    Returns:
        pd.DataFrame: 'plant_id_eia' and 'operating_capacity_mwh', the sum of each
            generator's nameplate capacity times the hours of the months in which it
            operated. Plants whose generators have no capacity are not included.
    """
    generators = load_data.load_pudl_table(
        "out_eia__yearly_generators",
        year,
        columns=[
            "plant_id_eia",
            "generator_id",
            "capacity_mw",
            "operational_status",
            "generator_operating_date",
            "generator_retirement_date",
        ],
    )
    generator_months = generators.merge(
        pd.DataFrame(
            {"report_date": pd.date_range(f"{year}-01-01", periods=12, freq="MS")}
        ),
        how="cross",
    )
    month_end = generator_months["report_date"] + pd.offsets.MonthEnd(0)
    in_operation = (
        generator_months["operational_status"].isin(["existing", "retired"])
        & (
            generator_months["generator_operating_date"].isna()
            | (generator_months["generator_operating_date"] <= month_end)
        )
        & (
            generator_months["generator_retirement_date"].isna()
            | (
                generator_months["generator_retirement_date"]
                >= generator_months["report_date"]
            )
        )
    ).fillna(False)

    # months in which each generator produced electricity
    generation = (
        eia923_allocated.assign(
            report_date=pd.to_datetime(eia923_allocated["report_date"])
        )
        .groupby(["plant_id_eia", "generator_id", "report_date"], dropna=False)[
            "net_generation_mwh"
        ]
        .sum(min_count=1)
        .reset_index()
    )
    generator_months = generator_months.merge(
        generation,
        how="left",
        on=["plant_id_eia", "generator_id", "report_date"],
        validate="1:1",
    )
    generating = generator_months["net_generation_mwh"].fillna(0) > 0

    generator_months["operating_capacity_mwh"] = (
        generator_months["capacity_mw"]
        * generator_months["report_date"].dt.days_in_month
        * 24
    ).where(in_operation | generating)
    return (
        generator_months.groupby("plant_id_eia", dropna=False)["operating_capacity_mwh"]
        .sum(min_count=1)
        .dropna()
        .reset_index()
    )


def calculate_nonbaseload_factor(
    capacity_factor: pd.Series, primary_fuel: pd.Series
) -> pd.Series:
    """Calculates the nonbaseload factor of each plant, which is the share of its
    generation and emissions that is considered nonbaseload, from its annual capacity
    factor, as in eGRID.

    The capacity factor stands in for where a plant falls in the dispatch order: plants
    with a capacity factor below 0.2 are fully nonbaseload (a factor of 1), plants
    above 0.8 are fully baseload (a factor of 0), and in between the factor is
    -5/3 * capacity factor + 4/3, which falls linearly from 1 to 0.

    The factor is blank for plants whose capacity factor is blank or zero, which have no
    nonbaseload generation, and for plants whose primary fuel is not dispatched to
    follow load, so that their capacity factor says nothing about whether they are on
    the margin. As in eGRID's code (but not its technical guide), biomass plants are not
    excluded.

    Args:
        capacity_factor (pd.Series): annual capacity factor of each plant.
        primary_fuel (pd.Series): primary fuel (energy source code) of each plant.

    Returns:
        pd.Series: nonbaseload factor of each plant, between 0 and 1.
    """
    excluded_fuels = [
        "GEO",  # geothermal
        "MWH",  # energy storage
        "NUC",  # nuclear
        "PUR",  # purchased steam
        "SUN",  # solar
        "WAT",  # water (hydro, including pumped storage)
        "WND",  # wind
    ]
    factor = pd.Series(
        np.select(
            [capacity_factor < 0.2, capacity_factor > 0.8],
            [1.0, 0.0],
            default=-5 / 3 * capacity_factor + 4 / 3,
        ),
        index=capacity_factor.index,
    )
    return factor.where((capacity_factor > 0) & ~primary_fuel.isin(excluded_fuels))


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


def write_plant_metadata(
    plant_static_attributes: pd.DataFrame,
    eia923_allocated: pd.DataFrame,
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
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
    Rounding for each column is based on the median non-zero value, except for the
    columns of `PLANT_FACTOR_COLUMNS`, which have four decimals.

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
        if c in PLANT_FACTOR_COLUMNS:
            decimals[c] = 4
            continue
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
    """Helper function to write combined data by BA and for the entire US.

    The national ("US") data is summed from the same data used for the BA files at
    each resolution, so it equals the sum of the BA files.

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
    # sum the data for all BAs to the national level
    national_keys = ["fuel_category", "report_date"] + (
        ["datetime_utc"] if include_hourly else []
    )
    national_columns = DATA_COLUMNS + [
        column for column in STORAGE_DATA_COLUMNS if column in fleet_data.columns
    ]
    national_table = (
        fleet_data[fleet_data["ba_code"].map(lambda ba: isinstance(ba, str))]
        .groupby(national_keys, dropna=False)[national_columns]
        .sum(min_count=1)
        .reset_index()
    )

    for ba in list(fleet_data.ba_code.unique()) + ["US"]:
        if not isinstance(ba, str):
            logger.warning(
                f"not aggregating {sum(fleet_data.ba_code.isna())} plants "
                f"with numeric BA {ba}"
            )
            continue

        # filter the data for a single BA
        if ba == "US":
            ba_table = national_table
        else:
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
                .sum(min_count=1)
                .reset_index()
            )
            ba_total["fuel_category"] = "total"
            ba_residual = calculate_residual_mix(
                ba_table, ba_total, ["datetime_utc", "report_date"]
            )

            # concat the totals and residual mix to the fuel-specific totals
            ba_table = pd.concat(
                [ba_table, ba_residual, ba_total], axis=0, ignore_index=True
            )

            # create a dataframe for the hourly values that groups duplicate
            # datetime_utc values
            ba_table_hourly = ba_table.copy().drop(columns=["report_date"])
            ba_table_hourly = (
                ba_table_hourly.groupby(["fuel_category", "datetime_utc"])
                .sum(min_count=1)
                .reset_index()
            )

            # output the hourly data
            ba_table_hourly = add_generated_emission_rate_columns(ba_table_hourly)

            # create a local datetime column. The US spans multiple timezones, so the
            # national data only has a UTC datetime
            datetime_columns = ["datetime_utc"]
            if ba != "US":
                try:
                    local_tz = load_data.ba_timezone(ba, "local")
                    ba_table_hourly["datetime_local"] = ba_table_hourly[
                        "datetime_utc"
                    ].dt.tz_convert(local_tz)
                except ValueError:
                    ba_table_hourly["datetime_local"] = pd.NaT
                datetime_columns = ["datetime_local", "datetime_utc"]

            # re-order columns
            ba_table_hourly = ba_table_hourly[
                ["fuel_category"]
                + datetime_columns
                + DATA_COLUMNS
                + GENERATED_EMISSION_RATE_COLS
            ]
            ba_table_hourly = sort_fuel_categories(ba_table_hourly, ["datetime_utc"])

            validation.validate_unique_datetimes(
                year,
                df=ba_table_hourly,
                df_name="power sector hourly ba table",
                keys=["fuel_category"],
            )
            # each BA's data covers its own local year, so the national data spans
            # the local years of all BAs, and the first and last few hours only
            # include data from some BAs
            if ba != "US":
                validation.check_for_complete_hourly_timeseries(
                    ba_table_hourly,
                    f"power sector hourly {ba} table",
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
            # include any energy storage data, which is only available at the monthly
            # level
            combined_data_columns = DATA_COLUMNS + [
                column for column in STORAGE_DATA_COLUMNS if column in ba_table.columns
            ]
            ba_total = (
                ba_table.groupby(["report_date"], dropna=False)[combined_data_columns]
                .sum(min_count=1)
                .reset_index()
            )
            ba_total["fuel_category"] = "total"
            ba_residual = calculate_residual_mix(ba_table, ba_total, ["report_date"])

            # concat the totals and residual mix to the fuel-specific totals
            ba_table = pd.concat(
                [ba_table, ba_residual, ba_total], axis=0, ignore_index=True
            )

            agg = {}
            if include_monthly:
                agg["monthly"] = ["fuel_category", "report_date"]

            if include_annual:
                agg["annual"] = ["fuel_category"]

            for agg_level, groupby_cols in agg.items():
                # aggregate data
                ba_table = (
                    ba_table.groupby(groupby_cols, dropna=False)
                    .sum(numeric_only=True, min_count=1)
                    .reset_index()
                )
                ba_table = add_generated_emission_rate_columns(ba_table)
                # re-order columns
                ba_table = ba_table[
                    groupby_cols + combined_data_columns + GENERATED_EMISSION_RATE_COLS
                ]
                ba_table = sort_fuel_categories(ba_table, groupby_cols[1:])
                # the annual national file is always written, even if skip_outputs is
                # True, because it is used to calculate consumption-based emissions
                output_to_results(
                    ba_table,
                    year,
                    ba,
                    f"power_sector_data/{agg_level}/",
                    path_prefix,
                    skip_outputs and not (ba == "US" and agg_level == "annual"),
                )


def sort_fuel_categories(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Sort power sector results by fuel category and `keys`.

    Individual fuel categories are sorted alphabetically, followed by the "residual"
    and "total" rows, so that the aggregated rows are always at the bottom of the
    file.

    Args:
        df (pd.DataFrame): power sector results with a fuel_category column.
        keys (list[str]): columns to sort by within each fuel category.

    Returns:
        pd.DataFrame: sorted data frame with a reset index.
    """
    aggregate_order = {"residual": 1, "total": 2}
    return (
        df.assign(fuel_order=df["fuel_category"].map(aggregate_order).fillna(0))
        .sort_values(["fuel_order", "fuel_category"] + keys)
        .drop(columns="fuel_order")
        .reset_index(drop=True)
    )


def calculate_residual_mix(
    ba_table: pd.DataFrame, ba_total: pd.DataFrame, keys: list[str]
) -> pd.DataFrame:
    """Calculate the fossil-only residual mix as a "residual" fuel category.

    The fossil-only residual mix approximates the residual mix by summing only the
    fuel categories in `RESIDUAL_MIX_FUEL_CATEGORIES`. The residual mix covers the
    same keys as the total data, with zeros where there is no generation from these
    fuel categories, so that every BA has a residual row even if it has no fossil
    generators. Energy storage columns are left blank since storage is not included
    in the residual mix.

    Args:
        ba_table (pd.DataFrame): data for a single BA, by fuel category and `keys`.
        ba_total (pd.DataFrame): total data for the BA, by `keys`.
        keys (list[str]): columns other than fuel_category identifying each record.

    Returns:
        pd.DataFrame: residual mix data with the same keys as `ba_total`.
    """
    ba_residual = (
        ba_table[ba_table["fuel_category"].isin(RESIDUAL_MIX_FUEL_CATEGORIES)]
        .groupby(keys, dropna=False)[DATA_COLUMNS]
        .sum(min_count=1)
        .reindex(ba_total.set_index(keys).index, fill_value=0)
        .reset_index()
    )
    ba_residual["fuel_category"] = "residual"
    return ba_residual


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
    eia_only_data = validation.identify_reporting_frequency(eia_only_data, year)
    partial_cems_subplant = validation.identify_reporting_frequency(
        partial_cems_subplant, year
    )
    partial_cems_plant = validation.identify_reporting_frequency(
        partial_cems_plant, year
    )

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
        if len(df) == 0:  # Empty df
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


def summarize_annually_reported_eia_data(
    eia923_allocated: pd.DataFrame,
    year: int,
    monthly_storage_data: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Creates table summarizing the percent of final data from annually-reported EIA data.

    Args:
        eia923_allocated (pd.DataFrame): allocated EIA-923 data by generator.
        year (int): the data year.
        monthly_storage_data (pd.DataFrame | None, optional): monthly energy storage
            charging and discharging data by subplant. If provided, the percent of
            energy storage data from annually-reported EIA data is also summarized.
            Defaults to None.

    Returns:
        pd.DataFrame: table summarizing the percent of each data column that comes from
            annually-reported EIA data.
    """

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

    eia_data = validation.identify_reporting_frequency(eia923_allocated, year)

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

    # add the percent of energy storage data from annually-reported EIA data. All energy
    # storage data comes from EIA-923, so the percent of input data and output data is
    # the same, and it is not combined with CEMS data
    if monthly_storage_data is not None:
        storage_data = validation.identify_reporting_frequency(
            monthly_storage_data, year
        )
        storage_from_annual = (
            storage_data.loc[
                storage_data["eia_data_resolution"] == "annual", STORAGE_DATA_COLUMNS
            ].sum(min_count=1)
            / storage_data[STORAGE_DATA_COLUMNS].sum(min_count=1)
            * 100
        ).round(2)
        summary_rows = [
            row
            for row in [
                "% of EIA-923 input data from EIA annual reporters",
                "% of output data from EIA annual reporters",
            ]
            if row in annual_data_summary.index
        ]
        for column in STORAGE_DATA_COLUMNS:
            annual_data_summary[column] = np.nan
            annual_data_summary.loc[summary_rows, column] = storage_from_annual[column]

    annual_data_summary = annual_data_summary.reset_index()

    return annual_data_summary
