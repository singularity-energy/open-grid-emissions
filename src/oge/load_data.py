import functools
from typing import Callable
import numpy as np
import pandas as pd
import sqlalchemy as sa
import warnings

from pathlib import Path

from oge import PUDL_ENGINE
from oge.column_checks import get_dtypes, apply_dtypes
from oge.constants import (
    CLEAN_FUELS,
    ConversionFactors,
    earliest_data_year,
    earliest_puerto_rico_data_year,
    earliest_mats_data_year,
    earliest_validated_year,
    latest_validated_year,
    current_early_release_year,
    cancelled_status_codes,
)
from oge.filepaths import (
    downloads_folder,
    find_downloaded_file,
    get_pudl_build_version,
    get_pudl_data_store,
    outputs_folder,
    pudl_folder,
    reference_table_folder,
)
from oge.download_data import get_pudl_release
from oge.logging_util import get_logger

logger = get_logger(__name__)


def load_cems_data(year: int) -> pd.DataFrame:
    """Loads CEMS data for the specified year from the PUDL database

    Args:
        year (int): a four-digut year.

    Returns:
        pd.DataFrame: hourly CEMS data.
    """
    # specify the columns to use from the CEMS database
    cems_columns = [
        "plant_id_epa",  # try to load this column to make sure it has been converted to plant_id_eia
        "plant_id_eia",
        "emissions_unit_id_epa",
        "operating_datetime_utc",
        "operating_time_hours",
        "gross_load_mw",
        "steam_load_lbs",
        "co2_mass_tons",
        "co2_mass_measurement_code",
        "nox_mass_lbs",
        "nox_mass_measurement_code",
        "so2_mass_lbs",
        "so2_mass_measurement_code",
        "heat_content_mmbtu",
    ]

    # load the CEMS data
    if PUDL_ENGINE is None:
        pudl_release = get_pudl_release(get_pudl_build_version(), get_pudl_data_store())
        logger.info(f"Loading core_epacems__hourly_emissions from PUDL {pudl_release}")
    cems = pd.read_parquet(
        pudl_folder("core_epacems__hourly_emissions.parquet"),
        filters=[["year", "==", year]],
        columns=cems_columns,
    )
    # convert to tz-naive datetime to allow for dtype application
    cems = apply_dtypes(cems)

    # strip leading zeros from emissions unit id
    cems["emissions_unit_id_epa"] = cems["emissions_unit_id_epa"].str.lstrip("0")

    # update the plant_id_eia column using manual matches
    cems = update_epa_to_eia_map(cems, year)

    cems = cems.rename(
        columns={
            "operating_datetime_utc": "datetime_utc",
            "heat_content_mmbtu": "fuel_consumed_mmbtu",
            "nox_mass_lbs": "nox_mass_lb",
            "so2_mass_lbs": "so2_mass_lb",
            "gross_load_mw": "gross_generation_mwh",  # we will convert this to mwh
        }
    )

    # fill any missing values for steam load with zero
    cems["steam_load_lbs"] = cems["steam_load_lbs"].fillna(0)

    # convert co2 mass in tons to lb
    cems["co2_mass_lb"] = cems["co2_mass_tons"] * ConversionFactors.short_ton_to_lbs

    # re-order columns
    cems = cems[
        [
            "plant_id_eia",
            "emissions_unit_id_epa",
            "datetime_utc",
            "operating_time_hours",
            "gross_generation_mwh",
            "steam_load_lbs",
            "fuel_consumed_mmbtu",
            "co2_mass_lb",
            "nox_mass_lb",
            "so2_mass_lb",
            "plant_id_epa",
            "co2_mass_measurement_code",
            "nox_mass_measurement_code",
            "so2_mass_measurement_code",
        ]
    ]

    # re-code mass_measurement_codes as categoricals to reduce file size
    cems = cems.astype(
        {
            "plant_id_eia": "int",
            "emissions_unit_id_epa": "str",
            "co2_mass_measurement_code": "category",
            "nox_mass_measurement_code": "category",
            "so2_mass_measurement_code": "category",
        }
    )

    return cems


def load_mats_data(year: int) -> pd.DataFrame:
    """Loads hourly mercury (Hg) data reported under the EPA Mercury and Air Toxics
    Standards (MATS) for the specified year from PUDL.

    MATS hourly data only exists for coal- and petroleum coke-fired units that monitor
    Hg with a CEMS or sorbent traps, so most CEMS units have no rows. Missing hours are
    not substituted. No MATS data exists before `constants.earliest_mats_data_year`, so
    an empty table is returned for earlier years.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: hourly Hg mass by unit, keyed on `plant_id_epa`,
            `emissions_unit_id_epa`, and `datetime_utc`.
    """
    output_columns = [
        "plant_id_epa",
        "emissions_unit_id_epa",
        "datetime_utc",
        "hg_mass_lb",
        "hg_mass_measurement_code",
    ]
    if year < earliest_mats_data_year:
        return apply_dtypes(pd.DataFrame(columns=output_columns))

    if PUDL_ENGINE is None:
        pudl_release = get_pudl_release(get_pudl_build_version(), get_pudl_data_store())
        logger.info(f"Loading core_epamats__hourly_emissions from PUDL {pudl_release}")
    mats = pd.read_parquet(
        pudl_folder("core_epamats__hourly_emissions.parquet"),
        filters=[["year", "==", year]],
        columns=[
            "plant_id_epa",
            "emissions_unit_id_epa",
            "operating_datetime_utc",
            "hg_mass_lbs",
            "hg_mass_measurement_code",
        ],
    ).rename(
        columns={
            "operating_datetime_utc": "datetime_utc",
            "hg_mass_lbs": "hg_mass_lb",
        }
    )
    mats = apply_dtypes(mats)

    # strip leading zeros from emissions unit id to match `load_cems_data()`
    mats["emissions_unit_id_epa"] = mats["emissions_unit_id_epa"].str.lstrip("0")

    return mats[output_columns]


def load_cems_ids(year: int | None = None) -> pd.DataFrame:
    """Loads a dataframe of all unique plant_id_eia, emissions_unit_id_epa combinations
    that exist from `constants.earliest_data_year` to `constants.latest_validated_year`.
    This is used in the process of creating subplant_ids to ensure complete coverage.

    If optionally run with a single year, only the cems ids for that year are loaded.

    Args:
        year (int, optional): the year to load the cems ids for. If None, the cems ids
            for all years from `constants.earliest_data_year` to
            `constants.latest_validated_year` are loaded.

    Returns:
        pd.DataFrame: a two column table relating plant_id_eia to emissions_unit_id_epa.
    """
    if year is None:
        min_year = earliest_data_year
        max_year = max(latest_validated_year, current_early_release_year)
    else:
        min_year = year
        max_year = year
    # although we could directly load all years at once from the cems parquet file,
    # this would lead to a memoryerror, so we load one year at a time and drop
    # duplicates before concatenating the next year to the dataframe
    cems_ids = []
    if PUDL_ENGINE is None:
        pudl_release = get_pudl_release(get_pudl_build_version(), get_pudl_data_store())
        logger.info(f"Loading core_epacems__hourly_emissions from PUDL {pudl_release}")
    for year_for_id in range(min_year, max_year + 1):
        cems_id_year = pd.read_parquet(
            pudl_folder("core_epacems__hourly_emissions.parquet"),
            filters=[["year", "==", year_for_id]],
            columns=["plant_id_epa", "plant_id_eia", "emissions_unit_id_epa"],
        ).drop_duplicates()
        cems_id_year = apply_dtypes(cems_id_year)
        # strip leading zeros from emissions unit id
        cems_id_year["emissions_unit_id_epa"] = cems_id_year[
            "emissions_unit_id_epa"
        ].str.lstrip("0")
        # update the plant_id_eia column using manual matches
        cems_id_year = update_epa_to_eia_map(cems_id_year, year_for_id)
        cems_ids.append(cems_id_year)
        # only keep new ids for the given year
        cems_ids = [pd.concat(cems_ids, axis=0).drop_duplicates()]

    cems_ids = (
        pd.concat(cems_ids, axis=0)
        .drop_duplicates()
        .sort_values(by=["plant_id_eia", "emissions_unit_id_epa"])
    )

    cems_ids = apply_dtypes(cems_ids)

    return cems_ids[["plant_id_eia", "emissions_unit_id_epa"]]


def load_complete_eia_generators_for_subplants() -> pd.DataFrame:
    """Loads a dataframe that contains a complete list of generators, including their
    unit ids, prime movers, operating dates, and operating status. This will be used
    when creating subplant IDs to ensure that a complete set of generators is
    represented. Because some of these values are incomplete or missing in the pudl data
    we load these values from the raw downloaded EIA-860 dataset to fill in any gaps.

    Returns:
        pd.DataFrame: the complete list of generators from PUDL with their unit IDs,
            prime movers, operating dates, and operating status.
    """
    complete_gens = load_pudl_table(
        "out_eia__yearly_generators",
        columns=[
            "report_date",
            "plant_id_eia",
            "generator_id",
            "unit_id_pudl",
            "prime_mover_code",
            "operational_status_code",
            "generator_operating_date",
            "generator_retirement_date",
            "original_planned_generator_operating_date",
            "current_planned_generator_operating_date",
        ],
    )

    # create a column that indicates the earliest year a generator reported data to EIA
    complete_gens["earliest_report_date"] = complete_gens.groupby(
        ["plant_id_eia", "generator_id"]
    )["report_date"].transform("min")

    # drop any data that was reported prior to the earliest year
    complete_gens = complete_gens[
        (complete_gens["report_date"].dt.year >= earliest_data_year)
    ]

    # for any retired gens, forward fill the most recently available unit_id_pudl to
    # the most recent available year
    complete_gens["unit_id_pudl"] = complete_gens.groupby(
        ["plant_id_eia", "generator_id"]
    )["unit_id_pudl"].ffill()

    # remove generators that retired prior to the earliest year
    complete_gens = complete_gens[
        ~(
            (complete_gens["operational_status_code"] == "RE")
            & (complete_gens["generator_retirement_date"].dt.year < earliest_data_year)
        )
    ]

    # remove generators that are proposed but not yet under construction, or cancelled
    complete_gens = complete_gens[
        ~complete_gens["operational_status_code"].isin(cancelled_status_codes)
    ]

    # only keep the most recent entry for each generator
    complete_gens = complete_gens.sort_values(
        by=["plant_id_eia", "generator_id", "report_date"], ascending=True
    ).drop_duplicates(subset=["plant_id_eia", "generator_id"], keep="last")

    # remove any generators that were under construction sometime after
    # `constants.earliest_data_year` but were cancelled or disappeared from the data
    # before `constants.earliest_validated_year``
    under_construction_status_codes = ["U", "V", "TS"]
    complete_gens = complete_gens[
        ~(
            (
                complete_gens["report_date"].dt.year
                < max(latest_validated_year, current_early_release_year)
            )
            & (
                complete_gens["operational_status_code"].isin(
                    under_construction_status_codes
                )
            )
        )
    ]

    # remove generators that have no operating or retirement date, and the last time
    # they reported data was prior to the earliest validated year. This is often
    # proposed plants that are assigned a new plant_id_eia once operational
    complete_gens = complete_gens[
        ~(
            (complete_gens["generator_operating_date"].isna())
            & (complete_gens["generator_retirement_date"].isna())
            & (complete_gens["report_date"].dt.year < earliest_validated_year)
        )
    ]

    # remove generators that have no operating or retirement date as of the latest
    # validated year and which did not have a status of testing.
    complete_gens = complete_gens[
        ~(
            (complete_gens["generator_operating_date"].isna())
            & (complete_gens["generator_retirement_date"].isna())
            & (
                complete_gens["report_date"].dt.year
                < max(latest_validated_year, current_early_release_year)
            )
            & (complete_gens["operational_status_code"] != "TS")
        )
    ]

    ####################
    # merge into complete_gens and fill missing operating dates with the EIA-860 data
    generator_data_from_eia860 = load_raw_eia860_generator_dates_and_unit_ids(
        max(latest_validated_year, current_early_release_year)
    )
    complete_gens = complete_gens.merge(
        generator_data_from_eia860,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="1:1",
    )
    complete_gens["generator_operating_date"] = complete_gens[
        "generator_operating_date"
    ].fillna(complete_gens["operating_date_eia"])
    complete_gens = complete_gens.drop(columns="operating_date_eia")

    #######################
    # update the unit_id_eia_numeric to be one higher than the highest existing
    # unit_id_pudl. If unit_id_eia_numeric is NA, the updated value should also still
    # be NA
    complete_gens["unit_id_eia_numeric"] = complete_gens[
        "unit_id_eia_numeric"
    ] + complete_gens.groupby("plant_id_eia")["unit_id_pudl"].transform("max").fillna(0)

    # fill in missing unit_id_pudl with the updated values
    complete_gens["unit_id_pudl"] = complete_gens["unit_id_pudl"].fillna(
        complete_gens["unit_id_eia_numeric"]
    )

    return complete_gens


def load_raw_eia860_plant_geographical_info(year: int) -> pd.DataFrame:
    """Loads plant geographical information from the raw EIA-860 to fill in missing
    information in the pudl data. PUDL deletes data for these fields if there are
    inconsistencies across the historical data.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: list of plants from EIA-860 with their geographical information.
    """
    # load geographic information from the raw EIA-860 file to supplement missing
    # information from pudl
    if (year == current_early_release_year) and (
        current_early_release_year != latest_validated_year
    ):
        filepath = find_downloaded_file(
            f"eia860/eia860{year}ER", f"2___Plant_Y{year}*.xlsx"
        )
        header_row = 2
    else:
        filepath = downloads_folder(f"eia860/eia860{year}/2___Plant_Y{year}.xlsx")
        header_row = 1
    plant_geographical_eia860 = pd.read_excel(
        filepath,
        header=header_row,
        usecols=[
            "Plant Code",
            "Plant Name",
            "City",
            "State",
            "County",
            "Latitude",
            "Longitude",
        ],
    ).rename(
        columns={
            "Plant Code": "plant_id_eia",
            "Plant Name": "plant_name_eia",
            "City": "city",
            "State": "state",
            "County": "county",
            "Latitude": "latitude",
            "Longitude": "longitude",
        }
    )

    for c, t in [
        ["latitude", "Float64"],
        ["longitude", "Float64"],
    ]:
        if plant_geographical_eia860[c].apply(lambda x: x == " ").sum() > 0:
            plant_geographical_eia860[c] = (
                plant_geographical_eia860[c].replace(" ", pd.NA).astype(t)
            )

    return plant_geographical_eia860


def load_raw_eia860_generator_dates_and_unit_ids(year: int) -> pd.DataFrame:
    """Loads generator operating dates and unit_id_eia codes from the raw EIA-860 to
    fill in missing dates and unit ids in the pudl data. PUDL deletes data for these
    fields if there are inconsistencies across the historical data

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: list of generators from EIA-860 with their operating dates and
            unit IDs
    """
    # load operating dates from the raw EIA-860 file to supplement missing operating
    # dates from pudl
    if (year == current_early_release_year) and (
        current_early_release_year != latest_validated_year
    ):
        filepath = find_downloaded_file(
            f"eia860/eia860{year}ER", f"3_1_Generator_Y{year}*.xlsx"
        )
        header_row = 2
        footer_row = 1
    else:
        filepath = downloads_folder(f"eia860/eia860{year}/3_1_Generator_Y{year}.xlsx")
        header_row = 1
        footer_row = 1

    generator_op_dates_eia860 = pd.read_excel(
        filepath,
        header=header_row,
        skipfooter=footer_row,
        sheet_name="Operable",
        usecols=[
            "Plant Code",
            "Generator ID",
            "Operating Month",
            "Operating Year",
            "Unit Code",
        ],
    ).rename(
        columns={
            "Plant Code": "plant_id_eia",
            "Generator ID": "generator_id",
            "Operating Month": "month",
            "Operating Year": "year",
            "Unit Code": "unit_id_eia",
        }
    )

    # create a datetime column from the month and year
    # Set " " to NA before converting to datetime
    generator_op_dates_eia860["month"] = generator_op_dates_eia860["month"].replace(
        " ", pd.NA
    )
    generator_op_dates_eia860["year"] = generator_op_dates_eia860["year"].replace(
        " ", pd.NA
    )
    generator_op_dates_eia860["operating_date_eia"] = pd.to_datetime(
        generator_op_dates_eia860[["year", "month"]].assign(Day=1)
    )
    generator_op_dates_eia860 = generator_op_dates_eia860.drop(
        columns=["month", "year"]
    )

    # load unit codes for proposed generators
    proposed_unit_ids_eia860 = (
        pd.read_excel(
            filepath,
            sheet_name="Proposed",
            header=header_row,
            usecols=["Plant Code", "Generator ID", "Unit Code"],
        )
        .dropna(subset="Unit Code")
        .rename(
            columns={
                "Plant Code": "plant_id_eia",
                "Generator ID": "generator_id",
                "Unit Code": "unit_id_eia",
            }
        )
    )

    # EIA-860 does not include Puerto Rico, so load the operating dates and unit codes
    # of Puerto Rico generators from the December EIA-860M
    # TODO: remove once PUDL keeps a yearly history of Puerto Rico generators
    puerto_rico_unit_ids = (
        load_raw_eia860m_puerto_rico_generators(year).rename(
            columns={"generator_operating_date": "operating_date_eia"}
        )[["plant_id_eia", "generator_id", "operating_date_eia", "unit_id_eia"]]
        if year >= earliest_puerto_rico_data_year
        else None
    )

    # concat the data together
    generator_data_from_eia860 = pd.concat(
        [generator_op_dates_eia860, proposed_unit_ids_eia860, puerto_rico_unit_ids],
        axis=0,
    )

    # create a numeric version of the ID, starting at 1
    generator_data_from_eia860["unit_id_eia_numeric"] = (
        generator_data_from_eia860.groupby(["plant_id_eia"])["unit_id_eia"].transform(
            lambda x: pd.factorize(x)[0] + 1
        )
    )

    # unit_id_eia_numeric of 0 represents missing unit_id_eia, so we want to replace
    # these with nan
    generator_data_from_eia860["unit_id_eia_numeric"] = generator_data_from_eia860[
        "unit_id_eia_numeric"
    ].replace(0, np.NaN)

    generator_data_from_eia860["unit_id_eia_numeric"] = pd.to_numeric(
        generator_data_from_eia860["unit_id_eia_numeric"]
    )

    return generator_data_from_eia860


def load_raw_eia860m_puerto_rico_generators(year: int) -> pd.DataFrame:
    """Loads the Puerto Rico generators reported in the December EIA-860M of a year.

    Puerto Rico generators are only reported in EIA-860M, which lists them from December
    2018 onward. PUDL only includes the most recent EIA-860M snapshot of them, so until
    PUDL keeps a history of these generators, we use the December file of each year, as
    eGRID does. Operating and planned generators are kept, along with generators that
    retired during the year.

    Args:
        year (int): a four-digit year, 2018 or later.

    Raises:
        ValueError: if `year` is before 2018.

    Returns:
        pd.DataFrame: one row per generator, using PUDL's EIA-860 column names, with
            the generator's operating status as of December and a report_date at the
            start of the year.
    """
    if year < earliest_puerto_rico_data_year:
        raise ValueError(
            f"EIA-860M does not report Puerto Rico generators for {year}, only for "
            f"{earliest_puerto_rico_data_year} onward."
        )
    filepath = downloads_folder(f"eia860m/december_generator{year}.xlsx")

    column_names = {
        "Entity ID": "utility_id_eia",
        "Entity Name": "utility_name_eia",
        "Plant ID": "plant_id_eia",
        "Plant Name": "plant_name_eia",
        "Plant State": "state",
        "County": "county",
        "Balancing Authority Code": "balancing_authority_code_eia",
        "Generator ID": "generator_id",
        "Unit Code": "unit_id_eia",
        "Nameplate Capacity (MW)": "capacity_mw",
        "Net Summer Capacity (MW)": "summer_capacity_mw",
        "Net Winter Capacity (MW)": "winter_capacity_mw",
        "Nameplate Energy Capacity (MWh)": "energy_storage_capacity_mwh",
        "DC Net Capacity (MW)": "net_capacity_mwdc",
        "Technology": "technology_description",
        "Energy Source Code": "energy_source_code_1",
        "Prime Mover Code": "prime_mover_code",
        "Status": "operational_status_code",
        "Latitude": "latitude",
        "Longitude": "longitude",
    }
    # each date is reported as a month column and a year column
    date_columns = {
        "generator_operating_date": ("Operating Month", "Operating Year"),
        "current_planned_generator_operating_date": (
            "Planned Operation Month",
            "Planned Operation Year",
        ),
        "planned_generator_retirement_date": (
            "Planned Retirement Month",
            "Planned Retirement Year",
        ),
        "generator_retirement_date": ("Retirement Month", "Retirement Year"),
    }

    tabs = []
    for sheet_name in ["Operating_PR", "Planned_PR", "Retired_PR"]:
        # the number of title rows above the column names changes between years
        header_row = (
            pd.read_excel(filepath, sheet_name=sheet_name, header=None, nrows=5)
            .eq("Entity ID")
            .any(axis=1)
            .idxmax()
        )
        tab = pd.read_excel(
            filepath,
            sheet_name=sheet_name,
            header=header_row,
            dtype={"Generator ID": str, "Unit Code": str},
            na_values=" ",
        )
        # drop the blank row and the notes below the table
        tab = tab[pd.to_numeric(tab["Plant ID"], errors="coerce").notna()]
        for date_column, (month_column, year_column) in date_columns.items():
            if year_column in tab.columns:
                tab[date_column] = pd.to_datetime(
                    pd.DataFrame(
                        {
                            "year": pd.to_numeric(tab[year_column], errors="coerce"),
                            "month": pd.to_numeric(tab[month_column], errors="coerce"),
                            "day": 1,
                        }
                    ),
                    errors="coerce",
                )
        tabs.append(
            tab.rename(columns=column_names).filter(
                items=list(column_names.values()) + list(date_columns)
            )
        )
    operating, planned, retired = tabs

    # the retired tab lists every generator that has retired, so only keep those that
    # retired during the year
    retired = retired[retired["generator_retirement_date"].dt.year == year].assign(
        operational_status_code="RE"
    )
    generators = pd.concat([operating, planned, retired], ignore_index=True)

    # statuses are reported like "(OP) Operating", so keep the code in parentheses
    generators["operational_status_code"] = generators[
        "operational_status_code"
    ].str.replace(r"^\((\w+)\).*", r"\1", regex=True)
    generators["operational_status"] = (
        generators["operational_status_code"]
        .map(
            {
                "RE": "retired",
                "OP": "existing",
                "SB": "existing",
                "OA": "existing",
                "OS": "existing",
            }
        )
        .fillna("proposed")
    )
    generators["report_date"] = pd.Timestamp(year=year, month=1, day=1)

    return apply_dtypes(generators)


def load_raw_eia923_puerto_rico_generation_fuel(year: int) -> pd.DataFrame:
    """Loads monthly generation and fuel consumption for Puerto Rico plants from the raw
    EIA-923, which reports them on a separate sheet from the other states. PUDL does not
    yet include this sheet.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: one row per plant, prime mover, energy source code, and month,
            using the column names of PUDL's generation and fuel table.
    """
    is_early_release = (year == current_early_release_year) and (
        current_early_release_year != latest_validated_year
    )
    directory = f"eia923/f923_{year}er" if is_early_release else f"eia923/f923_{year}"
    filepath = find_downloaded_file(
        directory, f"EIA923_Schedules_2_3_4_5_M_12_{year}*.xlsx"
    )

    # the number of title rows above the column names changes between years, and early
    # release files add a column of notes to the left of the table
    title_rows = pd.read_excel(
        filepath, sheet_name="Page 1 Puerto Rico", header=None, nrows=15
    )
    header_row = title_rows.eq("Plant Id").any(axis=1).idxmax()
    generation_fuel = pd.read_excel(
        filepath,
        sheet_name="Page 1 Puerto Rico",
        header=header_row,
        na_values=".",
    )
    generation_fuel.columns = generation_fuel.columns.str.replace("\n", " ")
    generation_fuel = generation_fuel.loc[:, "Plant Id":]

    # drop the rows that EIA uses to reconcile plant totals with state totals, which PUDL
    # also drops from the other states
    generation_fuel = generation_fuel[
        ~generation_fuel["Plant Id"].isin([99999, 999999])
    ]

    monthly_values = {
        "Quantity": "fuel_consumed_units",
        "Elec_Quantity": "fuel_consumed_for_electricity_units",
        "MMBtuPer_Unit": "fuel_mmbtu_per_unit",
        "Tot_MMBtu": "fuel_consumed_mmbtu",
        "Elec_MMBtu": "fuel_consumed_for_electricity_mmbtu",
        "Netgen": "net_generation_mwh",
    }
    id_columns = {
        "Plant Id": "plant_id_eia",
        "Plant Name": "plant_name_eia",
        "Operator Id": "utility_id_eia",
        "Operator Name": "utility_name_eia",
        "Reported Prime Mover": "prime_mover_code",
        "Reported Fuel Type Code": "energy_source_code",
    }
    months = pd.date_range(f"{year}-01-01", periods=12, freq="MS").strftime("%B")
    value_columns = [f"{value} {month}" for value in monthly_values for month in months]

    # a plant can report the same prime mover and fuel on more than one row (e.g.
    # Cervecera de Puerto Rico in 2020, once as CHP and once as non-CHP), so combine
    # these rows and recalculate their heat content, as PUDL does for the other states
    is_duplicated = generation_fuel.duplicated(
        subset=["Plant Id", "Reported Prime Mover", "Reported Fuel Type Code"],
        keep=False,
    )
    combined_rows = (
        generation_fuel[is_duplicated]
        .groupby(list(id_columns), dropna=False)[value_columns]
        .sum(min_count=1)
        .reset_index()
    )
    for month in months:
        combined_rows[f"MMBtuPer_Unit {month}"] = (
            combined_rows[f"Tot_MMBtu {month}"] / combined_rows[f"Quantity {month}"]
        ).replace([np.inf, -np.inf], np.nan)
    generation_fuel = pd.concat(
        [generation_fuel[~is_duplicated], combined_rows], ignore_index=True
    )

    # each value is reported in one column per month, so convert these to one row per
    # month
    generation_fuel = generation_fuel.set_index(list(id_columns))[value_columns]
    generation_fuel.columns = pd.MultiIndex.from_tuples(
        [column.rsplit(" ", 1) for column in value_columns], names=[None, "month"]
    )
    generation_fuel = (
        generation_fuel.stack("month", future_stack=True)
        .reset_index()
        .rename(columns=id_columns | monthly_values)
    )
    generation_fuel["report_date"] = pd.to_datetime(
        generation_fuel.pop("month") + f" {year}", format="%B %Y"
    )
    # drop months with no values reported, as PUDL does for the other states
    generation_fuel = generation_fuel.dropna(
        subset=list(monthly_values.values()), how="all"
    )

    # add the columns that PUDL derives from the EIA-923 data
    fuel_types = load_pudl_table(
        "core_eia__codes_energy_sources", columns=["code", "fuel_type_code_pudl"]
    )
    generation_fuel["fuel_type_code_pudl"] = generation_fuel["energy_source_code"].map(
        fuel_types.set_index("code")["fuel_type_code_pudl"]
    )
    generation_fuel["data_maturity"] = "provisional" if is_early_release else "final"

    return apply_dtypes(generation_fuel)


def update_epa_to_eia_map(cems_df: pd.DataFrame, year: int) -> pd.DataFrame:
    """Updates the plant_id_eia column in the CEMS data frame loaded from pudl based
    on the manual epa_eia_crosswalk_manual table.

    This is meant to be applied on a single year of CEMS data at a time, since EPA-EIA
    mappings can change from year to year

    Args:
        cems_df (pd.DataFrame): input CEMS data frame.
        year (int): used to filter to the correct mapping for the current year

    Returns:
        pd.DataFrame: CEMS data frame with updated plant_id_eia column.
    """
    # load the manual table
    manual_plant_map = pd.read_csv(
        reference_table_folder("epa_eia_crosswalk_manual.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["notes"])

    # if no start_year or end_year is specified, fill with numbers to help with filtering
    manual_plant_map["start_year"] = manual_plant_map["start_year"].fillna(0)
    manual_plant_map["end_year"] = manual_plant_map["end_year"].fillna(9999)

    # only keep mappings that are valid for the current year
    manual_plant_map = manual_plant_map[
        (year >= manual_plant_map["start_year"])
        & (year <= manual_plant_map["end_year"])
    ].drop(columns=["start_year", "end_year"])

    manual_plant_map = manual_plant_map[
        ["plant_id_epa", "emissions_unit_id_epa", "plant_id_eia"]
    ].drop_duplicates()

    # merge into the cems data
    cems_df = cems_df.merge(
        manual_plant_map,
        how="left",
        on=["plant_id_epa", "emissions_unit_id_epa"],
        suffixes=(None, "_manual"),
        validate="m:1",
    )

    # update the eia plant ids
    cems_df.update({"plant_id_eia": cems_df["plant_id_eia_manual"]})

    # drop the intermediate column
    cems_df = cems_df.drop(columns=["plant_id_eia_manual"])

    return cems_df


def add_report_date(
    df: pd.DataFrame, plant_timezone: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Add a report_date column to the a dataframe based on the plant's local timezone.

    Args:
        df (pd.DataFrame): data frame containing plant_id_eia and datetime_utc columns.
        plant_timezone: (pd.DataFrame | None), by default, if None, this function will
            load timezone data from the pudl plants entity table. However, in order to
            use this function where pudl cannot be read from s3, this provides the
            option to pass in your own plant timezone dataframe, for example loaded from
            outputs/plant_static_attributes

    Returns:
        pd.DataFrame: original data frame with a report_date column added.
    """
    if plant_timezone is None:
        plant_timezone = load_pudl_table(
            "core_eia__entity_plants", columns=["plant_id_eia", "timezone"]
        )

    # get timezone
    df = df.merge(
        plant_timezone,
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )

    # create a datetimeindex from the datetime_utc column
    datetime_utc = pd.DatetimeIndex(df["datetime_utc"])

    # create blank column to hold local datetimes
    df["report_date"] = pd.to_datetime(np.NaN)

    # get list of unique timezones
    timezones = list(df["timezone"].unique())

    # convert UTC to the local timezone
    for tz in timezones:
        tz_mask = df["timezone"] == tz  # find all rows where the tz matches
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Converting to PeriodArray/Index representation will drop timezone information.",
            )
            df.loc[tz_mask, "report_date"] = (
                datetime_utc[tz_mask]
                .tz_convert(tz)  # convert to local time
                .to_series(index=df[tz_mask].index)  # convert to a series
                .dt.to_period("M")
                .dt.to_timestamp()  # convert to a YYYY-MM-01 stamp
            )

    df["report_date"] = pd.to_datetime(df["report_date"]).astype("datetime64[s]")

    # drop the operating_datetime_local column
    df = df.drop(columns=["timezone"])

    return df


def load_ghg_emission_factors() -> pd.DataFrame:
    """Read in the table of emissions factors and convert to lb/mmbtu.

    Returns:
        pd.DataFrame: GHG emission factors table.
    """

    efs = pd.read_csv(
        reference_table_folder("emission_factors_for_co2_ch4_n2o.csv"),
        dtype=get_dtypes(),
    )

    # convert co2 mass in short tons to lb
    efs["co2_tons_per_mmbtu"] = (
        efs["co2_tons_per_mmbtu"] * ConversionFactors.short_ton_to_lbs
    )

    # rename the columns
    efs = efs.rename(columns={"co2_tons_per_mmbtu": "co2_lb_per_mmbtu"})

    return efs


def load_nox_emission_factors() -> pd.DataFrame:
    """Read in the NOx emission factors.

    Returns:
        pd.DataFrame: NOx emission factors table.
    """
    emission_factors = pd.read_csv(
        reference_table_folder("emission_factors_for_nox.csv"),
        dtype=get_dtypes(),
    )

    return emission_factors


def load_so2_emission_factors() -> pd.DataFrame:
    """Read in the SO2 emission factors.

    Returns:
        pd.DataFrame: SO2 emission factors table.
    """
    df = pd.read_csv(
        reference_table_folder("emission_factors_for_so2.csv"),
        dtype=get_dtypes(),
    )

    return df


def load_pudl_table_from_local(
    table_name: str,
    columns: list[str] | None = None,
    datetime_column: str | None = None,
    start_datetime: pd.Timestamp | None = None,
    end_datetime: pd.Timestamp | None = None,
) -> pd.DataFrame:
    """Load a table from the PUDL SQLite database.

    Args:
        table_name (str): the name of the table to load.
        columns (list[str] | None, optional): the columns to select. Defaults to None.
        datetime_column (str, optional): the name of the datetime column to filter on.
            Defaults to None.
        start_datetime (pd.Timestamp, optional): the start datetime for filtering.
            Defaults to None.
        end_datetime (pd.Timestamp, optional): the end datetime for filtering.
            Defaults to None.

    Returns:
        pd.DataFrame: The loaded table.
    """
    columns_to_select = "*" if columns is None else ", ".join(columns)

    # Validate datetime filtering parameters
    if (datetime_column is None) != (start_datetime is None and end_datetime is None):
        raise ValueError(
            "datetime_column must be specified if and only if start_datetime and "
            "end_datetime are provided"
        )

    # Check if table exists in the database
    inspector = sa.inspect(PUDL_ENGINE)

    if datetime_column is None:
        # No filtering on datetime
        table = (
            pd.read_sql(
                f"SELECT {columns_to_select} FROM {table_name}",
                PUDL_ENGINE,
            )
            if inspector.has_table(table_name)
            else pd.read_parquet(
                pudl_folder(f"{table_name}.parquet"),
                columns=columns,
            )
        )
    else:
        # Filter on datetime
        start = start_datetime.strftime("%Y-%m-%d")
        end = end_datetime.strftime("%Y-%m-%d")
        table = (
            pd.read_sql(
                f"SELECT {columns_to_select} FROM {table_name} WHERE \
                {datetime_column} >= '{start}' AND {datetime_column} <= '{end}'",
                PUDL_ENGINE,
            )
            if inspector.has_table(table_name)
            else pd.read_parquet(
                pudl_folder(f"{table_name}.parquet"),
                filters=[
                    (datetime_column, ">=", start_datetime),
                    (datetime_column, "<=", end_datetime),
                ],
                columns=columns,
            )
        )
    return table


def load_pudl_table_from_s3(
    table_name: str,
    datetime_column: str | None = None,
    columns: list[str] | None = None,
    start_datetime: pd.Timestamp | None = None,
    end_datetime: pd.Timestamp | None = None,
) -> pd.DataFrame:
    """Load a table from a Parquet file.

    Logs the PUDL release being used for transparency, since OGE only specifies whether
    to use "stable" or "nightly" builds, and the specific release version can change
    over time.

    Args:
        table_name (str): the name of the table to load.
        datetime_column (str, optional): the name of the datetime column to filter on.
            Defaults to None.
        columns (list[str] | None, optional): the columns to select. Defaults to None.
        start_datetime (pd.Timestamp, optional): the start datetime for filtering.
            Defaults to None.
        end_datetime (pd.Timestamp, optional): the end datetime for filtering.
            Defaults to None.

    Returns:
        pd.DataFrame: The loaded table.
    """

    # Validate datetime filtering parameters
    if (datetime_column is None) != (start_datetime is None and end_datetime is None):
        raise ValueError(
            "datetime_column must be specified if and only if start_datetime and "
            "end_datetime are provided"
        )
    pudl_release = get_pudl_release(get_pudl_build_version(), get_pudl_data_store())
    logger.info(f"Loading {table_name} from PUDL {pudl_release}")
    if datetime_column is None:
        # No filtering on datetime
        table = pd.read_parquet(
            pudl_folder(f"{table_name}.parquet"),
            columns=columns,
        )
    else:
        # Filter on datetime
        filters = [
            (datetime_column, ">=", start_datetime),
            (datetime_column, "<=", end_datetime),
        ]
        table = pd.read_parquet(
            pudl_folder(f"{table_name}.parquet"),
            filters=filters,
            columns=columns,
        )
    return table


def load_pudl_table(
    table_name: str,
    year: int = None,
    columns: list[str] = None,
    end_year: int = None,
    dt: pd.Timestamp = None,
    end_dt: pd.Timestamp = None,
) -> pd.DataFrame:
    """Loads a table from the downloaded SQLite database or s3 parquet files.

    Args:
        table_name (str): one of the options specified in the data dictionary:
            https://catalystcoop-pudl.readthedocs.io/en/latest/data_dictionaries/pudl_db.html
        year (int, optional): a four-digit year.
            If specified, but not `end_year`, only a single year will be loaded.
            If not specified, all years will be loaded.
            if both are specified, all years in that range inclusive will be loaded.
            Defaults to None.
        columns (list[str], optional): Columns to return. If not specified, all columns
            will be returned. Defaults to None.
        end_year (int, optional): a four-digit year >= `year`. Defaults to None.
        dt (pd.Timestamp, optional): a specific datetime to filter on. Defaults to None.
        end_dt (pd.Timestamp, optional): a specific end datetime to filter on. Defaults
            to None.

    Returns:
        pd.DataFrame: the required table.
    """

    # Validate that year/end_year and dt/end_dt are not both specified
    if (year is not None or end_year is not None) and (
        dt is not None or end_dt is not None
    ):
        raise ValueError(
            "Cannot specify both year-based ('year', 'end_year') and "
            "datetime-based ('dt', 'end_dt') parameters. They are mutually exclusive."
        )

    # Validate parameter dependencies
    if (end_year is not None and year is None) or (end_dt is not None and dt is None):
        raise ValueError(
            "'end_year' cannot be specified without 'year' and "
            "'end_dt' cannot be specified without 'dt'."
        )

    datetime_column = None
    start_datetime = None
    end_datetime = None

    if year is not None:
        # filter on report date
        datetime_column = "report_date"
        start_datetime = pd.Timestamp(year, 1, 1)
        end_datetime = (
            pd.Timestamp(end_year, 12, 1)
            if end_year is not None
            else pd.Timestamp(year, 12, 1)
        )
    if dt is not None:
        datetime_column = "datetime_utc"
        start_datetime = dt
        end_datetime = end_dt if end_dt is not None else dt

    table = (
        load_pudl_table_from_local(
            table_name,
            datetime_column=datetime_column,
            columns=columns,
            start_datetime=start_datetime,
            end_datetime=end_datetime,
        )
        if PUDL_ENGINE
        else load_pudl_table_from_s3(
            table_name,
            datetime_column=datetime_column,
            columns=columns,
            start_datetime=start_datetime,
            end_datetime=end_datetime,
        )
    )

    # Puerto Rico
    # Puerto Rico generators are only reported in EIA-860M, and PUDL only includes the most
    # recent EIA-860M snapshot of them, with a report_date in the latest year. Until PUDL
    # keeps a yearly history of these generators, we replace PUDL's snapshot with the
    # generators in the December EIA-860M of each year, as eGRID does. PUDL also does not
    # yet include the Puerto Rico generation and fuel data, which EIA-923 reports on a
    # separate sheet, so we load it from the raw EIA-923.
    # TODO: remove once PUDL keeps a yearly history of Puerto Rico generators and
    # includes the EIA-923 Puerto Rico data
    puerto_rico_loaders = {
        "core_eia860__scd_generators": load_puerto_rico_generators,
        "out_eia__yearly_generators": load_puerto_rico_generators,
        "out_eia923__monthly_generation_fuel_combined": (
            load_raw_eia923_puerto_rico_generation_fuel
        ),
        "out_eia923__generation_fuel_combined": (
            load_raw_eia923_puerto_rico_generation_fuel
        ),
    }
    if table_name in puerto_rico_loaders:
        table = replace_puerto_rico_data(
            table, puerto_rico_loaders[table_name], year, end_year, columns
        )

    table = apply_dtypes(table)

    if table.empty:
        logger.warning(f"{table_name} is empty")

    return table


def replace_puerto_rico_data(
    table: pd.DataFrame,
    load_puerto_rico_year: Callable[[int], pd.DataFrame],
    year: int | None,
    end_year: int | None,
    columns: list[str] | None,
) -> pd.DataFrame:
    """Replaces the Puerto Rico data in a PUDL table with data loaded for each year
    that was loaded. Any rows that PUDL has for these plants are dropped.

    Args:
        table (pd.DataFrame): a PUDL table, which must include plant_id_eia.
        load_puerto_rico_year (Callable[[int], pd.DataFrame]): loads the Puerto Rico
            data for one year, in the format of the PUDL table.
        year (int | None): the first year loaded, or None if all years were loaded.
        end_year (int | None): the last year loaded, or None if only `year` was
            loaded.
        columns (list[str] | None): the columns loaded, or None if all were loaded.

    Returns:
        pd.DataFrame: the table, with Puerto Rico data for each year loaded from 2018
            onward.
    """
    latest_year = max(latest_validated_year, current_early_release_year)
    first_year = max(year or earliest_data_year, earliest_puerto_rico_data_year)
    last_year = min(end_year or year or latest_year, latest_year)

    if first_year > last_year:
        return table
    puerto_rico_data = pd.concat(
        [
            load_puerto_rico_year(puerto_rico_year)
            for puerto_rico_year in range(first_year, last_year + 1)
        ],
        ignore_index=True,
    )
    # drop any PUDL data for these plants, such as PUDL's EIA-860M snapshot of Puerto
    # Rico generators, which is only loaded with the latest year
    table = table[~table["plant_id_eia"].isin(puerto_rico_data["plant_id_eia"])]
    if columns is not None:
        puerto_rico_data = puerto_rico_data.reindex(columns=columns)
    return pd.concat([table, puerto_rico_data], ignore_index=True)


@functools.cache
def load_puerto_rico_generators(year: int) -> pd.DataFrame:
    """Loads the Puerto Rico generators for a year in the format of PUDL's generator
    tables.

    Generators come from the December EIA-860M of the year. Plants that report
    generation to EIA-923 in the year but are missing from that file get their
    generators from PUDL's EIA-860M snapshot instead, keeping those that were operating
    during the year.

    The result is cached for each year, since the generator tables are loaded many times
    during a pipeline run. Callers should not modify the returned table in place.

    Args:
        year (int): a four-digit year, 2018 or later.

    Returns:
        pd.DataFrame: one row per Puerto Rico generator.
    """
    generators = load_raw_eia860m_puerto_rico_generators(year)

    # Some plants report generation to EIA-923 before their generators appear in
    # EIA-860M. For these, use the generators in PUDL's EIA-860M snapshot, which is the
    # most recent one available.
    plants_with_generation = load_raw_eia923_puerto_rico_generation_fuel(year)[
        "plant_id_eia"
    ].unique()
    missing_plants = set(plants_with_generation) - set(generators["plant_id_eia"])
    if missing_plants:
        # read the table directly, since loading it with load_pudl_table would replace
        # its Puerto Rico generators by calling this function again
        snapshot_columns = generators.columns.drop(
            ["unit_id_eia", "report_date"]
        ).tolist()
        snapshot = apply_dtypes(
            load_pudl_table_from_local(
                "out_eia__yearly_generators", columns=snapshot_columns
            )
            if PUDL_ENGINE
            else load_pudl_table_from_s3(
                "out_eia__yearly_generators", columns=snapshot_columns
            )
        )
        snapshot = snapshot[
            snapshot["plant_id_eia"].isin(missing_plants)
            & (snapshot["generator_operating_date"].dt.year <= year)
            & ~(snapshot["generator_retirement_date"].dt.year < year)
        ].copy()
        # generators that retired after the year were still operating during it
        retired_later = snapshot["generator_retirement_date"].dt.year > year
        snapshot.loc[retired_later, "generator_retirement_date"] = pd.NaT
        snapshot.loc[retired_later, "operational_status_code"] = "OP"
        snapshot.loc[retired_later, "operational_status"] = "existing"
        snapshot["report_date"] = pd.Timestamp(year=year, month=1, day=1)
        logger.info(
            f"Using PUDL's EIA-860M snapshot for {snapshot['plant_id_eia'].nunique()} "
            f"Puerto Rico plants that report generation in {year} but are missing from "
            "that year's December EIA-860M"
        )
        generators = pd.concat([generators, snapshot], ignore_index=True)

    # add the columns that PUDL derives from the EIA-860 data
    fuel_types = load_pudl_table(
        "core_eia__codes_energy_sources", columns=["code", "fuel_type_code_pudl"]
    )
    generators["fuel_type_code_pudl"] = generators["energy_source_code_1"].map(
        fuel_types.set_index("code")["fuel_type_code_pudl"]
    )
    generators["fuel_type_count"] = (
        generators["energy_source_code_1"].notna().astype(int)
    )
    generators["data_maturity"] = "monthly_update"

    return apply_dtypes(generators)


def load_epa_eia_crosswalk_from_raw() -> pd.DataFrame:
    """Read in the EPA-EIA Crosswalk table downloaded from the EPA website.

    This is only used in data_cleaning.remove_non_grid_connected_plants() to access the
    MATCH_TYPE_GEN column, which identifies the units that CAMD identified as NGC.

    Returns:
        pd.DataFrame: EPA-EIA crosswalk table.
    """

    crosswalk = pd.read_csv(
        downloads_folder("epa/epa_eia_crosswalk.csv"),
        usecols=["CAMD_PLANT_ID", "EIA_PLANT_ID", "CAMD_UNIT_ID", "MATCH_TYPE_GEN"],
    )

    # remove leading zeros from the emissions_unit_id_epa
    crosswalk["CAMD_UNIT_ID"] = crosswalk["CAMD_UNIT_ID"].str.lstrip("0")

    # some eia plant ids are missing. Let us assume that the EIA and EPA plant ids
    # match in this case
    crosswalk["EIA_PLANT_ID"] = crosswalk["EIA_PLANT_ID"].fillna(
        crosswalk["CAMD_PLANT_ID"]
    )
    # change the id to an int now that no plant ids are missing
    crosswalk["EIA_PLANT_ID"] = crosswalk["EIA_PLANT_ID"].astype(int)

    # rename the columns
    crosswalk = crosswalk.rename(
        columns={
            "CAMD_PLANT_ID": "plant_id_epa",
            "EIA_PLANT_ID": "plant_id_eia",
            "CAMD_UNIT_ID": "emissions_unit_id_epa",
            "MATCH_TYPE_GEN": "epa_match_type",
        }
    )

    return crosswalk


def load_epa_eia_crosswalk(year: int) -> pd.DataFrame:
    """Read in the manual EPA-EIA Crosswalk table.

    This is currently only used for subplant identification

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: manual EIA-EPA crosswalk table.
    """

    crosswalk = load_pudl_table("core_epa__assn_eia_epacamd")
    crosswalk["emissions_unit_id_epa"] = crosswalk["emissions_unit_id_epa"].str.lstrip(
        "0"
    )

    # load manually inputted data
    crosswalk_manual = pd.read_csv(
        reference_table_folder("epa_eia_crosswalk_manual.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["notes"])

    # concat this data with the main table
    crosswalk = pd.concat(
        [crosswalk, crosswalk_manual],
        axis=0,
    )

    # load eGRID plant mapping
    egrid_plant_map = pd.read_csv(
        reference_table_folder("eGRID_crosswalk_of_EIA_ID_to_EPA_ID.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["plant_id_egrid"])

    # concat this data with the main table
    crosswalk = pd.concat(
        [crosswalk, egrid_plant_map],
        axis=0,
    )

    # drop duplicate crosswalks
    # pudl now includes crosswalks related to multiple report years
    # keep the newest mapping
    crosswalk = crosswalk.drop_duplicates(
        subset=[
            "plant_id_epa",
            "emissions_unit_id_epa",
            "generator_id_epa",
            "plant_id_eia",
            "boiler_id",
            "generator_id",
        ],
        keep="last",
    )

    # load EIA-860 data
    gen_ids = load_pudl_table(
        "core_eia860__scd_generators", year, columns=["plant_id_eia", "generator_id"]
    )

    # merge in any plants that are missing from the EPA crosswalk but appear in EIA-860
    crosswalk = crosswalk.merge(
        gen_ids, how="outer", on=["plant_id_eia", "generator_id"], validate="m:1"
    )
    crosswalk["plant_id_epa"] = crosswalk["plant_id_epa"].fillna(
        crosswalk["plant_id_eia"]
    )

    return crosswalk


def load_campd_facility(year: int) -> pd.DataFrame:
    """Loads the unit attributes from the EPA's Clean Air Markets Program Data (CAMPD)
    facility file. These are used to fill unit attributes that are missing from the EIA
    data, and to report the attributes of CEMS units in the eGRID format.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table with one row per CAMPD unit, giving its EIA plant ID, unit
            type, boiler firing type and wet/dry bottom type, fuel types and the EIA
            energy source code of its primary fuel, source category, associated
            generators, program codes, operating status, emission controls, and
            commercial operation date.
    """
    facility = pd.read_csv(
        downloads_folder(f"epa/facility-{year}.csv"),
        usecols=[
            "Facility ID",
            "Unit ID",
            "Facility Name",
            "State",
            "Source Category",
            "Unit Type",
            "Primary Fuel Type",
            "Secondary Fuel Type",
            "Associated Generators & Nameplate Capacity (MWe)",
            "Program Code",
            "Operating Status",
            "SO2 Controls",
            "NOx Controls",
            "Hg Controls",
            "Commercial Operation Date",
        ],
        dtype={"Unit ID": "str"},
    ).rename(
        columns={
            "Facility ID": "plant_id_epa",
            "Unit ID": "emissions_unit_id_epa",
            "Facility Name": "plant_name_epa",
            "State": "state",
            "Source Category": "source_category_epa",
            "Unit Type": "unit_type_epa",
            "Primary Fuel Type": "primary_fuel_type_epa",
            "Secondary Fuel Type": "secondary_fuel_type_epa",
            "Associated Generators & Nameplate Capacity (MWe)": (
                "associated_generators_and_capacity_epa"
            ),
            "Program Code": "program_code_epa",
            "Operating Status": "operating_status_epa",
            "SO2 Controls": "so2_controls_epa",
            "NOx Controls": "nox_controls_epa",
            "Hg Controls": "hg_controls_epa",
            "Commercial Operation Date": "commercial_operation_date_epa",
        }
    )

    # remove leading zeros from the emissions_unit_id_epa to match the CEMS data
    facility["emissions_unit_id_epa"] = facility["emissions_unit_id_epa"].str.lstrip(
        "0"
    )

    # units that start or stop operating during the year have the date in parentheses,
    # and units that change type list each type, e.g. "Integrated gasification combined
    # cycle (Ended Jun 07, 2025), Combustion turbine (Started Jun 07, 2025)". Keep the
    # type in use at the end of the year.
    facility["unit_type_epa"] = (
        facility["unit_type_epa"]
        .str.replace(r"\s*\([^)]*\)", "", regex=True)
        .str.split(",")
        .str[-1]
        .str.strip()
    )

    # boiler firing type and wet/dry bottom type of each unit type that is a boiler,
    # using the categories in the emission factor tables. CAMPD does not report the
    # bottom type of some boiler types.
    boiler_design = {
        "Arch-fired boiler": ("vertical_firing", pd.NA),
        "Bubbling fluidized bed boiler": ("fluidized_bed_firing", pd.NA),
        "Cell burner boiler": ("cell_burner", pd.NA),
        "Circulating fluidized bed boiler": ("fluidized_bed_firing", pd.NA),
        "Cyclone boiler": ("cyclone_firing", pd.NA),
        "Dry bottom turbo-fired boiler": ("other", "dry"),
        "Dry bottom vertically-fired boiler": ("vertical_firing", "dry"),
        "Dry bottom wall-fired boiler": ("wall_fired", "dry"),
        "Other boiler": ("other", pd.NA),
        "Pressurized fluidized bed boiler": ("fluidized_bed_firing", pd.NA),
        "Stoker": ("stoker", pd.NA),
        "Tangentially-fired": ("tangential_firing", pd.NA),
        "Wet bottom turbo-fired boiler": ("other", "wet"),
        "Wet bottom vertically-fired boiler": ("vertical_firing", "wet"),
        "Wet bottom wall-fired boiler": ("wall_fired", "wet"),
    }
    unit_boiler_design = facility["unit_type_epa"].map(boiler_design)
    facility["boiler_firing_type"] = unit_boiler_design.str[0]
    facility["wet_dry_bottom"] = unit_boiler_design.str[1]

    # EIA energy source code of the primary fuel. Some units list several primary
    # fuels, in which case the first one is used.
    camd_to_eia_fuel_type = {
        "Pipeline Natural Gas": "NG",
        "Coal": "SUB",  # assume that generic coal is subbituminous to be conservative
        "Coal Refuse": "WC",
        "Residual Oil": "RFO",
        "Other Oil": "WO",
        "Diesel Oil": "DFO",
        "Natural Gas": "NG",
        "Wood": "WDS",
        "Process Gas": "PRG",
        "Other Gas": "OG",
        "Petroleum Coke": "PC",
        "Other Solid Fuel": "OBS",
        "Tire Derived Fuel": "TDF",
    }
    facility["energy_source_code_epa"] = (
        facility["primary_fuel_type_epa"]
        .str.split(",")
        .str[0]
        .map(camd_to_eia_fuel_type)
    )

    # add the EIA plant ID, which is the same as the EPA plant ID unless it is mapped
    # differently in the manual crosswalk
    facility["plant_id_eia"] = facility["plant_id_epa"]
    facility = update_epa_to_eia_map(apply_dtypes(facility), year)

    return facility


def load_gross_to_net_data(
    level: str,
    conversion_type: str,
    threshold_column: str,
    lower_threshold: float,
    upper_threshold: float,
    year: int,
) -> pd.DataFrame:
    """Loads gross-to-net generation conversion factors calculated in
    `gross_to_net_generation`.

    Args:
        level (str): aggregation level, either 'plant' or 'subplant'
        conversion_type (str): which data to load, either 'regression' or 'ratio'.
        threshold_column (str): name of the column to be used to filter values.
        lower_threshold (float): the value below which in threshold_column the
            conversion factors should be removed
        upper_threshold (float): the value above which in threshold_column the
            conversion factors should be removed
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table containing relevant keys and conversion factors.

    Example:
        If you wanted to load subplant regression results with an adjusted r2 greater
        than 0.9, you would specify:

        >>> load_gross_to_net_data('plant', 'regression', 'rsqaured_adjusted', None)
    """
    gtn_data = pd.read_csv(
        outputs_folder(f"gross_to_net/{level}_gross_to_net_{conversion_type}.csv.zip"),
        dtype=get_dtypes(),
    )

    # filter the data based on the upper and lower thresholds, if specified
    if lower_threshold is not None:
        gtn_data = gtn_data[gtn_data[threshold_column] >= lower_threshold]

    if upper_threshold is not None:
        gtn_data = gtn_data[gtn_data[threshold_column] <= upper_threshold]

    # if loading regression data, add a count of units in each subplant to the
    # regression results
    if conversion_type == "regression":
        subplant_crosswalk = pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv.zip"),
            dtype=get_dtypes(),
        )
        subplant_crosswalk = subplant_crosswalk[
            ["plant_id_eia", "emissions_unit_id_epa", "subplant_id"]
        ].drop_duplicates()

        if level == "plant":
            groupby_columns = ["plant_id_eia"]
        elif level == "subplant":
            groupby_columns = ["plant_id_eia", "subplant_id"]
        subplant_crosswalk = (
            subplant_crosswalk.groupby(groupby_columns, dropna=False)
            .count()
            .reset_index()
            .rename(columns={"emissions_unit_id_epa": f"units_in_{level}"})
        )

        gtn_data = gtn_data.merge(
            subplant_crosswalk, how="left", on=groupby_columns, validate="many_to_one"
        )

        # divide the intercept by the number of units in each subplant to evenly
        # distribute this to each unit
        gtn_data["intercept"] = gtn_data["intercept"] / gtn_data[f"units_in_{level}"]

    # make sure the report date column is a datetime if loading ratios
    if conversion_type == "ratio":
        gtn_data["report_date"] = pd.to_datetime(gtn_data["report_date"]).astype(
            "datetime64[s]"
        )

    return gtn_data


def load_ipcc_gwp() -> pd.DataFrame:
    """Loads table containing global warming potential (GWP) values for CO2, CH4,
    and N2O.


    Returns:
        pd.DataFrame: table containing the global warming potential for GHG.
    """
    return pd.read_csv(
        reference_table_folder("ipcc_gwp.csv"),
        dtype=get_dtypes(),
    )


def load_raw_eia930_data(year: int, description: str) -> pd.DataFrame:
    """Loads raw balance or interchange EIA-930 file for the full specified year.

    Args:
        year (int): a four digit year.
        description (str): either 'INTERCHANGE' or 'BALANCE'.

    Returns:
        pd.DataFrame: the eia930 table.
    """
    eia_930 = pd.concat(
        [
            pd.read_csv(
                downloads_folder(f"eia930/EIA930_{description}_{year}_Jan_Jun.csv"),
                thousands=",",
                parse_dates=["UTC Time at End of Hour"],
            ),
            pd.read_csv(
                downloads_folder(f"eia930/EIA930_{description}_{year}_Jul_Dec.csv"),
                thousands=",",
                parse_dates=["UTC Time at End of Hour"],
            ),
        ]
    )

    # convert from end of hour timestamp to beginning of hour timestamp
    eia_930["UTC Time at End of Hour"] = eia_930[
        "UTC Time at End of Hour"
    ] - pd.Timedelta(hours=1)

    # localize the timezone as UTC time and rename the column
    eia_930["UTC Time at End of Hour"] = eia_930[
        "UTC Time at End of Hour"
    ].dt.tz_localize("UTC")
    eia_930 = eia_930.rename(
        columns={"UTC Time at End of Hour": "operating_datetime_utc"}
    )

    # Make sure that the columns use a consistent naming scheme!
    # Defends against EIA suddenly adding underscores (which they have done before).
    eia_930.columns = eia_930.columns.str.replace("_", " ")

    return eia_930


def load_ba_reference() -> pd.DataFrame:
    """Loads the balancing authority information table.

    Returns:
        pd.DataFrame: the BA table.
    """
    return pd.read_csv(
        reference_table_folder("ba_reference.csv"),
        dtype=get_dtypes(),
        parse_dates=["activation_date", "retirement_date"],
    )


def load_diba_data(year: int) -> pd.DataFrame:
    """Loads information on the directly interconnected balancing authorities (DIBAs).

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: the table with DIBAs information
    """
    # load information about directly interconnected balancing authorities (DIBAs)
    dibas = load_pudl_table(
        "core_eia930__hourly_interchange", dt=pd.Timestamp(year, 1, 1)
    )

    # Set datetime_utc to the beginning of the hour and assign UTC timezone
    dibas["datetime_utc"] = dibas["datetime_utc"] - pd.Timedelta(hours=1)

    # load ba information to get regions, dropping Canada and Mexico.
    ba_info = (
        (
            load_pudl_table(
                "core_eia__codes_balancing_authorities",
                columns=["code", "balancing_authority_region_code_eia"],
            )
            .dropna()
            .query("balancing_authority_region_code_eia not in ['CAN', 'MEX']")
        )
        .set_index("code")["balancing_authority_region_code_eia"]
        .to_dict()
    )

    # add region information
    dibas["ba_region"] = dibas["balancing_authority_code_eia"].map(ba_info)
    dibas["diba_region"] = dibas["balancing_authority_code_adjacent_eia"].map(ba_info)

    # filter to relevant columns and drop duplicates
    dibas = dibas[
        [
            "balancing_authority_code_eia",
            "balancing_authority_code_adjacent_eia",
            "ba_region",
            "diba_region",
        ]
    ].drop_duplicates()
    dibas = dibas.rename(
        columns={
            "balancing_authority_code_eia": "ba_code",
            "balancing_authority_code_adjacent_eia": "diba_code",
        }
    )

    # add information about the local timezone of each ba/diba
    ba_tz = load_ba_reference()[["ba_code", "timezone_local"]]
    dibas = dibas.merge(ba_tz, how="left", on="ba_code", validate="m:1")
    dibas = dibas.merge(
        ba_tz,
        how="left",
        left_on="diba_code",
        right_on="ba_code",
        suffixes=(None, "_diba"),
        validate="m:1",
    ).drop(columns="ba_code_diba")

    return dibas


def ba_timezone(ba: str, type: str) -> str:
    """Retrieves the timezone for a single balancing area.

    Args:
        ba (str): the BA code.
        type (str): either 'reporting_eia930' or 'local'.
            'reporting_eia930' will return the timezone used by the BA when reporting
            to EIA-930.
            'local' will return the actual local timezone.

    Raises:
        UserWarning: if BA does not have a timezone in the BA reference table.

    Returns:
        str: the timezone of the BA.
    """

    tz = pd.read_csv(
        reference_table_folder("ba_reference.csv"),
        usecols=["ba_code", f"timezone_{type}"],
    )
    tz = tz.loc[tz["ba_code"] == ba, f"timezone_{type}"]

    if len(tz) == 0:
        raise UserWarning(
            f"The BA {ba} does not have a timezone specified in reference_tables/ba_reference.csv. Please add."
        )
    else:
        tz = tz.item()

    return tz


def load_emissions_controls_eia923(year: int) -> pd.DataFrame:
    """Loads emission controls information from EIA-923 for a specified year.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: the emission controls table. Note that the data frame will
    """
    emissions_controls_eia923_names = [
        "report_date",
        "plant_id_eia",
        "equipment_tech_description",
        "particulate_control_id_eia",
        "so2_control_id_eia",
        "nox_control_id_eia",
        "mercury_control_id_eia",
        "operational_status",
        "hours_in_service",  # not yet in pudl
        "annual_nox_emission_rate_lb_per_mmbtu",
        "ozone_season_nox_emission_rate_lb_per_mmbtu",  # not yet in pudl
        "particulate_emission_rate_lb_per_mmbtu",
        "particulate_removal_efficiency_annual",
        "particulate_removal_efficiency_at_full_load",
        "particulate_test_date",
        "so2_removal_efficiency_annual",
        "so2_removal_efficiency_at_full_load",
        "so2_test_date",
        "fgd_sorbent_consumption_1000_tons",
        "fgd_electricity_consumption_mwh",
        "mercury_removal_efficiency",
        "mercury_emission_rate_lb_per_trillion_btu",
        "acid_gas_removal_efficiency",
    ]

    # For 2012-2015 and earlier, mercury emission rate is not reported in EIA923, so we
    # need to remove that column to avoid an error.
    if year <= 2015:
        emissions_controls_eia923_names.remove(
            "mercury_emission_rate_lb_per_trillion_btu"
        )

    if year >= 2012:
        if year <= latest_validated_year:
            # Handle filename changes across years.
            schedule_8_filename = {
                2012: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final_Revision.xlsx"
                ),
                2013: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_PartsA-D_EnvData_2013_Final_Revision.xlsx"
                ),
                2014: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final_Revision.xlsx"
                ),
                2015: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final_Revision.xlsx"
                ),
                2016: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final_Revision.xlsx"
                ),
                2017: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Envir_Infor_{year}_Final.xlsx"
                ),
                2018: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final.xlsx"
                ),
                2019: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final_Revision.xlsx"
                ),
                2020: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final_Revision.xlsx"
                ),
                2021: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final_Revision.xlsx"
                ),
                2022: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Environmental_Information_{year}_Final.xlsx"
                ),
                2023: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Envir_Infor_{year}_Final.xlsx"
                ),
                2024: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Envir_Infor_{year}_Final.xlsx"
                ),
                2025: downloads_folder(
                    f"eia923/f923_{year}/EIA923_Schedule_8_Annual_Envir_Infor_{year}_Final.xlsx"
                ),
            }[year]
            header_row = 4
        elif (year == current_early_release_year) and (
            current_early_release_year != latest_validated_year
        ):
            schedule_8_filename = find_downloaded_file(
                f"eia923/f923_{year}er",
                f"EIA923_Schedule_8*{year}*Early*Release*.xlsx",
            )
            header_row = 5

        emissions_controls_eia923 = pd.read_excel(
            io=schedule_8_filename,
            sheet_name="8C Air Emissions Control Info",
            header=header_row,
            names=emissions_controls_eia923_names,
            dtype=get_dtypes(),
            na_values=".",
            parse_dates=["report_date"],
            date_format={"particulate_test_date": "%m-%Y", "so2_test_date": "%m-%Y"},
        )
        emissions_controls_eia923["report_date"] = emissions_controls_eia923[
            "report_date"
        ].astype("datetime64[s]")
        emissions_controls_eia923["equipment_tech_description"] = (
            emissions_controls_eia923["equipment_tech_description"].str.lower()
        )
    else:
        logger.warning(
            "Emissions control data prior to 2012 has not been integrated into the data pipeline."
        )
        logger.warning(
            "This may overestimate SO2 and NOx emissions calculated from EIA-923 data."
        )
        emissions_controls_eia923 = pd.DataFrame(
            columns=emissions_controls_eia923_names
        )

    return emissions_controls_eia923


def load_unit_to_boiler_associations(year: int) -> pd.DataFrame:
    """Creates a table that associates EPA units with EIA boilers

    Args:
        year (int): a four-digit year

    Returns:
        pd.DataFrame: table relating EPA emission unit IDs to EIA boiler IDs.
    """
    subplant_crosswalk = pd.read_csv(
        outputs_folder(f"{year}/subplant_crosswalk_{year}.csv.zip"),
        dtype=get_dtypes(),
    )
    boiler_generator_assn = load_pudl_table("core_eia860__assn_boiler_generator", year)
    unit_boiler_assn = subplant_crosswalk.merge(
        boiler_generator_assn, how="left", on=["plant_id_eia", "generator_id"]
    )
    unit_boiler_assn = unit_boiler_assn[
        ["plant_id_eia", "emissions_unit_id_epa", "boiler_id"]
    ]

    return unit_boiler_assn


def load_default_gtn_ratios() -> pd.DataFrame:
    """Loads the default gross to net generation ratios.

    Returns:
        pd.DataFrame: gross-to-net ratios table for prime movers.
    """
    default_gtn = pd.read_csv(
        reference_table_folder("default_gross_to_net_ratios.csv"),
        dtype=get_dtypes(),
    )[["prime_mover_code", "default_gtn_ratio"]]

    return default_gtn


def load_egrid_plant_file(year: int) -> pd.DataFrame:
    """Loads plant level data from eGRID.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: the eGRID annual plant level data frame.
    """
    egrid_columns = [
        "BACODE",
        "PSTATABB",
        "PLPRMFL",
        "ORISPL",
        "PNAME",
        "PLGENATN",
        "PLGENATR",
        "PLHTIANT",
        "UNNOX",
        "UNSO2",
        "UNCO2",
        "UNCH4",
        "UNN2O",
        "UNHTIT",
        "UNHTIOZT",
        "UNHTISRC",
        "UNHOZSRC",
        "PLCO2AN",
        "PLCO2EQA",
        "PLNOXAN",
        "PLSO2AN",
        "CHPFLAG",
        "ELCALLOC",
    ]

    # load plant level data from egrid
    egrid_plant = pd.read_excel(
        downloads_folder(f"egrid/egrid{year}_data.xlsx"),
        sheet_name=f"PLNT{str(year)[-2:]}",
        header=1,
        usecols=egrid_columns,
    )

    # calculate total net generation from reported renewable and nonrenewable generation
    egrid_plant["net_generation_mwh"] = (
        egrid_plant["PLGENATN"] + egrid_plant["PLGENATR"]
    )
    egrid_plant = egrid_plant.drop(columns=["PLGENATN", "PLGENATR"])

    # rename the columns
    egrid_plant = egrid_plant.rename(
        columns={
            "BACODE": "ba_code",
            "PSTATABB": "state",
            "PLPRMFL": "plant_primary_fuel",
            "ORISPL": "plant_id_egrid",
            "PNAME": "plant_name_eia",
            "UNHTIT": "fuel_consumed_mmbtu",
            "PLHTIANT": "fuel_consumed_for_electricity_mmbtu",
            "UNCO2": "co2_mass_lb",  # this is actually in tons, but we are converting in the next step
            "UNCH4": "ch4_mass_lb",
            "UNN2O": "n2o_mass_lb",
            "UNNOX": "nox_mass_lb",  # this is actually in tons, but we are converting in the next step
            "UNSO2": "so2_mass_lb",  # this is actually in tons, but we are converting in the next step
            "PLCO2AN": "co2_mass_lb_for_electricity_adjusted",  # this is actually in tons, but we are converting in the next step
            "PLCO2EQA": "co2e_mass_lb_for_electricity_adjusted",  # this is actually in tons, but we are converting in the next step
            "PLNOXAN": "nox_mass_lb_for_electricity_adjusted",  # this is actually in tons, but we are converting in the next step
            "PLSO2AN": "so2_mass_lb_for_electricity_adjusted",  # this is actually in tons, but we are converting in the next step
            "CHPFLAG": "chp_flag",
            "ELCALLOC": "chp_electric_allocation_factor",
            "UNHTIOZT": "fuel_consumed_mmbtu_ozone_season",
            "UNHTISRC": "fuel_data_source_annual",
            "UNHOZSRC": "fuel_data_source_ozone",
        }
    )

    # convert mass tons to lb
    egrid_plant["co2_mass_lb"] = (
        egrid_plant["co2_mass_lb"] * ConversionFactors.short_ton_to_lbs
    )
    egrid_plant["nox_mass_lb"] = (
        egrid_plant["nox_mass_lb"] * ConversionFactors.short_ton_to_lbs
    )
    egrid_plant["so2_mass_lb"] = (
        egrid_plant["so2_mass_lb"] * ConversionFactors.short_ton_to_lbs
    )
    egrid_plant["co2_mass_lb_for_electricity_adjusted"] = (
        egrid_plant["co2_mass_lb_for_electricity_adjusted"]
        * ConversionFactors.short_ton_to_lbs
    )
    egrid_plant["co2e_mass_lb_for_electricity_adjusted"] = (
        egrid_plant["co2e_mass_lb_for_electricity_adjusted"]
        * ConversionFactors.short_ton_to_lbs
    )
    egrid_plant["nox_mass_lb_for_electricity_adjusted"] = (
        egrid_plant["nox_mass_lb_for_electricity_adjusted"]
        * ConversionFactors.short_ton_to_lbs
    )
    egrid_plant["so2_mass_lb_for_electricity_adjusted"] = (
        egrid_plant["so2_mass_lb_for_electricity_adjusted"]
        * ConversionFactors.short_ton_to_lbs
    )

    # if egrid has a missing value for co2 for a clean plant, replace with zero
    egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS),
        "co2_mass_lb_for_electricity_adjusted",
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS),
        "co2_mass_lb_for_electricity_adjusted",
    ].fillna(0)
    egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS), "co2_mass_lb"
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS), "co2_mass_lb"
    ].fillna(0)

    # reorder the columns
    egrid_plant = egrid_plant[
        [
            "ba_code",
            "state",
            "plant_id_egrid",
            "plant_name_eia",
            "plant_primary_fuel",
            "chp_flag",
            "chp_electric_allocation_factor",
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
            "co2_mass_lb",
            "co2_mass_lb_for_electricity_adjusted",
            "ch4_mass_lb",
            "n2o_mass_lb",
            "co2e_mass_lb_for_electricity_adjusted",
            "nox_mass_lb",
            "nox_mass_lb_for_electricity_adjusted",
            "so2_mass_lb",
            "so2_mass_lb_for_electricity_adjusted",
            "fuel_consumed_mmbtu_ozone_season",
            "fuel_data_source_annual",
            "fuel_data_source_ozone",
        ]
    ]

    # create a column for eia id
    egrid_plant = add_egrid_plant_id(egrid_plant, from_id="egrid", to_id="eia")

    return egrid_plant


def load_egrid_ba_file(year: int) -> pd.DataFrame:
    """Loads in eGRID balancing authority totals.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table of total net generation, fuel consumed and CO2 emission by
            balancing authority.
    """
    # load egrid BA totals
    egrid_ba = pd.read_excel(
        downloads_folder(f"egrid/egrid{year}_data.xlsx"),
        sheet_name=f"BA{str(year)[-2:]}",
        header=1,
        usecols=["BANAME", "BACODE", "BAHTIANT", "BANGENAN", "BACO2AN"],
    )
    # rename the columns
    egrid_ba = egrid_ba.rename(
        columns={
            "BANAME": "ba_name",
            "BACODE": "ba_code",
            "BAHTIANT": "fuel_consumed_for_electricity_mmbtu",
            "BANGENAN": "net_generation_mwh",
            "BACO2AN": "co2_mass_lb_adjusted",
        }
    )
    egrid_ba = egrid_ba.sort_values(by="ba_code", ascending=True)
    egrid_ba["co2_mass_lb_adjusted"] = (
        egrid_ba["co2_mass_lb_adjusted"] * ConversionFactors.short_ton_to_lbs
    )

    return egrid_ba


def add_egrid_plant_id(df: pd.DataFrame, from_id: str, to_id: str) -> pd.DataFrame:
    """Add a plant_id_{to_id} column in the input plant data frame.

    Args:
        df (pd.DataFrame): input plant data.
        from_id (str): ID to convert from.
        to_id (str): ID to add.

    Returns:
        pd.DataFrame: original data frame with additional target ID column
    """
    egrid_crosswalk = pd.read_csv(
        reference_table_folder("eGRID_crosswalk_of_EIA_ID_to_EPA_ID.csv"),
        dtype=get_dtypes(),
    )
    id_map = dict(
        zip(
            list(egrid_crosswalk[f"plant_id_{from_id}"]),
            list(egrid_crosswalk[f"plant_id_{to_id}"]),
        )
    )

    df[f"plant_id_{to_id}"] = df[f"plant_id_{from_id}"]
    df[f"plant_id_{to_id}"] = df[f"plant_id_{to_id}"].replace(id_map)

    return df


def test():
    return Path.cwd()
