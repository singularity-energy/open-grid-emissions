import pandas as pd
import numpy as np
import sqlalchemy as sa
import warnings
from pathlib import Path

from oge.column_checks import get_dtypes, apply_dtypes
from oge.filepaths import downloads_folder, reference_table_folder, outputs_folder
import oge.validation as validation
from oge.logging_util import get_logger
from oge.constants import (
    CLEAN_FUELS,
    ConversionFactors,
    earliest_data_year,
    earliest_validated_year,
    latest_validated_year,
)

logger = get_logger(__name__)

# initialize the pudl_engine
PUDL_ENGINE = sa.create_engine("sqlite:///" + downloads_folder("pudl/pudl.sqlite"))


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
        "steam_load_1000_lbs",
        "co2_mass_tons",
        "co2_mass_measurement_code",
        "nox_mass_lbs",
        "nox_mass_measurement_code",
        "so2_mass_lbs",
        "so2_mass_measurement_code",
        "heat_content_mmbtu",
    ]

    # load the CEMS data
    cems = pd.read_parquet(
        downloads_folder("pudl/core_epacems__hourly_emissions.parquet"),
        filters=[["year", "==", year]],
        columns=cems_columns,
    )
    # convert to tz-naive datetime to allow for dtype application
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(None)
    cems = apply_dtypes(cems)
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(
        "UTC"
    )

    # update the plant_id_eia column using manual matches
    cems = update_epa_to_eia_map(cems)

    cems = cems.rename(
        columns={
            "operating_datetime_utc": "datetime_utc",
            "heat_content_mmbtu": "fuel_consumed_mmbtu",
            "steam_load_1000_lbs": "steam_load_1000_lb",
            "nox_mass_lbs": "nox_mass_lb",
            "so2_mass_lbs": "so2_mass_lb",
            "gross_load_mw": "gross_generation_mwh",  # we will convert this to mwh
        }
    )

    # fill any missing values for steam load with zero
    cems["steam_load_1000_lb"] = cems["steam_load_1000_lb"].fillna(0)

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
            "steam_load_1000_lb",
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

    validation.validate_unique_datetimes(
        year, cems, "cems", ["plant_id_eia", "emissions_unit_id_epa"]
    )

    return cems


def load_cems_ids() -> pd.DataFrame:
    """Loads a dataframe of all unique plant_id_eia, emissions_unit_id_epa combinations
    that exist from `constants.earliest_data_year` to `constants.latest_validated_year`.
    This is used in the process of creating subplant_ids to ensure complete coverage.

    Returns:
        pd.DataFrame: a two column table relating plant_id_eia to emissions_unit_id_epa.
    """
    # although we could directly load all years at once from the cems parquet file,
    # this would lead to a memoryerror, so we load one year at a time and drop
    # duplicates before concatenating the next year to the dataframe
    cems_ids = []
    # The `constants.earliest_data_year` is 2005
    for year in range(earliest_data_year, latest_validated_year + 1):
        cems_id_year = pd.read_parquet(
            downloads_folder("pudl/core_epacems__hourly_emissions.parquet"),
            filters=[["year", "==", year]],
            columns=["plant_id_epa", "plant_id_eia", "emissions_unit_id_epa"],
        ).drop_duplicates()
        cems_ids.append(cems_id_year)
        cems_ids = [pd.concat(cems_ids, axis=0).drop_duplicates()]

    cems_ids = (
        pd.concat(cems_ids, axis=0)
        .drop_duplicates()
        .sort_values(by=["plant_id_eia", "emissions_unit_id_epa"])
    )

    cems_ids = apply_dtypes(cems_ids)

    # update the plant_id_eia column using manual matches
    cems_ids = update_epa_to_eia_map(cems_ids)

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
    # only keep data for years <= the year
    # this avoids using potentially preliminary early-release data
    complete_gens = complete_gens[
        (complete_gens["report_date"].dt.year >= earliest_data_year)
        & (complete_gens["report_date"].dt.year <= latest_validated_year)
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
    cancelled_or_proposed_status_codes = ["CN", "IP", "P", "L", "T"]
    complete_gens = complete_gens[
        ~complete_gens["operational_status_code"].isin(
            cancelled_or_proposed_status_codes
        )
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
            (complete_gens["report_date"].dt.year < earliest_validated_year)
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
            & (complete_gens["report_date"].dt.year < latest_validated_year)
            & (complete_gens["operational_status_code"] != "TS")
        )
    ]

    ####################
    # merge into complete_gens and fill missing operating dates with the EIA-860 data
    generator_data_from_eia860 = load_raw_eia860_generator_dates_and_unit_ids(
        latest_validated_year
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
    plant_geographical_eia860 = pd.read_excel(
        downloads_folder(f"eia860/eia860{year}/2___Plant_Y{year}.xlsx"),
        header=1,
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
    generator_op_dates_eia860 = pd.read_excel(
        downloads_folder(f"eia860/eia860{year}/3_1_Generator_Y{year}.xlsx"),
        header=1,
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
            "Operating Month": "Month",
            "Operating Year": "Year",
            "Unit Code": "unit_id_eia",
        }
    )

    # create a datetime column from the month and year
    generator_op_dates_eia860["operating_date_eia"] = pd.to_datetime(
        generator_op_dates_eia860[["Year", "Month"]].assign(Day=1)
    )
    generator_op_dates_eia860 = generator_op_dates_eia860.drop(
        columns=["Month", "Year"]
    )

    # load unit codes for proposed generators
    proposed_unit_ids_eia860 = (
        pd.read_excel(
            downloads_folder(f"eia860/eia860{year}/3_1_Generator_Y{year}.xlsx"),
            sheet_name="Proposed",
            header=1,
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

    # concat the data together
    generator_data_from_eia860 = pd.concat(
        [generator_op_dates_eia860, proposed_unit_ids_eia860], axis=0
    )

    # create a numeric version of the ID, starting at 1
    generator_data_from_eia860["unit_id_eia_numeric"] = (
        generator_data_from_eia860.groupby(
            ["plant_id_eia"]
        )["unit_id_eia"].transform(lambda x: pd.factorize(x)[0] + 1)
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


def load_cems_gross_generation(start_year: int, end_year: int) -> pd.DataFrame:
    """Loads hourly CEMS gross generation data for multiple years.

    Args:
        start_year (int): a four-digit year.
        end_year (int): a four-digit year

    Returns:
        pd.DataFrame: CEMS gross generation data over the required year aggregated by
            plant, unit and month.
    """

    # specify the columns to use from the CEMS database
    cems_columns = [
        "plant_id_epa",
        "plant_id_eia",
        "emissions_unit_id_epa",
        "operating_datetime_utc",
        "operating_time_hours",
        "gross_load_mw",
    ]

    # load cems data
    cems = pd.read_parquet(
        downloads_folder("pudl/core_epacems__hourly_emissions.parquet"),
        filters=[["year", ">=", start_year], ["year", "<=", end_year]],
        columns=cems_columns,
    )
    # convert to tz-naive datetime to allow for dtype application
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(None)
    cems = apply_dtypes(cems)
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(
        "UTC"
    )

    # update the plant_id_eia column using manual matches
    cems = update_epa_to_eia_map(cems)

    # only keep values when the plant was operating
    # this will help speed up calculations and allow us to add this data back later
    cems = cems[(cems["gross_load_mw"] > 0) | (cems["operating_time_hours"] > 0)]

    # rename the heat content column to use the convention used in the EIA data
    cems = cems.rename(
        columns={
            "operating_datetime_utc": "datetime_utc",
            "gross_load_mw": "gross_generation_mwh",
        }
    )

    # add a report date
    cems = add_report_date(cems)

    cems = cems[
        [
            "plant_id_eia",
            "emissions_unit_id_epa",
            "report_date",
            "gross_generation_mwh",
        ]
    ]

    # group data by plant, unit, month
    cems = cems.groupby(
        ["plant_id_eia", "emissions_unit_id_epa", "report_date"], dropna=False
    ).sum()

    return cems


def update_epa_to_eia_map(cems_df: pd.DataFrame) -> pd.DataFrame:
    """Updates the plant_id_eia column in the CEMS data frame loaded from pudl based
    on the manual epa_eia_crosswalk_manual table

    Args:
        cems_df (pd.DataFrame): input CEMS data frame.

    Returns:
        pd.DataFrame: CEMS data frame with updated plant_id_eia column.
    """
    # load the manual table
    manual_plant_map = pd.read_csv(
        reference_table_folder("epa_eia_crosswalk_manual.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["notes"])

    # only keep rows where the epa and eia plant ids don't match
    manual_plant_map = manual_plant_map.loc[
        manual_plant_map["plant_id_epa"] != manual_plant_map["plant_id_eia"],
        ["plant_id_epa", "emissions_unit_id_epa", "plant_id_eia"],
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
    """Read in the NOx emission factors from eGRID Table C2.

    Returns:
        pd.DataFrame: NOx emission factors table.
    """
    emission_factors = pd.read_csv(
        reference_table_folder("emission_factors_for_nox.csv"),
        dtype=get_dtypes(),
    )

    # standardize units as lower case
    emission_factors["emission_factor_denominator"] = emission_factors[
        "emission_factor_denominator"
    ].str.lower()

    return emission_factors


def load_so2_emission_factors() -> pd.DataFrame:
    """Read in the SO2 emission factors from eGRID Table C3.

    The SO2 emission rate depends on the sulfur content of fuel, so it is
    reported in Table C3 as a formula like `123*S`.

    Returns:
        pd.DataFrame: SO2 emission factors table.
    """
    df = pd.read_csv(
        reference_table_folder("emission_factors_for_so2.csv"),
        dtype=get_dtypes(),
    )

    # Add a boolean column that reports whether the emission factor is a formula or
    # value.
    df["multiply_by_sulfur_content"] = (
        df["emission_factor"].str.contains("*", regex=False).astype(int)
    )

    # Extract the numeric coefficient from the emission factor.
    df["emission_factor"] = (
        df["emission_factor"].str.replace("*S", "", regex=False).astype(float)
    )

    # standardize units as lower case
    df["emission_factor_denominator"] = df["emission_factor_denominator"].str.lower()

    return df


def load_pudl_table(
    table_name: str, year: int = None, columns: list[str] = None, end_year: int = None
) -> pd.DataFrame:
    """Loads a table from the pudl database.

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

    Returns:
        pd.DataFrame: the required table.
    """
    if columns is None:
        columns_to_select = "*"
    else:
        columns_to_select = ", ".join(columns)

    if year is None:
        # load the table without filtering dates
        table = pd.read_sql(
            f"SELECT {columns_to_select} FROM {table_name}",
            PUDL_ENGINE,
        )
    elif year is not None and end_year is None:
        # load the table for a single year
        table = pd.read_sql(
            f"SELECT {columns_to_select} FROM {table_name} WHERE \
                report_date >= '{year}-01-01' AND report_date < '{year + 1}-01-01'",
            PUDL_ENGINE,
        )
    else:
        # load the for the specified years
        table = pd.read_sql(
            f"SELECT {columns_to_select} FROM {table_name} WHERE \
                report_date >= '{year}-01-01' AND report_date < '{end_year + 1}-01-01'",
            PUDL_ENGINE,
        )

    table = apply_dtypes(table)

    if table.empty:
        logger.warning(f"{table_name} is empty")

    return table


def load_epa_eia_crosswalk_from_raw(year: int) -> pd.DataFrame:
    """Read in the EPA-EIA Crosswalk table downloaded from the EPA website.

    This is only used in OGE to access the CAMD_FUEL_TYPE column, which is dropped from
    the PUDL version of the table.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: EPA-EIA crosswalk table.
    """

    crosswalk = pd.read_csv(
        downloads_folder("epa/epa_eia_crosswalk.csv"),
        usecols=[
            "CAMD_PLANT_ID",
            "EIA_PLANT_ID",
            "CAMD_UNIT_ID",
            "EIA_GENERATOR_ID",
            "EIA_BOILER_ID",
            "EIA_FUEL_TYPE",
            "CAMD_FUEL_TYPE",
        ],
    )

    # remove leading zeros from the generator id and emissions_unit_id_epa
    crosswalk["EIA_GENERATOR_ID"] = crosswalk["EIA_GENERATOR_ID"].str.lstrip("0")

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
            "EIA_GENERATOR_ID": "generator_id",
            "EIA_BOILER_ID": "boiler_id",
            "EIA_FUEL_TYPE": "energy_source_code_eia",
            "CAMD_FUEL_TYPE": "energy_source_code_epa",
        }
    )

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
    crosswalk["energy_source_code_epa"] = crosswalk["energy_source_code_epa"].replace(
        camd_to_eia_fuel_type
    )

    # fill missing values in the energy_source_code_eia column with values from the
    # energy_source_code_epa column
    crosswalk["energy_source_code_eia"] = crosswalk["energy_source_code_eia"].fillna(
        crosswalk["energy_source_code_epa"]
    )

    # load manually inputted data
    crosswalk_manual = pd.read_csv(
        reference_table_folder("epa_eia_crosswalk_manual.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["notes"])

    # load EIA-860 data
    gen_esc_860 = load_pudl_table(
        "core_eia860__scd_generators",
        year,
        columns=["plant_id_eia", "generator_id", "energy_source_code_1"],
    )

    # merge the energy source code from EIA-860
    crosswalk_manual = crosswalk_manual.merge(
        gen_esc_860, how="left", on=["plant_id_eia", "generator_id"], validate="m:1"
    ).rename(columns={"energy_source_code_1": "energy_source_code_eia"})

    # concat this data with the main table
    crosswalk = pd.concat(
        [crosswalk, crosswalk_manual],
        axis=0,
    )

    # merge in any plants that are missing from the EPA crosswalk but appear in EIA-860
    crosswalk = crosswalk.merge(
        gen_esc_860, how="outer", on=["plant_id_eia", "generator_id"], validate="m:1"
    )
    crosswalk["plant_id_epa"] = crosswalk["plant_id_epa"].fillna(
        crosswalk["plant_id_eia"]
    )
    crosswalk["energy_source_code_eia"] = crosswalk["energy_source_code_eia"].fillna(
        crosswalk["energy_source_code_1"]
    )
    crosswalk["energy_source_code_epa"] = crosswalk["energy_source_code_epa"].fillna(
        crosswalk["energy_source_code_1"]
    )
    crosswalk = crosswalk.drop(columns=["energy_source_code_1"])

    return crosswalk


def load_epa_eia_crosswalk(year: int) -> pd.DataFrame:
    """Read in the manual EPA-EIA Crosswalk table.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: manual EIA-EPA crosswalk table.
    """

    crosswalk = load_pudl_table("core_epa__assn_eia_epacamd")

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
    dibas = load_raw_eia930_data(year, "INTERCHANGE")
    dibas = dibas[
        [
            "Balancing Authority",
            "Directly Interconnected Balancing Authority",
            "Region",
            "DIBA Region",
        ]
    ].drop_duplicates()
    dibas = dibas.rename(
        columns={
            "Balancing Authority": "ba_code",
            "Directly Interconnected Balancing Authority": "diba_code",
            "Region": "ba_region",
            "DIBA Region": "diba_region",
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
        }[year]

        emissions_controls_eia923 = pd.read_excel(
            io=schedule_8_filename,
            sheet_name="8C Air Emissions Control Info",
            header=4,
            names=emissions_controls_eia923_names,
            dtype=get_dtypes(),
            na_values=".",
            parse_dates=["report_date"],
            date_format={"particulate_test_date": "%m-%Y", "so2_test_date": "%m-%Y"},
        )
        emissions_controls_eia923["report_date"] = emissions_controls_eia923[
            "report_date"
        ].astype("datetime64[s]")
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

    # We also want to remove any plants that are located in Puerto Rico
    egrid_plant = egrid_plant[(egrid_plant["state"] != "PR")]

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
