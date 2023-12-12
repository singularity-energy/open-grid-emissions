import pandas as pd
import numpy as np
import sqlalchemy as sa
import warnings
from pathlib import Path

from oge.column_checks import get_dtypes
from oge.filepaths import downloads_folder, manual_folder, outputs_folder
from oge.validation import validate_unique_datetimes
from oge.logging_util import get_logger

from pudl.metadata.fields import apply_pudl_dtypes

logger = get_logger(__name__)

# initialize the pudl_engine
PUDL_ENGINE = sa.create_engine("sqlite:///" + downloads_folder("pudl/pudl.sqlite"))


def load_cems_data(year):
    """
    Loads CEMS data for the specified year from the PUDL database
    Inputs:
        year: the year for which data should be retrieved (YYYY)
    Returns:
        cems: pandas dataframe with hourly CEMS data
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
        downloads_folder("pudl/hourly_emissions_epacems.parquet"),
        filters=[["year", "==", year]],
        columns=cems_columns,
    )
    # convert to tz-naive datetime to allow for dtype application
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(None)
    cems = apply_pudl_dtypes(cems)
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(
        "UTC"
    )

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
    cems["co2_mass_lb"] = cems["co2_mass_tons"] * 2000

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

    validate_unique_datetimes(cems, "cems", ["plant_id_eia", "emissions_unit_id_epa"])

    return cems


def load_cems_ids(start_year, end_year):
    """Loads CEMS ids for multiple years."""

    # load cems data
    cems = pd.read_parquet(
        downloads_folder("pudl/hourly_emissions_epacems.parquet"),
        filters=[["year", ">=", start_year], ["year", "<=", end_year]],
        columns=["plant_id_eia", "emissions_unit_id_epa"],
    ).drop_duplicates()
    cems = apply_pudl_dtypes(cems)

    return cems


def load_cems_gross_generation(start_year, end_year):
    """Loads hourly CEMS gross generation data for multiple years."""

    # specify the columns to use from the CEMS database
    cems_columns = [
        "plant_id_eia",
        "emissions_unit_id_epa",
        "operating_datetime_utc",
        "operating_time_hours",
        "gross_load_mw",
    ]

    # load cems data
    cems = pd.read_parquet(
        downloads_folder("pudl/hourly_emissions_epacems.parquet"),
        filters=[["year", ">=", start_year], ["year", "<=", end_year]],
        columns=cems_columns,
    )
    # convert to tz-naive datetime to allow for dtype application
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(None)
    cems = apply_pudl_dtypes(cems)
    cems["operating_datetime_utc"] = cems["operating_datetime_utc"].dt.tz_localize(
        "UTC"
    )

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


def add_report_date(df):
    """
    Add a report date column to the cems data based on the plant's local timezone

    Args:
        df (pd.Dataframe): dataframe containing 'plant_id_eia' and 'datetime_utc' columns
    Returns:
        Original dataframe with 'report_date' column added
    """
    plant_timezone = load_pudl_table(
        "plants_entity_eia", columns=["plant_id_eia", "timezone"]
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
    df["report_date"] = np.NaN

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


def load_ghg_emission_factors():
    """
    Read in the table of emissions factors and convert to lb/mmbtu
    """

    efs = pd.read_csv(
        manual_folder("emission_factors_for_co2_ch4_n2o.csv"),
        dtype=get_dtypes(),
    )

    # convert co2 mass in short tons to lb
    efs["co2_tons_per_mmbtu"] = efs["co2_tons_per_mmbtu"] * 2000

    # rename the columns
    efs = efs.rename(columns={"co2_tons_per_mmbtu": "co2_lb_per_mmbtu"})

    return efs


def load_nox_emission_factors():
    """Read in the NOx emission factors from eGRID Table C2."""
    emission_factors = pd.read_csv(
        manual_folder("emission_factors_for_nox.csv"),
        dtype=get_dtypes(),
    )

    # standardize units as lower case
    emission_factors["emission_factor_denominator"] = emission_factors[
        "emission_factor_denominator"
    ].str.lower()

    return emission_factors


def load_so2_emission_factors():
    """
    Read in the SO2 emission factors from eGRID Table C3.

    The SO2 emission rate depends on the sulfur content of fuel, so it is
    reported in Table C3 as a formula like `123*S`.
    """
    df = pd.read_csv(
        manual_folder("emission_factors_for_so2.csv"),
        dtype=get_dtypes(),
    )

    # Add a boolean column that reports whether the emission factor is a formula or value.
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
):
    """
    Loads a table from the pudl database.

    `table_name` must be one of the options specified in the data dictionary:
    https://catalystcoop-pudl.readthedocs.io/en/latest/data_dictionaries/pudl_db.html

    There are multiple options for date filtering:
        - if `year` is not specified, all years will be loaded
        - if `year` is specified, but not `end_year`, only a single year will be loaded
        - if both `year` and `end_year` are specified, all years in that range inclusive
          will be loaded. `end_year` must be >= `year`.

    If a list of `columns` is passed, only those columns will be returned. Otherwise,
    all columns will be returned.
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

    table = apply_pudl_dtypes(table)

    return table


def load_epa_eia_crosswalk_from_raw(year):
    """
    Read in the manual EPA-EIA Crosswalk table downloaded from the EPA website.

    This is only used in OGE to access the CAMD_FUEL_TYPE column, which is dropped from
    the PUDL version of the table.
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
    crosswalk["CAMD_UNIT_ID"] = crosswalk["CAMD_UNIT_ID"].str.lstrip("0")

    # some eia plant ids are missing. Let us assume that the EIA and EPA plant ids match in this case
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

    # fill missing values in the energy_source_code_eia column with values from the energy_source_code_epa column
    crosswalk["energy_source_code_eia"] = crosswalk["energy_source_code_eia"].fillna(
        crosswalk["energy_source_code_epa"]
    )

    # load manually inputted data
    crosswalk_manual = pd.read_csv(
        manual_folder("epa_eia_crosswalk_manual.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["notes"])

    # load EIA-860 data
    gen_esc_860 = load_pudl_table(
        "generators_eia860",
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


def load_epa_eia_crosswalk(year):
    """
    Read in the manual EPA-EIA Crosswalk table.
    """

    crosswalk = load_pudl_table("epacamd_eia")

    # load manually inputted data
    crosswalk_manual = pd.read_csv(
        manual_folder("epa_eia_crosswalk_manual.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["notes"])

    # concat this data with the main table
    crosswalk = pd.concat(
        [crosswalk, crosswalk_manual],
        axis=0,
    )

    # load EIA-860 data
    gen_ids = load_pudl_table(
        "generators_eia860", year, columns=["plant_id_eia", "generator_id"]
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
    level, conversion_type, threshold_column, lower_threshold, upper_threshold, year
):
    """
    Loads gross-to-net generation conversion factors calculated in `gross_to_net_generation`.

    If you wanted to load subplant regression results with an adjusted r2 greater than 0.9, you would specify:
        level = 'subplant'
        conversion_type = 'regression'
        threshold_column = 'rsqaured_adjusted'
        lower_threshold = 0.9
        upper_threshold = None

    Args:
        level: Aggregation level, either 'plant' or 'subplant'
        conversion_type: which data to load, either 'regression' or 'ratio'
        threshold_column: name of the column to be used to filter values
        lower_threshold: the value below which in threshold_column the conversion factors should be removed
        upper_threshold: the value above which in threshold_column the conversion factors should be removed
    Returns:
        gtn_data: pandas dataframe containing revevant keys and conversion factors
    """
    gtn_data = pd.read_csv(
        outputs_folder(f"gross_to_net/{level}_gross_to_net_{conversion_type}.csv"),
        dtype=get_dtypes(),
    )

    # filter the data based on the upper and lower thresholds, if specified
    if lower_threshold is not None:
        gtn_data = gtn_data[gtn_data[threshold_column] >= lower_threshold]

    if upper_threshold is not None:
        gtn_data = gtn_data[gtn_data[threshold_column] <= upper_threshold]

    # if loading regression data, add a count of units in each subplant to the regression results
    if conversion_type == "regression":
        subplant_crosswalk = pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
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

        # divide the intercept by the number of units in each subplant to evenly distribute this to each unit
        gtn_data["intercept"] = gtn_data["intercept"] / gtn_data[f"units_in_{level}"]

    # make sure the report date column is a datetime if loading ratios
    if conversion_type == "ratio":
        gtn_data["report_date"] = pd.to_datetime(gtn_data["report_date"]).astype(
            "datetime64[s]"
        )

    return gtn_data


def load_ipcc_gwp():
    """Load a table containing global warming potential (GWP) values for CO2, CH4, and N2O."""
    return pd.read_csv(manual_folder("ipcc_gwp.csv"), dtype=get_dtypes())


def load_raw_eia930_data(year, description):
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


def load_ba_reference():
    return pd.read_csv(
        manual_folder("ba_reference.csv"),
        dtype=get_dtypes(),
        parse_dates=["activation_date", "retirement_date"],
    )


def load_diba_data(year):
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


def ba_timezone(ba, type):
    """
    Retrieves the timezone for a single balancing area.
    Args:
        ba: string containing the ba_code
        type: either 'reporting_eia930' or 'local'. Reporting will return the TZ used by the BA when reporting to EIA-930, local will return the actual local tz
    """

    tz = pd.read_csv(
        manual_folder("ba_reference.csv"),
        usecols=["ba_code", f"timezone_{type}"],
    )
    tz = tz.loc[tz["ba_code"] == ba, f"timezone_{type}"]

    if len(tz) == 0:
        raise UserWarning(
            f"The BA {ba} does not have a timezone specified in data/manual/ba_reference.csv. Please add."
        )
    else:
        tz = tz.item()

    return tz


def load_emissions_controls_eia923(year: int):
    emissions_controls_eia923_names = [
        "report_date",
        "plant_id_eia",
        "equipment_tech_description",
        "particulate_control_id_eia",
        "so2_control_id_eia",
        "nox_control_id_eia",
        "mercury_control_id_eia",
        "operational_status",
        "hours_in_service",
        "annual_nox_emission_rate_lb_per_mmbtu",
        "ozone_season_nox_emission_rate_lb_per_mmbtu",
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

    datetime_columns = ["report_date", "particulate_test_date", "so2_test_date"]

    # For 2012-2015 and earlier, mercury emission rate is not reported in EIA923, so we need
    # to remove that column to avoid an error.
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
        }[year]

        emissions_controls_eia923 = pd.read_excel(
            io=schedule_8_filename,
            sheet_name="8C Air Emissions Control Info",
            header=4,
            names=emissions_controls_eia923_names,
            dtype=get_dtypes(),
            na_values=".",
            parse_dates=datetime_columns,
            date_format={"particulate_test_date": "%m-%Y", "so2_test_date": "%m-%Y"},
        )
        emissions_controls_eia923["report_date"] = emissions_controls_eia923[
            "report_date"
        ].astype("datetime64[s]")
    else:
        logger.warning(
            "Emissions control data prior to 2014 has not been integrated into the data pipeline."
        )
        logger.warning(
            "This may overestimate SO2 and NOx emissions calculated from EIA-923 data."
        )
        emissions_controls_eia923 = pd.DataFrame(
            columns=emissions_controls_eia923_names
        )

    return emissions_controls_eia923


def load_unit_to_boiler_associations(year):
    """Creates a table that associates EPA units with EIA boilers"""
    subplant_crosswalk = pd.read_csv(
        outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
        dtype=get_dtypes(),
    )
    boiler_generator_assn = load_pudl_table("boiler_generator_assn_eia860", year)
    unit_boiler_assn = subplant_crosswalk.merge(
        boiler_generator_assn, how="left", on=["plant_id_eia", "generator_id"]
    )
    unit_boiler_assn = unit_boiler_assn[
        ["plant_id_eia", "emissions_unit_id_epa", "boiler_id"]
    ]

    return unit_boiler_assn


def load_default_gtn_ratios():
    """Read in the default gross to net generation ratios."""
    default_gtn = pd.read_csv(
        manual_folder("default_gross_to_net_ratios.csv"),
        dtype=get_dtypes(),
    )[["prime_mover_code", "default_gtn_ratio"]]

    return default_gtn


def test():
    return Path.cwd()
