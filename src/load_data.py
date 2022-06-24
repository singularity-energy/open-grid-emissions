import pandas as pd
import sqlalchemy as sa

import pudl.output.pudltabl

from src.column_checks import get_dtypes


def load_cems_data(year):
    """
    Loads CEMS data for the specified year from the PUDL database
    Inputs:
        year: the year for which data should be retrieved (YYYY)
    Returns:
        cems: pandas dataframe with hourly CEMS data
    """
    # specify the path to the CEMS data
    cems_path = f"../data/downloads/pudl/pudl_data/parquet/epacems/year={year}"

    # specify the columns to use from the CEMS database
    cems_columns = [
        "plant_id_eia",
        "unitid",
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
    cems = pd.read_parquet(cems_path, columns=cems_columns)

    # rename cems plant_id_eia to plant_id_epa (PUDL simply renames the ORISPL_CODE column from the raw CEMS data as 'plant_id_eia' without actually crosswalking to the EIA id)
    # rename the heat content column to use the convention used in the EIA data
    cems = cems.rename(
        columns={
            "operating_datetime_utc": "datetime_utc",
            "plant_id_eia": "plant_id_epa",
            "heat_content_mmbtu": "fuel_consumed_mmbtu",
            "steam_load_1000_lbs": "steam_load_1000_lb",
            "nox_mass_lbs": "nox_mass_lb",
            "so2_mass_lbs": "so2_mass_lb",
            "gross_load_mw": "gross_generation_mwh",  # we are going to convert this in a later step
        }
    )

    # if the unitid has any leading zeros, remove them
    cems["unitid"] = cems["unitid"].str.lstrip("0")

    # crosswalk the plant IDs and add a plant_id_eia column
    cems = crosswalk_epa_eia_plant_ids(cems, year)

    # fill any missing values for operating time or steam load with zero
    cems["operating_time_hours"] = cems["operating_time_hours"].fillna(0)
    cems["steam_load_1000_lb"] = cems["steam_load_1000_lb"].fillna(0)

    # NOTE: The co2 and heat content data are reported as rates (e.g. tons/hr) rather than absolutes
    # Thus they need to be multiplied by operating_time_hours
    # See https://github.com/catalyst-cooperative/pudl/issues/1581
    # calculate gross generation by multiplying gross_load_mw by operating_time_hours
    for col in [
        "gross_generation_mwh",
        "co2_mass_tons",
        "fuel_consumed_mmbtu",
        "nox_mass_lb",
        "so2_mass_lb",
    ]:
        cems[col] = cems[col] * cems["operating_time_hours"]

    # convert co2 mass in tons to lb
    cems["co2_mass_lb"] = cems["co2_mass_tons"] * 2000

    # re-order columns
    cems = cems[
        [
            "plant_id_eia",
            "unitid",
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
            "unitid": "str",
            "co2_mass_measurement_code": "category",
            "nox_mass_measurement_code": "category",
            "so2_mass_measurement_code": "category",
        }
    )

    return cems


def crosswalk_epa_eia_plant_ids(cems, year):
    """
    Adds a column to the CEMS data that matches the EPA plant ID to the EIA plant ID
    Inputs:
        cems: pandas dataframe with hourly emissions data and columns for "plant_id_epa" and "unitid"
    Returns:
        cems: pandas dataframe with an additional column for "plant_id_eia"
    """

    psdc = load_epa_eia_crosswalk(year)

    # create a table that matches EPA plant and unit IDs to an EIA plant ID
    plant_id_crosswalk = psdc[
        ["plant_id_epa", "unitid", "plant_id_eia", "generator_id"]
    ].drop_duplicates()

    # only keep plant ids where the two are different
    plant_id_crosswalk = plant_id_crosswalk[
        plant_id_crosswalk["plant_id_epa"] != plant_id_crosswalk["plant_id_eia"]
    ].dropna()

    # match plant_id_eia on plant_id_epa and unitid
    cems = cems.merge(plant_id_crosswalk, how="left", on=["plant_id_epa", "unitid"])

    # if the merge resulted in any missing plant_id associations, fill with the plant_id_epa, assuming that they are the same
    cems["plant_id_eia"] = cems["plant_id_eia"].fillna(cems["plant_id_epa"])

    # change the id column from float dtype to int
    cems["plant_id_eia"] = cems["plant_id_eia"].astype(int)

    return cems


def load_pudl_table(table_name, year=None):
    """
    Loads a table from the PUDL SQL database.
    Inputs:
        sql_query: string containing table name (to retrieve entire table) or SQL query (to select part of table)
    Returns:
        table: pandas dataframe containing requested query
    """
    # specify the relative path to the sqllite database, and create an sqalchemy engine
    pudl_db = "sqlite:///../data/downloads/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)

    if year is not None:
        sql_query = f"SELECT * FROM {table_name} WHERE report_date >= '{year}-01-01' AND report_date <= '{year}-12-01'"
    else:
        sql_query = table_name

    table = pd.read_sql(sql_query, pudl_engine)

    return table


def load_ghg_emission_factors():
    """
    Read in the table of emissions factors and convert to lb/mmbtu
    """

    efs = pd.read_csv(
        "../data/manual/egrid_static_tables/table_C1_emission_factors_for_CO2_CH4_N2O.csv",
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
        "../data/manual/egrid_static_tables/table_C2_emission_factors_for_NOx.csv",
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
        "../data/manual/egrid_static_tables/table_C3_emission_factors_for_SO2.csv",
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


def initialize_pudl_out(year=None):
    """
    Initializes a `pudl_out` object used to create tables for EIA and FERC Form 1 analysis.

    If `year` is set to `None`, all years of data are returned.
    """
    pudl_db = "sqlite:///../data/downloads/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)

    if year is None:
        pudl_out = pudl.output.pudltabl.PudlTabl(pudl_engine)
    else:

        pudl_out = pudl.output.pudltabl.PudlTabl(
            pudl_engine,
            freq="MS",
            start_date=f"{year}-01-01",
            end_date=f"{year}-12-31",
            fill_tech_desc=False,
        )
    return pudl_out


def load_epa_eia_crosswalk(year):
    """
    Read in the manual EPA-EIA Crosswalk table.
    """
    """
    map_eia_epa_file = importlib.resources.open_binary(
        'pudl.package_data.glue', 'eia_epa_id_crosswalk.csv')

    return pd.read_csv(
        map_eia_epa_file,
        usecols=['plant_id_epa', 'plant_id_eia', 'unitid',
                 'generator_id', 'boiler_id', 'energy_source_code'],
        dtype={'plant_id_epa': 'int32', 'plant_id_eia': 'int32'})"""

    crosswalk = pd.read_csv(
        "../data/downloads/epa/epa_eia_crosswalk.csv",
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

    # remove leading zeros from the generator id and unitid
    crosswalk["EIA_GENERATOR_ID"] = crosswalk["EIA_GENERATOR_ID"].str.lstrip("0")
    crosswalk["CAMD_UNIT_ID"] = crosswalk["CAMD_UNIT_ID"].str.lstrip("0")

    # change the id to an int
    # NOTE: because of NA values, the column is a float, and cannot be
    # converted to int unless the NAs are removed
    # crosswalk["EIA_PLANT_ID"] = crosswalk["EIA_PLANT_ID"].astype(int)

    # rename the columns
    crosswalk = crosswalk.rename(
        columns={
            "CAMD_PLANT_ID": "plant_id_epa",
            "EIA_PLANT_ID": "plant_id_eia",
            "CAMD_UNIT_ID": "unitid",
            "EIA_GENERATOR_ID": "generator_id",
            "EIA_BOILER_ID": "boiler_id",
            "EIA_FUEL_TYPE": "energy_source_code_eia",
            "CAMD_FUEL_TYPE": "energy_source_code_epa",
        }
    )

    camd_to_eia_fuel_type = {
        "Pipeline Natural Gas": "NG",
        "Coal": "SUB",  # assume that generic coal is subbituminous to be conservative
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

    # **** manual adjustments ****
    # The EPA's crosswalk document incorrectly maps plant_id_epa 55248 to plant_id_eia 55248
    # the correct plant_id_eia is 2847
    crosswalk.loc[crosswalk["plant_id_epa"] == 55248, "plant_id_eia"] = 2847

    # load manually inputted data
    crosswalk_manual = pd.read_csv(
        "../data/manual/epa_eia_crosswalk_manual.csv", dtype=get_dtypes()
    ).drop(columns=["notes"])

    # merge the energy source code from EIA-860
    pudl_out = initialize_pudl_out(year=year)
    crosswalk_manual = crosswalk_manual.merge(
        pudl_out.gens_eia860()[
            ["plant_id_eia", "generator_id", "energy_source_code_1"]
        ],
        how="left",
        on=["plant_id_eia", "generator_id"],
    ).rename(columns={"energy_source_code_1": "energy_source_code_eia"})

    # concat this data with the main table
    crosswalk = pd.concat([crosswalk, crosswalk_manual], axis=0,)

    return crosswalk


def load_gross_to_net_data(
    level, conversion_type, threshold_column, lower_threshold, upper_threshold
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
        f"../data/outputs/gross_to_net/{level}_gross_to_net_{conversion_type}.csv",
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
            "../data/outputs/subplant_crosswalk.csv", dtype=get_dtypes()
        )
        subplant_crosswalk = subplant_crosswalk[
            ["plant_id_eia", "unitid", "subplant_id"]
        ].drop_duplicates()

        if level == "plant":
            groupby_columns = ["plant_id_eia"]
        elif level == "subplant":
            groupby_columns = ["plant_id_eia", "subplant_id"]
        subplant_crosswalk = (
            subplant_crosswalk.groupby(groupby_columns, dropna=False)
            .count()
            .reset_index()
            .rename(columns={"unitid": f"units_in_{level}"})
        )

        gtn_data = gtn_data.merge(
            subplant_crosswalk, how="left", on=groupby_columns, validate="many_to_one"
        )

        # divide the intercept by the number of units in each subplant to evenly distribute this to each unit
        gtn_data["intercept"] = gtn_data["intercept"] / gtn_data[f"units_in_{level}"]

    # make sure the report date column is a datetime if loading ratios
    if conversion_type == "ratio":
        gtn_data["report_date"] = pd.to_datetime(gtn_data["report_date"])

    return gtn_data


def load_ipcc_gwp():
    """Load a table containing global warming potential (GWP) values for CO2, CH4, and N2O."""
    return pd.read_csv(
        "../data/manual/ipcc_gwp.csv", index_col="ipcc_version", dtype=get_dtypes()
    )


def load_raw_eia930_data(year, description):

    eia_930 = pd.concat(
        [
            pd.read_csv(
                f"../data/downloads/eia930/EIA930_{description}_{year}_Jan_Jun.csv",
                thousands=",",
                parse_dates=["UTC Time at End of Hour"],
            ),
            pd.read_csv(
                f"../data/downloads/eia930/EIA930_{description}_{year}_Jul_Dec.csv",
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

    # TODO re-localize the timezones for the BAs that report in a different timezone
    # ba_reference = load_ba_reference()
    # bas_to_convert_tz = list(ba_reference.loc[ba_reference.timezone_reporting_eia930 != ba_reference.timezone_local, 'ba_code'])

    return eia_930


def load_ba_reference():
    return pd.read_csv(
        "../data/manual/ba_reference.csv",
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
    dibas = dibas.merge(ba_tz, how="left", on="ba_code")
    dibas = dibas.merge(
        ba_tz,
        how="left",
        left_on="diba_code",
        right_on="ba_code",
        suffixes=(None, "_diba"),
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
        "../data/manual/ba_reference.csv", usecols=["ba_code", f"timezone_{type}"]
    )
    tz = tz.loc[tz["ba_code"] == ba, f"timezone_{type}"]

    if len(tz) == 0:
        raise UserWarning(
            f"The BA {ba} does not have a timezone specified in data/manual/ba_reference.csv. Please add."
        )
    else:
        tz = tz.item()

    return tz
