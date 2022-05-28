import os
import requests
import tarfile
import pandas as pd
import sqlalchemy as sa

import src.data_cleaning as data_cleaning
import shutil
import gzip

import pudl.output.pudltabl


def download_pudl_data(zenodo_url):
    """
    Downloads the archived PUDL data release. The most recent version can be found at https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#zenodo-archives
    Inputs:
        zenodo_url: the url to the .tgz file hosted on zenodo
    """

    # get the version number
    pudl_version = zenodo_url.split("/")[-1].replace(".tgz", "")

    def download_pudl(zenodo_url, pudl_version):
        r = requests.get(zenodo_url, params={"download": "1"}, stream=True)
        # specify parameters for progress bar
        total_size_in_bytes = int(r.headers.get("content-length", 0))
        block_size = 1024 * 1024 * 10  # 10 MB
        downloaded = 0
        with open("../data/pudl.tgz", "wb") as fd:
            for chunk in r.iter_content(chunk_size=block_size):
                print(
                    f"Downloading PUDL. Progress: {(round(downloaded/total_size_in_bytes*100,2))}%   \r",
                    end="",
                )
                fd.write(chunk)
                downloaded += block_size

        # extract the tgz file
        print("Extracting PUDL data...")
        with tarfile.open("../data/pudl.tgz") as tar:
            tar.extractall("../data/")

        # rename the extracted directory to pudl so that we don't have to update this for future versions
        os.rename(f"../data/{pudl_version}", "../data/pudl")

        # add a version file
        with open("../data/pudl/pudl_version.txt", "w+") as v:
            v.write(pudl_version)

        # delete the downloaded tgz file
        os.remove("../data/pudl.tgz")

        print("PUDL download complete")

    # if the pudl data already exists, do not re-download
    if os.path.exists("../data/pudl"):
        with open("../data/pudl/pudl_version.txt", "r") as f:
            existing_version = f.readlines()[0]
        if pudl_version == existing_version:
            print("PUDL data already downloaded")
        else:
            print("Downloading new version of pudl")
            shutil.rmtree("../data/pudl")
            download_pudl(zenodo_url, pudl_version)
    else:
        download_pudl(zenodo_url, pudl_version)


def download_chalendar_files():
    """
    download_chalendar_files
    Download raw and cleaned files. Eventually we'll do our own processing to get our own version of chalendar,
    but still will be useful to use this raw file and compare to this cleaned file.

    TODO: download functions share a lot of code, could refactor
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists("../data/eia930"):
        os.mkdir("../data/eia930")
    # if there is not a directory for chalendar-formatted files, make it
    if not os.path.exists("../data/eia930/chalendar"):
        os.mkdir("../data/eia930/chalendar")

    # download the cleaned and raw files
    urls = [
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_elec.csv.gz",
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_raw.csv.gz",
    ]
    for url in urls:
        filename = url.split("/")[-1].replace(".gz", "")
        # if the file already exists, do not re-download it
        if os.path.exists(f"../data/eia930/chalendar/{filename}"):
            print(f"{filename} already downloaded")
        else:
            r = requests.get(url, stream=True)

            with open(f"../data/eia930/chalendar/{filename}.gz", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)

            # Unzip
            with gzip.open(f"../data/eia930/chalendar/{filename}.gz", "rb") as f_in:
                with open(f"../data/eia930/chalendar/{filename}", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f"../data/eia930/chalendar/{filename}.gz")


def download_egrid_files(egrid_files_to_download):
    """
    Downloads the egrid excel files
    Inputs:
        egrid_files_to_download: a list of urls for the egrid excel files that you want to download
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists("../data/egrid"):
        os.mkdir("../data/egrid")

    # download the egrid files
    for url in egrid_files_to_download:
        filename = url.split("/")[-1]
        # if the file already exists, do not re-download it
        if os.path.exists(f"../data/egrid/{filename}"):
            print(f"{filename} already downloaded")
        else:
            r = requests.get(url, stream=True)

            with open(f"../data/egrid/{filename}", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)


def download_eia930_data(years_to_download):
    """
    Downloads the six month csv files from the EIA-930 website
    Inputs:
        years_to_download: list of four-digit year numbers to download from EIA-930
    """
    # if there is not yet a directory for EIA-930, make it
    if not os.path.exists("../data/eia930"):
        os.mkdir("../data/eia930")

    # download the egrid files
    for year in years_to_download:
        for period in ["Jan_Jun", "Jul_Dec"]:
            if os.path.exists(f"../data/eia930/EIA930_BALANCE_{year}_{period}.csv"):
                print(f"{year}_{period} data already downloaded")
            else:
                print(f"downloading {year}_{period} data")
                r = requests.get(
                    f"https://www.eia.gov/electricity/gridmonitor/sixMonthFiles/EIA930_BALANCE_{year}_{period}.csv",
                    stream=True,
                )

                with open(
                    f"../data/eia930/EIA930_BALANCE_{year}_{period}.csv", "wb"
                ) as fd:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        fd.write(chunk)


def download_epa_psdc(psdc_url):
    """
    Downloads the EPA's Power Sector Data Crosswalk
    Check for new releases at https://github.com/USEPA/camd-eia-crosswalk
    Inputs:
        psdc_url: the url to the csv file hosted on github
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists("../data/epa"):
        os.mkdir("../data/epa")

    filename = psdc_url.split("/")[-1]
    # if the file already exists, do not re-download it
    if os.path.exists(f"../data/epa/{filename}"):
        print(f"{filename} already downloaded")
    else:
        r = requests.get(psdc_url, stream=True)

        with open(f"../data/epa/{filename}", "wb") as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)


def load_cems_data(year):
    """
    Loads CEMS data for the specified year from the PUDL database
    Inputs:
        year: the year for which data should be retrieved (YYYY)
    Returns:
        cems: pandas dataframe with hourly CEMS data
    """
    # specify the path to the CEMS data
    cems_path = f"../data/pudl/pudl_data/parquet/epacems/year={year}"

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
        "heat_content_mmbtu",
        "unit_id_epa",
    ]

    # load the CEMS data
    cems = pd.read_parquet(cems_path, columns=cems_columns)

    # rename cems plant_id_eia to plant_id_epa (PUDL simply renames the ORISPL_CODE column from the raw CEMS data as 'plant_id_eia' without actually crosswalking to the EIA id)
    # rename the heat content column to use the convention used in the EIA data
    cems = cems.rename(
        columns={
            "plant_id_eia": "plant_id_epa",
            "heat_content_mmbtu": "fuel_consumed_mmbtu",
        }
    )

    # if the unitid has any leading zeros, remove them
    cems["unitid"] = cems["unitid"].str.lstrip("0")

    # crosswalk the plant IDs and add a plant_id_eia column
    cems = data_cleaning.crosswalk_epa_eia_plant_ids(cems, year)

    # fill any missing values for operating time or steam load with zero
    cems["operating_time_hours"] = cems["operating_time_hours"].fillna(0)
    cems["steam_load_1000_lbs"] = cems["steam_load_1000_lbs"].fillna(0)

    # NOTE: The co2 and heat content data are reported as rates (e.g. tons/hr) rather than absolutes
    # Thus they need to be multiplied by operating_time_hours
    # See https://github.com/catalyst-cooperative/pudl/issues/1581
    # calculate gross generation by multiplying gross_load_mw by operating_time_hours
    cems["gross_generation_mwh"] = cems["gross_load_mw"] * cems["operating_time_hours"]
    cems["co2_mass_tons"] = cems["co2_mass_tons"] * cems["operating_time_hours"]
    cems["fuel_consumed_mmbtu"] = (
        cems["fuel_consumed_mmbtu"] * cems["operating_time_hours"]
    )

    # add a unique unit id
    cems["cems_id"] = (
        cems["plant_id_eia"].astype(str) + "_" + cems["unitid"].astype(str)
    )

    # re-order columns
    cems = cems[
        [
            "plant_id_eia",
            "unitid",
            "cems_id",
            "operating_datetime_utc",
            "operating_time_hours",
            "gross_load_mw",
            "gross_generation_mwh",
            "steam_load_1000_lbs",
            "fuel_consumed_mmbtu",
            "co2_mass_tons",
            "co2_mass_measurement_code",
            "plant_id_epa",
            "unit_id_epa",
        ]
    ]

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
    pudl_db = "sqlite:///../data/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)

    if year is not None:
        sql_query = f"SELECT * FROM {table_name} WHERE report_date >= '{year}-01-01' AND report_date <= '{year}-12-01'"
    else:
        sql_query = table_name

    table = pd.read_sql(sql_query, pudl_engine)

    return table


def load_emission_factors():
    """
    Read in the table of emissions factors
    """
    return pd.read_csv(
        "../data/egrid/egrid_static_tables/table_C1_emission_factors_for_CO2_CH4_N2O.csv"
    )


def load_emission_factors_nox():
    """Read in the NOx emission factors from eGRID Table C2."""
    return pd.read_csv(
        "../data/egrid/egrid_static_tables/table_C2_emission_factors_for_NOx.csv"
    )


def load_emission_factors_so2():
    """
    Read in the SO2 emission factors from eGRID Table C3.

    The SO2 emission rate depends on the sulfur content of fuel, so it is
    reported in Table C3 as a formula like `123*S`.
    """
    df = pd.read_csv( "../data/egrid/egrid_static_tables/table_C3_emission_factors_for_SO2.csv")

    # Add a boolean column that reports whether the emission factor is a formula or value.
    df['multiply_by_sulfur_content'] = df['Emission Factor'].str.contains('*', regex=False).astype(int)

    # Extract the numeric coefficient from the emission factor.
    df['emission_factor_coeff'] = df['Emission Factor'].str.replace('*S', '', regex=False).astype(float)

    return df


def initialize_pudl_out(year=None):
    """
    Initializes a `pudl_out` object used to create tables for EIA and FERC Form 1 analysis.

    If `year` is set to `None`, all years of data are returned.
    """
    pudl_db = "sqlite:///../data/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)

    if year is None:
        pudl_out = pudl.output.pudltabl.PudlTabl(pudl_engine)
    else:

        pudl_out = pudl.output.pudltabl.PudlTabl(
            pudl_engine, freq="MS", start_date=f"{year}-01-01", end_date=f"{year}-12-31", fill_tech_desc=False
        )
    return pudl_out


def create_flat_eia930_placeholder_profiles(
    monthly_eia_data_to_distribute, energy_source_groups, year
):
    # For now, create a synthetic flat profile for each resource for now until we have shapes to distribute to from EIA-930

    # create lists of all bas and fuel types
    ba_list = list(monthly_eia_data_to_distribute["ba_code"].dropna().unique())
    fuel_list = list(energy_source_groups["fuel_category"].unique())

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles = []

    for ba in ba_list:
        for fuel in fuel_list:
            # create a dataframe
            df_temp = pd.DataFrame(
                index=pd.date_range(
                    start=f"{year}-01-01 00:00:00",
                    end=f"{year}-12-31 23:00:00",
                    freq="H",
                    tz=data_cleaning.ba_timezone(ba=ba, type="local"),
                    name="datetime_local",
                ),
                columns=["ba_code", "fuel_category", "profile"],
            ).reset_index()
            df_temp["ba_code"] = ba
            df_temp["fuel_category"] = fuel
            df_temp["profile"] = 1.0
            df_temp["report_date"] = df_temp["datetime_local"].astype(str).str[:7]
            df_temp["report_date"] = pd.to_datetime(df_temp["report_date"])
            hourly_profiles.append(df_temp)

    hourly_profiles = pd.concat(hourly_profiles, axis=0, ignore_index=True)


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
        "../data/epa/epa_eia_crosswalk.csv",
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
        "../data/manual/epa_eia_crosswalk_manual.csv", dtype={"generator_id": str}
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
