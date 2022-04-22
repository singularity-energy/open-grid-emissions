import os
import requests
import tarfile
import pandas as pd
import sqlalchemy as sa
import time

import src.data_cleaning as data_cleaning
import shutil


def download_pudl_data(zenodo_url):
    """
    Downloads the archived PUDL data release. The most recent version can be found at https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#zenodo-archives
    Inputs:
        zenodo_url: the url to the .tgz file hosted on zenodo
    """

    # get the version number
    pudl_version = zenodo_url.split('/')[-1].replace('.tgz','')

    def download_pudl(zenodo_url, pudl_version):
        r = requests.get(zenodo_url, params={"download":"1"}, stream=True)
        # specify parameters for progress bar
        total_size_in_bytes= int(r.headers.get('content-length', 0))
        block_size = 1024 * 1024 * 10 # 10 MB
        downloaded = 0
        with open("../data/pudl.tgz", 'wb') as fd:
            for chunk in r.iter_content(chunk_size=block_size):
                print(f'Downloading PUDL. Progress: {(round(downloaded/total_size_in_bytes*100,2))}%   \r', end='')
                fd.write(chunk)
                downloaded += block_size

        # extract the tgz file
        print('Extracting PUDL data...')
        with tarfile.open("../data/pudl.tgz") as tar:
            tar.extractall('../data/')

        # rename the extracted directory to pudl so that we don't have to update this for future versions
        os.rename(f'../data/{pudl_version}', '../data/pudl')

        # add a version file
        with open('../data/pudl/pudl_version.txt', 'w+') as v:
            v.write(pudl_version)

        # delete the downloaded tgz file
        os.remove("../data/pudl.tgz")

        print('PUDL download complete')

    # if the pudl data already exists, do not re-download
    if os.path.exists(f'../data/pudl'):
        with open('../data/pudl/pudl_version.txt', 'r') as f:
            existing_version = f.readlines()[0]
        if pudl_version == existing_version:
            print('PUDL data already downloaded')
        else:
            print('Downloading new version of pudl')
            shutil.rmtree('../data/pudl')
            download_pudl(zenodo_url, pudl_version)
    else:
        download_pudl(zenodo_url, pudl_version)
        

def download_egrid_files(egrid_files_to_download):
    """
    Downloads the egrid excel files
    Inputs: 
        egrid_files_to_download: a list of urls for the egrid excel files that you want to download
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists('../data/egrid'):
        os.mkdir('../data/egrid')

    # download the egrid files
    for url in egrid_files_to_download:
        filename = url.split("/")[-1]
        # if the file already exists, do not re-download it
        if os.path.exists(f'../data/egrid/{filename}'):
            print(f'{filename} already downloaded')
        else:
            r = requests.get(url, stream=True)
            
            with open(f'../data/egrid/{filename}', 'wb') as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)

def download_eia930_data(years_to_download):
    """
    Downloads the six month csv files from the EIA-930 website
    Inputs:
        years_to_download: list of four-digit year numbers to download from EIA-930
    """
    # if there is not yet a directory for EIA-930, make it
    if not os.path.exists('../data/eia930'):
        os.mkdir('../data/eia930')

    # download the egrid files
    for year in years_to_download:
        for period in ['Jan_Jun','Jul_Dec']:
            if os.path.exists(f'../data/eia930/EIA930_BALANCE_{year}_{period}.csv'):
                print(f'{year}_{period} data already downloaded')
            else:
                print(f'downloading {year}_{period} data')
                r = requests.get(f"https://www.eia.gov/electricity/gridmonitor/sixMonthFiles/EIA930_BALANCE_{year}_{period}.csv", stream=True)
            
                with open(f'../data/eia930/EIA930_BALANCE_{year}_{period}.csv', 'wb') as fd:
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
    if not os.path.exists('../data/epa'):
        os.mkdir('../data/epa')

    filename = psdc_url.split("/")[-1]
    # if the file already exists, do not re-download it
    if os.path.exists(f'../data/epa/{filename}'):
        print(f'{filename} already downloaded')
    else:
        r = requests.get(psdc_url, stream=True)
        
        with open(f'../data/epa/{filename}', 'wb') as fd:
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
    #specify the path to the CEMS data
    cems_path = f'../data/pudl/pudl_data/parquet/epacems/year={year}' 

    # specify the columns to use from the CEMS database
    cems_columns = ['plant_id_eia', 'unitid', 'operating_datetime_utc',
    'operating_time_hours', 'gross_load_mw', 'steam_load_1000_lbs',
    'co2_mass_tons', 'co2_mass_measurement_code', 'heat_content_mmbtu',
    'unit_id_epa']

    # load the CEMS data
    cems = pd.read_parquet(cems_path, columns=cems_columns)

    # rename cems plant_id_eia to plant_id_epa (PUDL simply renames the ORISPL_CODE column from the raw CEMS data as 'plant_id_eia' without actually crosswalking to the EIA id)
    cems = cems.rename(columns={'plant_id_eia': 'plant_id_epa'})

    # if the unitid has any leading zeros, remove them
    cems['unitid'] = cems['unitid'].str.lstrip('0')

    # crosswalk the plant IDs and add a plant_id_eia column
    cems = data_cleaning.crosswalk_epa_eia_plant_ids(cems)

    # calculate gross generation by multiplying gross_load_mw by operating_time_hours
    cems['gross_generation_mwh'] = cems['gross_load_mw'] * cems['operating_time_hours']

    # NOTE: The co2 and heat content data are reported as rates (e.g. tons/hr) rather than absolutes
    # Thus they need to be multiplied by operating_time_hours
    # See https://github.com/catalyst-cooperative/pudl/issues/1581
    cems['co2_mass_tons'] = cems['co2_mass_tons'] * cems['operating_time_hours']
    cems['heat_content_mmbtu'] = cems['heat_content_mmbtu'] * cems['operating_time_hours']

    # add a unique unit id
    cems['cems_id'] = cems['plant_id_eia'].astype(str) + "_" + cems['unitid'].astype(str)

    # re-order columns
    cems = cems[['plant_id_eia', 'unitid', 'cems_id',  'operating_datetime_utc',
    'operating_time_hours', 'gross_load_mw', 'gross_generation_mwh', 'steam_load_1000_lbs','heat_content_mmbtu',
    'co2_mass_tons', 'co2_mass_measurement_code', 'plant_id_epa','unit_id_epa']]

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
    pudl_db = 'sqlite:///../data/pudl/pudl_data/sqlite/pudl.sqlite'
    pudl_engine = sa.create_engine(pudl_db)

    if year!=None:
        sql_query = f"SELECT * FROM {table_name} WHERE report_date >= '{year}-01-01' AND report_date <= '{year}-12-01'"
    else:
        sql_query = table_name
    
    table = pd.read_sql(sql_query, pudl_engine)
    
    return table

def load_emission_factors():
    """
    Read in the table of emissions factors
    """
    return pd.read_csv('../data/egrid/egrid_static_tables/table_C1_emission_factors_for_CO2_CH4_N2O.csv')

def plants_eia860(pudl_engine, start_date=None, end_date=None):
    """Pull all fields from the EIA Plants tables.
    Args:
        pudl_engine (sqlalchemy.engine.Engine): SQLAlchemy connection engine
            for the PUDL DB.
        start_date (date-like): date-like object, including a string of the
            form 'YYYY-MM-DD' which will be used to specify the date range of
            records to be pulled.  Dates are inclusive.
        end_date (date-like): date-like object, including a string of the
            form 'YYYY-MM-DD' which will be used to specify the date range of
            records to be pulled.  Dates are inclusive.
    Returns:
        pandas.DataFrame: A DataFrame containing all the fields of the EIA 860
        Plants table.
    """
    pt = get_table_meta(pudl_engine)
    # grab the entity table
    plants_eia_tbl = pt["plants_entity_eia"]
    plants_eia_select = sa.sql.select(plants_eia_tbl)
    plants_eia_df = pd.read_sql(plants_eia_select, pudl_engine)

    # grab the annual table select
    plants_eia860_tbl = pt["plants_eia860"]
    plants_eia860_select = sa.sql.select(plants_eia860_tbl)
    if start_date is not None:
        start_date = pd.to_datetime(start_date)
        plants_eia860_select = plants_eia860_select.where(
            plants_eia860_tbl.c.report_date >= start_date
        )
    if end_date is not None:
        end_date = pd.to_datetime(end_date)
        plants_eia860_select = plants_eia860_select.where(
            plants_eia860_tbl.c.report_date <= end_date
        )
    plants_eia860_df = pd.read_sql(plants_eia860_select, pudl_engine).assign(
        report_date=lambda x: pd.to_datetime(x.report_date)
    )

    # plant glue table
    plants_g_eia_tbl = pt["plants_eia"]
    plants_g_eia_select = sa.sql.select(
        plants_g_eia_tbl.c.plant_id_eia,
        plants_g_eia_tbl.c.plant_id_pudl,
    )
    plants_g_eia_df = pd.read_sql(plants_g_eia_select, pudl_engine)

    out_df = pd.merge(plants_eia_df, plants_eia860_df, how="left", on=["plant_id_eia"])
    out_df = pd.merge(out_df, plants_g_eia_df, how="left", on=["plant_id_eia"])

    utils_eia_tbl = pt["utilities_eia"]
    utils_eia_select = sa.sql.select(utils_eia_tbl)
    utils_eia_df = pd.read_sql(utils_eia_select, pudl_engine)

    out_df = (
        pd.merge(out_df, utils_eia_df, how="left", on=["utility_id_eia"])
        .dropna(subset=["report_date", "plant_id_eia"])
        #.pipe(apply_dtype)
    )
    return out_df

def utilities_eia860(pudl_engine, start_date=None, end_date=None):
    """Pull all fields from the EIA860 Utilities table. NOTE: copied from pudl.output.eia860
    Args:
        pudl_engine (sqlalchemy.engine.Engine): SQLAlchemy connection engine
            for the PUDL DB.
        start_date (date-like): date-like object, including a string of the
            form 'YYYY-MM-DD' which will be used to specify the date range of
            records to be pulled.  Dates are inclusive.
        end_date (date-like): date-like object, including a string of the
            form 'YYYY-MM-DD' which will be used to specify the date range of
            records to be pulled.  Dates are inclusive.
    Returns:
        pandas.DataFrame: A DataFrame containing all the fields of the EIA 860
        Utilities table.
    """
    pt = get_table_meta(pudl_engine)
    # grab the entity table
    utils_eia_tbl = pt["utilities_entity_eia"]
    utils_eia_select = sa.sql.select(utils_eia_tbl)
    utils_eia_df = pd.read_sql(utils_eia_select, pudl_engine)

    # grab the annual eia entity table
    utils_eia860_tbl = pt["utilities_eia860"]
    utils_eia860_select = sa.sql.select(utils_eia860_tbl)

    if start_date is not None:
        start_date = pd.to_datetime(start_date)
        utils_eia860_select = utils_eia860_select.where(
            utils_eia860_tbl.c.report_date >= start_date
        )
    if end_date is not None:
        end_date = pd.to_datetime(end_date)
        utils_eia860_select = utils_eia860_select.where(
            utils_eia860_tbl.c.report_date <= end_date
        )
    utils_eia860_df = pd.read_sql(utils_eia860_select, pudl_engine)

    # grab the glue table for the utility_id_pudl
    utils_g_eia_tbl = pt["utilities_eia"]
    utils_g_eia_select = sa.sql.select(
        utils_g_eia_tbl.c.utility_id_eia,
        utils_g_eia_tbl.c.utility_id_pudl,
    )
    utils_g_eia_df = pd.read_sql(utils_g_eia_select, pudl_engine)

    out_df = pd.merge(utils_eia_df, utils_eia860_df, how="left", on=["utility_id_eia"])
    out_df = pd.merge(out_df, utils_g_eia_df, how="left", on=["utility_id_eia"])
    out_df = (
        out_df.assign(report_date=lambda x: pd.to_datetime(x.report_date))
        .dropna(subset=["report_date", "utility_id_eia"])
        #.pipe(apply_dtype)
    )
    first_cols = [
        "report_date",
        "utility_id_eia",
        "utility_id_pudl",
        "utility_name_eia",
    ]

    out_df = organize_cols(out_df, first_cols)
    return out_df

def boiler_generator_assn_eia860(pudl_engine, start_date=None, end_date=None):
    """Pull all fields from the EIA 860 boiler generator association table. NOTE: Copied from pudl.output.eia860
    Args:
        pudl_engine (sqlalchemy.engine.Engine): SQLAlchemy connection engine
            for the PUDL DB.
        start_date (date-like): date-like object, including a string of the
            form 'YYYY-MM-DD' which will be used to specify the date range of
            records to be pulled.  Dates are inclusive.
        end_date (date-like): date-like object, including a string of the
            form 'YYYY-MM-DD' which will be used to specify the date range of
            records to be pulled.  Dates are inclusive.
    Returns:
        pandas.DataFrame: A DataFrame containing all the fields from the EIA
        860 boiler generator association table.
    """
    pt = get_table_meta(pudl_engine)
    bga_eia860_tbl = pt["boiler_generator_assn_eia860"]
    bga_eia860_select = sa.sql.select(bga_eia860_tbl)

    if start_date is not None:
        start_date = pd.to_datetime(start_date)
        bga_eia860_select = bga_eia860_select.where(
            bga_eia860_tbl.c.report_date >= start_date
        )
    if end_date is not None:
        end_date = pd.to_datetime(end_date)
        bga_eia860_select = bga_eia860_select.where(
            bga_eia860_tbl.c.report_date <= end_date
        )
    out_df = pd.read_sql(bga_eia860_select, pudl_engine).assign(
        report_date=lambda x: pd.to_datetime(x.report_date)
    )
    return out_df

def organize_cols(df, cols):
    """
    Organize columns into key ID & name fields & alphabetical data columns.
    For readability, it's nice to group a few key columns at the beginning
    of the dataframe (e.g. report_year or report_date, plant_id...) and then
    put all the rest of the data columns in alphabetical order.
    Args:
        df: The DataFrame to be re-organized.
        cols: The columns to put first, in their desired output ordering.
    Returns:
        pandas.DataFrame: A dataframe with the same columns as the input
        DataFrame df, but with cols first, in the same order as they
        were passed in, and the remaining columns sorted alphabetically.
    """
    # Generate a list of all the columns in the dataframe that are not
    # included in cols
    data_cols = sorted([c for c in df.columns.tolist() if c not in cols])
    organized_cols = cols + data_cols
    return df[organized_cols]

def get_table_meta(pudl_engine):
    """Grab the pudl sqlitie database table metadata."""
    md = sa.MetaData()
    md.reflect(pudl_engine)
    return md.tables