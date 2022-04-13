import os
import requests
import tarfile
import pandas as pd
import sqlalchemy as sa

import src.data_cleaning as data_cleaning


def download_pudl_data(zenodo_url):
    """
    Downloads the archived PUDL data release. The most recent version can be found at https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#zenodo-archives
    Inputs:
        zenodo_url: the url to the .tgz file hosted on zenodo
    """
    zenodo_url = 'https://zenodo.org/record/5701406/files/pudl-v0.5.0-2021-11-14.tgz'
    # get the version number
    pudl_version = zenodo_url.split('/')[-1].replace('.tgz','')

    # if the pudl data already exists, do not re-download
    if os.path.exists(f'../data/pudl'):
        print('PUDL data already downloaded')
    else:
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
        os.rename(f'../data/{pudl_version}', 'pudl')

        # delete the downloaded tgz file
        os.remove("../data/pudl.tgz")

        print('PUDL download complete')

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

def load_pudl_table(sql_query):
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

    table = pd.read_sql(sql_query, pudl_engine)

    return table

def load_emission_factors():
    """
    Read in the table of emissions factors
    """
    return pd.read_csv('../data/egrid/egrid_static_tables/table_C1_emission_factors_for_CO2_CH4_N2O.csv')