import gzip
import os
import requests
import shutil
import tarfile
import zipfile

from filepaths import *


def download_pudl_data(zenodo_url):
    """
    Downloads the archived PUDL data release. The most recent version can be found at https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#zenodo-archives
    Inputs:
        zenodo_url: the url to the .tgz file hosted on zenodo
    """

    # get the version number
    pudl_version = zenodo_url.split("/")[-1].replace(".tgz", "")

    # if the pudl data already exists, do not re-download
    if os.path.exists(f"{downloads_folder()}pudl"):
        with open(f"{downloads_folder()}pudl/pudl_version.txt", "r") as f:
            existing_version = f.readlines()[0].replace('\n', '')
        if pudl_version == existing_version:
            print("    PUDL data already downloaded")
        else:
            print("    Downloading new version of pudl")
            shutil.rmtree(f"{downloads_folder()}pudl")
            download_pudl(zenodo_url, pudl_version)
            download_updated_pudl_database(download=True)
    else:
        download_pudl(zenodo_url, pudl_version)
        download_updated_pudl_database(download=True)


def download_pudl(zenodo_url, pudl_version):
    r = requests.get(zenodo_url, params={"download": "1"}, stream=True)
    # specify parameters for progress bar
    total_size_in_bytes = int(r.headers.get("content-length", 0))
    block_size = 1024 * 1024 * 10  # 10 MB
    downloaded = 0
    with open(f"{downloads_folder()}pudl.tgz", "wb") as fd:
        for chunk in r.iter_content(chunk_size=block_size):
            print(
                f"    Downloading PUDL. Progress: {(round(downloaded/total_size_in_bytes*100,2))}%   \r",
                end="",
            )
            fd.write(chunk)
            downloaded += block_size

    # extract the tgz file
    print("    Extracting PUDL data...")
    with tarfile.open(f"{downloads_folder()}pudl.tgz") as tar:
        tar.extractall(f"{data_folder()}")

    # rename the extracted directory to pudl so that we don't have to update this for future versions
    os.rename(f"{data_folder()}{pudl_version}", f"{downloads_folder()}pudl")

    # add a version file
    with open(f"{downloads_folder()}pudl/pudl_version.txt", "w+") as v:
        v.write(pudl_version)

    # delete the downloaded tgz file
    os.remove(f"{downloads_folder()}pudl.tgz")

    print("    PUDL download complete")


def download_updated_pudl_database(download=True):
    """
    Downloaded the updated `pudl.sqlite` file from datasette.

    This is temporary until a new version of the data is published on zenodo
    """
    if download is True:
        print("    Downloading updated pudl.sqlite from Datasette...")
        # remove the existing file from zenodo
        os.remove(f"{downloads_folder()}pudl/pudl_data/sqlite/pudl.sqlite")

        r = requests.get("https://data.catalyst.coop/pudl.db", stream=True)
        with open(f"{downloads_folder()}pudl/pudl_data/sqlite/pudl.sqlite", "wb") as fd:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                fd.write(chunk)

    else:
        pass


def download_chalendar_files():
    """
    download_chalendar_files
    Download raw and cleaned files. Eventually we'll do our own processing to get our own version of chalendar,
    but still will be useful to use this raw file and compare to this cleaned file.

    TODO: download functions share a lot of code, could refactor
    """
    # if there is not yet a directory for egrid, make it
    os.makedirs(f"{downloads_folder()}eia930/chalendar", exist_ok=True)

    # download the cleaned and raw files
    urls = [
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_elec.csv.gz",
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_raw.csv.gz",
    ]
    for url in urls:
        filename = url.split("/")[-1].replace(".gz", "")
        # if the file already exists, do not re-download it
        if os.path.exists(f"{downloads_folder()}eia930/chalendar/{filename}"):
            print(f"    {filename} already downloaded")
        else:
            print(f"    Downloading {filename}")
            r = requests.get(url, stream=True)

            with open(f"{downloads_folder()}eia930/chalendar/{filename}.gz", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)

            # Unzip
            with gzip.open(
                f"{downloads_folder()}eia930/chalendar/{filename}.gz", "rb"
            ) as f_in:
                with open(
                    f"{downloads_folder()}eia930/chalendar/{filename}", "wb"
                ) as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f"{downloads_folder()}eia930/chalendar/{filename}.gz")


def download_egrid_files(egrid_files_to_download):
    """
    Downloads the egrid excel files
    Inputs:
        egrid_files_to_download: a list of urls for the egrid excel files that you want to download
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists(f"{downloads_folder()}egrid"):
        os.mkdir(f"{downloads_folder()}egrid")

    # download the egrid files
    for url in egrid_files_to_download:
        filename = url.split("/")[-1]
        # if the file already exists, do not re-download it
        if os.path.exists(f"{downloads_folder()}egrid/{filename}"):
            print(f"    {filename} already downloaded")
        else:
            print(f"    Downloading {filename}")
            r = requests.get(url, stream=True)

            with open(f"{downloads_folder()}egrid/{filename}", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)


def download_eia930_data(years_to_download):
    """
    Downloads the six month csv files from the EIA-930 website
    Inputs:
        years_to_download: list of four-digit year numbers to download from EIA-930
    """
    # if there is not yet a directory for EIA-930, make it
    if not os.path.exists(f"{downloads_folder()}eia930"):
        os.mkdir(f"{downloads_folder()}eia930")

    # download the egrid files
    for year in years_to_download:
        for description in ["BALANCE", "INTERCHANGE"]:
            for months in ["Jan_Jun", "Jul_Dec"]:
                if os.path.exists(
                    f"{downloads_folder()}eia930/EIA930_{description}_{year}_{months}.csv"
                ):
                    print(f"    {description}_{year}_{months} data already downloaded")
                else:
                    print(f"    downloading {description}_{year}_{months} data")
                    r = requests.get(
                        f"https://www.eia.gov/electricity/gridmonitor/sixMonthFiles/EIA930_{description}_{year}_{months}.csv",
                        stream=True,
                    )

                    with open(
                        f"{downloads_folder()}eia930/EIA930_{description}_{year}_{months}.csv",
                        "wb",
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
    if not os.path.exists(f"{downloads_folder()}epa"):
        os.mkdir(f"{downloads_folder()}epa")

    filename = psdc_url.split("/")[-1]
    # if the file already exists, do not re-download it
    if os.path.exists(f"{downloads_folder()}epa/{filename}"):
        print(f"    {filename} already downloaded")
    else:
        print(f"    Downloading {filename}")
        r = requests.get(psdc_url, stream=True)

        with open(f"{downloads_folder()}epa/{filename}", "wb") as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)


def download_raw_eia923(year):
    """
    Downloads the egrid excel files
    Inputs:
        egrid_files_to_download: a list of urls for the egrid excel files that you want to download
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists(f"{downloads_folder()}eia923"):
        os.mkdir(f"{downloads_folder()}eia923")

    url = f"https://www.eia.gov/electricity/data/eia923/archive/xls/f923_{year}.zip"

    filename = url.split("/")[-1].split(".")[0]
    # if the file already exists, do not re-download it
    if os.path.exists(f"{downloads_folder()}eia923/{filename}"):
        print(f"    {year} EIA-923 already downloaded")
    else:
        print(f"    Downloading {year} EIA-923 data")
        r = requests.get(url, stream=True)

        with open(f"{downloads_folder()}eia923/{filename}.zip", "wb") as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)

        # Unzip
        with zipfile.ZipFile(
            f"{downloads_folder()}eia923/{filename}.zip", "r"
        ) as zip_to_unzip:
            zip_to_unzip.extractall(f"{downloads_folder()}eia923/{filename}")
        os.remove(f"{downloads_folder()}eia923/{filename}.zip")


def download_raw_eia860(year):
    """
    Downloads the egrid excel files
    Inputs:
        egrid_files_to_download: a list of urls for the egrid excel files that you want to download
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists(f"{downloads_folder()}eia860"):
        os.mkdir(f"{downloads_folder()}eia860")

    url = f"https://www.eia.gov/electricity/data/eia860/xls/eia860{year}.zip"
    archive_url = (
        f"https://www.eia.gov/electricity/data/eia860/archive/xls/eia860{year}.zip"
    )

    filename = url.split("/")[-1].split(".")[0]
    # if the file already exists, do not re-download it
    if os.path.exists(f"{downloads_folder()}eia860/{filename}"):
        print(f"    {year} EIA-860 already downloaded")
    else:
        print(f"    Downloading {year} EIA-860 data")
        try:
            r = requests.get(url, stream=True)

            with open(f"{downloads_folder()}eia860/{filename}.zip", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)

            # Unzip
            with zipfile.ZipFile(
                f"{downloads_folder()}eia860/{filename}.zip", "r"
            ) as zip_to_unzip:
                zip_to_unzip.extractall(f"{downloads_folder()}eia860/{filename}")
            os.remove(f"{downloads_folder()}eia860/{filename}.zip")
        except:
            r = requests.get(archive_url, stream=True)

            with open(f"{downloads_folder()}eia860/{filename}.zip", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)

            # Unzip
            with zipfile.ZipFile(
                f"{downloads_folder()}eia860/{filename}.zip", "r"
            ) as zip_to_unzip:
                zip_to_unzip.extractall(f"{downloads_folder()}eia860/{filename}")
            os.remove(f"{downloads_folder()}eia860/{filename}.zip")


def format_raw_eia860(year: int):
    """
    Makes sure that a folder of raw EIA-860 data has filenames that are consistent
    with the rest of the data pipeline.
    """
    if year < 2005:
        raise NotImplementedError(f'We haven\'t implemented data cleaning for {year} yet.')

    raw_folder = f"{downloads_folder()}eia860/eia860{year}"

    # For 2009-2012, the filenames for EnviroEquip and EnviroAssoc are different.
    if year >= 2009 and year <= 2012:
        # 2009 uses two digit years.
        year_format = str(year)[-2:] if year == 2009 else str(year)
        file_format = 'xls' if year <= 2010 else 'xlsx'

        # For some dumb reason, the 2011 data release doesn't include a year in the filename.
        equip_name = 'EnviroEquip.{}'.format(file_format) if year == 2011 else 'EnviroEquipY{}.{}'.format(year_format, file_format)
        enviro_assoc_input_filename = os.path.join(raw_folder, 'EnviroAssocY{}.{}'.format(year_format, file_format))
        enviro_equip_input_filename = os.path.join(raw_folder, equip_name)
        if not os.path.exists(enviro_assoc_input_filename):
            raise FileNotFoundError(enviro_assoc_input_filename)
        if not os.path.exists(enviro_equip_input_filename):
            raise FileNotFoundError(enviro_equip_input_filename)
        
        # Copy to a new file with our expected naming convention.
        shutil.copy(enviro_assoc_input_filename,
            os.path.join(raw_folder, InputDataFilenames.EIA_860_ENVIRO_ASSOC_FILE_FMT.format(year)))
        shutil.copy(enviro_assoc_input_filename,
            os.path.join(raw_folder, InputDataFilenames.EIA_860_ENVIRO_EQUIP_FILE_FMT.format(year)))


def check_required_files_raw_eia860(year: int):
    """
    Checks that all required files are present in an EIA-860 data folder.
    """
    raw_folder = f"{downloads_folder()}eia860/eia860{year}"
    if not os.path.exists(raw_folder):
        raise FileNotFoundError(f'The supplied EIA-86 data folder \'{raw_folder}\' does not exist.')

    enviro_assoc_file = InputDataFilenames.EIA_860_ENVIRO_ASSOC_FILE_FMT.format(year)
    enviro_equip_file = InputDataFilenames.EIA_860_ENVIRO_EQUIP_FILE_FMT.format(year)
    if not os.path.exists(os.path.join(raw_folder, enviro_assoc_file)):
        raise FileNotFoundError(f'EIA-860 for year {year} is missing {enviro_assoc_file}')
    if not os.path.exists(os.path.join(raw_folder, enviro_equip_file)):
        raise FileNotFoundError(f'EIA-860 for year {year} is missing {enviro_equip_file}')
