from typing import Optional
import gzip
import os
import requests
import shutil
import tarfile
import zipfile

from filepaths import downloads_folder, data_folder


class InputDataFilenames:
    # EIA-860 filenames.
    EIA_860_ENVIRO_ASSOC_FILE_FMT = '6_1_EnviroAssoc_Y{}.xlsx'  # .format(year)
    EIA_860_ENVIRO_EQUIP_FILE_FMT = '6_2_EnviroEquip_Y{}.xlsx'  # .format(year)

    # EIA-923 filenames.
    EIA_923_ENVIRONMENTAL_INFO_FMT = 'EIA923_Schedule_8_Annual_Environmental_Information_{}_Final_Revision.xlsx'


def download_helper(input_url: str,
                    download_path: str,
                    output_path: Optional[str] = None,
                    requires_unzip: bool = False,
                    requires_untar: bool = False,
                    requires_gzip: bool = False,
                    should_clean: bool = False,
                    chunk_size: int = 1024) -> bool:
    """
    Downloads a file or archive and optionally unzips/untars/copies it to a destination.

    Inputs:
        `input_url`: Where to download data from.
        `download_path`: An absolute filepath to download to.
        `output_path`: The final destination where the downloaded data should end up.
        `requires_unzip`: Should we unzip the file after downloading?
        `requires_untar`: Should we untar the file after downloading?
        `requires_gzip`: Should we un-gzip the file after downloading?
        `should_clean`: Should we delete the temporary downloaded file when finished?
        `chunk_size`: The chunk size for downloading.
    
    Returns:
        (bool) Whether the file was downloaded (it might be skipped if found).
    """
    # If the file already exists, do not re-download it.
    final_destination = output_path if output_path is not None else download_path
    if os.path.exists(final_destination):
        print(f"    {final_destination} already downloaded, skipping.")
        return False

    # Otherwise, download to the file in chunks.
    print(f"    Downloading from {input_url}")
    r = requests.get(input_url, stream=True)
    with open(download_path, "wb") as fd:
        for chunk in r.iter_content(chunk_size=chunk_size):
            fd.write(chunk)
    # Optionally unzip the downloaded file.
    if requires_unzip:
        if output_path is None:
            raise ValueError('Unzipping requires an output_path destination.')
        with zipfile.ZipFile(download_path, "r") as zip_to_unzip:
            zip_to_unzip.extractall(output_path)
    # Optionally un-tar the downloaded file.
    elif requires_untar:
        if output_path is None:
            raise ValueError('Extracting a tar requires an output_path destination.')
        with tarfile.open(download_path) as tar:
            tar.extractall(output_path)
    elif requires_gzip:
        if output_path is None:
            raise ValueError('Extracting a gzip requires an output_path destination.')
        with gzip.open(download_path, "rb") as f_in:
            with open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    # If the user didn't ask for unzip/untar, but specified a different output_path,
    # copy the downloaded file to there.
    elif output_path is not None and output_path != download_path:
        shutil.copy(download_path, output_path)
    # Finally, optionally clean up the downloaded temporary file.
    if should_clean and output_path != download_path:
        os.remove(download_path)
    return True


def download_pudl_data(zenodo_url: str):
    """
    Downloads an archived PUDL data release.
    
    The most recent version can be found at:
    https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#zenodo-archives

    Inputs:
        `zenodo_url`: the url to the .tgz file hosted on zenodo
    """
    # get the version number
    pudl_version = zenodo_url.split("/")[-1].replace(".tgz", "")

    # if the pudl data already exists, do not re-download
    if os.path.exists(downloads_folder("pudl")):
        pudl_version_file = downloads_folder("pudl/pudl_version.txt")
        with open(pudl_version_file, "r") as f:
            existing_version = f.readlines()[0].replace('\n', '')
        if pudl_version == existing_version:
            print("    PUDL version already downloaded")
            return
        else:
            print("    Downloading new version of pudl")
            shutil.rmtree(downloads_folder("pudl"))

    download_pudl(zenodo_url, pudl_version)
    download_updated_pudl_database(download=True)


def download_pudl(zenodo_url, pudl_version):
    r = requests.get(zenodo_url, params={"download": "1"}, stream=True)
    # specify parameters for progress bar
    total_size_in_bytes = int(r.headers.get("content-length", 0))
    block_size = 1024 * 1024 * 10  # 10 MB
    downloaded = 0
    with open(downloads_folder("pudl.tgz"), "wb") as fd:
        for chunk in r.iter_content(chunk_size=block_size):
            print(
                f"    Downloading PUDL. Progress: {(round(downloaded/total_size_in_bytes*100,2))}%   \r",
                end="",
            )
            fd.write(chunk)
            downloaded += block_size

    # extract the tgz file
    print("    Extracting PUDL data...")
    with tarfile.open(downloads_folder("pudl.tgz")) as tar:
        tar.extractall(data_folder())

    # rename the extracted directory to pudl so that we don't have to update this for future versions
    os.rename(data_folder(pudl_version), downloads_folder("pudl"))

    # add a version file
    with open(downloads_folder("pudl/pudl_version.txt"), "w+") as v:
        v.write(pudl_version)

    # delete the downloaded tgz file
    os.remove(downloads_folder("pudl.tgz"))
    print("    PUDL download complete")


def download_updated_pudl_database(download=True):
    """
    Downloaded the updated `pudl.sqlite` file from datasette.

    This is temporary until a new version of the data is published on zenodo.
    """
    if download is True:
        print("    Downloading updated pudl.sqlite from Datasette...")
        # remove the existing file from zenodo
        os.remove(downloads_folder("pudl/pudl_data/sqlite/pudl.sqlite"))

        r = requests.get("https://data.catalyst.coop/pudl.db", stream=True)
        with open(downloads_folder("pudl/pudl_data/sqlite/pudl.sqlite"), "wb") as fd:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                fd.write(chunk)


def download_chalendar_files():
    """
    Download raw and cleaned files. Eventually we'll do our own processing to get our
    own version of chalendar, but still will be useful to use this raw file and compare
    to this cleaned file.
    """
    os.makedirs(downloads_folder("eia930/chalendar"), exist_ok=True)

    # download the cleaned and raw files
    urls = [
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_elec.csv.gz",
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_raw.csv.gz",
    ]
    for url in urls:
        output_filename = url.split("/")[-1].replace(".gz", "")
        output_filepath = f"{downloads_folder()}eia930/chalendar/{output_filename}"
        download_helper(url, output_filepath + ".gz", output_filepath, requires_gzip=True, should_clean=True)


def download_egrid_files(urls_to_download: list[str]):
    """
    Downloads the egrid excel files.

    Inputs:
        `urls_to_download`: a list of urls for the excel files that you want to download
    """
    os.makedirs(downloads_folder("egrid"), exist_ok=True)

    for url in urls_to_download:
        filename = url.split("/")[-1]
        filepath = f"{downloads_folder()}egrid/{filename}"
        download_helper(url, filepath)


def download_eia930_data(years_to_download: list[int]):
    """
    Downloads the six month csv files from the EIA-930 website.

    Inputs:
        `years_to_download`: list of four-digit year numbers to download from EIA-930.
    """
    os.makedirs(downloads_folder("eia930"), exist_ok=True)

    for year in years_to_download:
        for description in ["BALANCE", "INTERCHANGE"]:
            for months in ["Jan_Jun", "Jul_Dec"]:
                download_url = f"https://www.eia.gov/electricity/gridmonitor/sixMonthFiles/EIA930_{description}_{year}_{months}.csv"
                download_filepath = f"{downloads_folder()}eia930/EIA930_{description}_{year}_{months}.csv"
                download_helper(download_url, download_filepath, chunk_size=1024 * 1024)


def download_epa_psdc(psdc_url: str):
    """
    Downloads the EPA's Power Sector Data Crosswalk.

    Check for new releases at https://github.com/USEPA/camd-eia-crosswalk.

    Inputs:
        `psdc_url`: the url to the csv file hosted on github
    """
    os.makedirs(downloads_folder("epa"), exist_ok=True)
    filename = psdc_url.split("/")[-1]
    download_path = downloads_folder(f"epa/{filename}")
    download_helper(psdc_url, download_path)


def download_raw_eia923(year: int):
    """
    Downloads raw EIA-923 data (zip files), and unzips them to the downloads folder.

    Inputs:
        `year`: A four-digit year.
    """
    if year < 2008:
        raise NotImplementedError(f'EIA-923 data is unavailable for \'{year}\'.')
    os.makedirs(downloads_folder("eia923"), exist_ok=True)
    url = f"https://www.eia.gov/electricity/data/eia923/archive/xls/f923_{year}.zip"
    filename = url.split("/")[-1].split(".")[0]
    download_filepath = downloads_folder(f"eia923/{filename}.zip")
    output_filepath = downloads_folder(f"eia923/{filename}")
    download_helper(url, download_filepath, output_filepath, requires_unzip=True, should_clean=True)


def download_raw_eia_906_920(year):
    """
    For years before 2008, the EIA releases Form 906 and 920 instead of 923.

    Inputs:
        `year`: A four-digit year.
    """
    if year < 2005 or year > 2007:
        raise NotImplementedError(f'EIA-906/920 data is unavailable for \'{year}\'.')
    output_folder = f'f906920_{year}'
    download_helper(
        f'https://www.eia.gov/electricity/data/eia923/archive/xls/f906920_{year}.zip',
        downloads_folder(os.path.join('eia923', output_folder + '.zip')),
        downloads_folder(os.path.join('eia923', output_folder)),
        requires_unzip=True, should_clean=True)


def download_raw_eia860(year):
    """
    Downloads raw EIA-860 data (zip files), and unzips them to the downloads folder.
    """
    if year < 2005:
        raise NotImplementedError(f'WARNING: We haven\'t tested EIA-860 for \'{year}\'.')
    os.makedirs(downloads_folder("eia860"), exist_ok=True)
    url = f"https://www.eia.gov/electricity/data/eia860/xls/eia860{year}.zip"
    archive_url = (
        f"https://www.eia.gov/electricity/data/eia860/archive/xls/eia860{year}.zip"
    )
    filename = url.split("/")[-1].split(".")[0]
    zip_filepath = downloads_folder(f"eia860/{filename}.zip")
    output_filepath = downloads_folder(f"eia860/{filename}")
    try:
        download_helper(url, zip_filepath, output_filepath, requires_unzip=True)
    except Exception:
        download_helper(archive_url, zip_filepath, output_filepath, requires_unzip=True)


def download_eia_electric_power_annual():
    """
    Downloads EIA Electric Power Annual uncontrolled emission factors.
    
    See: https://www.eia.gov/electricity/annual/
    """
    os.makedirs(downloads_folder('eia_electric_power_annual'), exist_ok=True)
    urls = [
        "https://www.eia.gov/electricity/annual/xls/epa_a_01.xlsx",
        "https://www.eia.gov/electricity/annual/xls/epa_a_02.xlsx",
        "https://www.eia.gov/electricity/annual/xls/epa_a_03.xlsx",
        "https://www.eia.gov/electricity/annual/xls/epa_a_04.xlsx"
    ]
    output_filenames = [
        'epa_a_01_so2_uncontrolled_efs.xlsx',
        'epa_a_02_nox_uncontrolled_efs.xlsx',
        'epa_a_03_co2_uncontrolled_efs.xlsx',
        'epa_a_04_nox_reduction_factors.xlsx',
    ]
    for url, output_filename in zip(urls, output_filenames):
        output_filepath = downloads_folder(f"eia_electric_power_annual/{output_filename}")
        download_helper(url, output_filepath)


def format_raw_eia860(year: int):
    """
    Makes sure that a folder of raw EIA-860 data has filenames that are consistent
    with the rest of the data pipeline.
    """
    if year < 2009:
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
        shutil.copy(
                enviro_assoc_input_filename,
                os.path.join(raw_folder, InputDataFilenames.EIA_860_ENVIRO_ASSOC_FILE_FMT.format(year)))
        shutil.copy(
                enviro_assoc_input_filename,
                os.path.join(raw_folder, InputDataFilenames.EIA_860_ENVIRO_EQUIP_FILE_FMT.format(year)))


def format_raw_eia923(year: int):
    """Makes sure that a folder of raw EIA-923 has files with consistent names."""
    if year < 2008:
        raise NotImplementedError(f'We haven\'t implemented data cleaning for {year} yet.')

    raw_folder = f"{downloads_folder()}eia923/f923_{year}"
    consistent_filename = os.path.join(raw_folder, InputDataFilenames.EIA_923_ENVIRONMENTAL_INFO_FMT.format(year))

    # Lots of annoying filename changes for each year.
    if year == 2008:
        base_filename = "SCHEDULE 3A 5A 8A 8B 8C 8D 8E 8F 2008.xlsm"
    elif year == 2009:
        base_filename = "SCHEDULE 3A 5A 8A 8B 8C 8D 8E 8F REVISED 2009 04112011.xls"
    elif year == 2010:
        base_filename = 'SCHEDULE 3A 5A 8A 8B 8C 8D 8E 8F 2010 on NOV 30 2011.xls'
    elif year == 2011:
        base_filename = 'EIA923_Schedule_8_PartsA-F_EnvData_' + str(year) + '{}.xlsx'
    elif year == 2013:
        base_filename = 'EIA923_Schedule_8_PartsA-D_EnvData_' + str(year) + '{}.xlsx'
    elif year == 2017:
        base_filename = 'EIA923_Schedule_8_Annual_Envir_Infor_' + str(year) + '{}.xlsx'
    else:
        base_filename = 'EIA923_Schedule_8_Annual_Environmental_Information_' + str(year) + '{}.xlsx'

    # Sometimes the filename ends in 'Final' instead of 'Final_Revision'.
    final_revision_filename = os.path.join(raw_folder, base_filename.format('_Final_Revision'))
    final_filename = os.path.join(raw_folder, base_filename.format('_Final'))

    if os.path.exists(final_revision_filename):
        # Only copy if the filename is different.
        if consistent_filename != final_revision_filename:
            shutil.copy(final_revision_filename, consistent_filename)
    elif os.path.exists(final_filename):
        # Only copy if the filename is different.
        if consistent_filename != final_filename:
            shutil.copy(final_filename, consistent_filename)
    else:
        raise NotImplementedError()


def check_required_files_raw_eia860(year: int):
    """
    Checks that all required files are present in an EIA-860 data folder.
    """
    raw_folder = f"{downloads_folder()}eia860/eia860{year}"
    if not os.path.exists(raw_folder):
        raise FileNotFoundError(f'The supplied EIA-860 data folder \'{raw_folder}\' does not exist.')

    enviro_assoc_file = InputDataFilenames.EIA_860_ENVIRO_ASSOC_FILE_FMT.format(year)
    enviro_equip_file = InputDataFilenames.EIA_860_ENVIRO_EQUIP_FILE_FMT.format(year)
    if not os.path.exists(os.path.join(raw_folder, enviro_assoc_file)):
        raise FileNotFoundError(f'EIA-860 for year {year} is missing {enviro_assoc_file}')
    if not os.path.exists(os.path.join(raw_folder, enviro_equip_file)):
        raise FileNotFoundError(f'EIA-860 for year {year} is missing {enviro_equip_file}')


def check_required_files_raw_eia923(year: int):
    """
    Checks that all required files are present in an EIA-923 data folder.
    """
    raw_folder = f"{downloads_folder()}eia923/f923_{year}"
    if not os.path.exists(raw_folder):
        raise FileNotFoundError(f'The supplied EIA-923 data folder \'{raw_folder}\' does not exist.')

    environmental_info_file = InputDataFilenames.EIA_923_ENVIRONMENTAL_INFO_FMT.format(year)
    if not os.path.exists(os.path.join(raw_folder, environmental_info_file)):
        raise FileNotFoundError(f'EIA-923 for year {year} is missing {environmental_info_file}')
