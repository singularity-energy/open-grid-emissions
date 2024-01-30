from typing import Optional
import datetime
import gzip
import os
import requests
import shutil
import tarfile
import zipfile

from oge.filepaths import downloads_folder, data_folder
from oge.logging_util import get_logger

logger = get_logger(__name__)


def download_helper(
    input_url: str,
    download_path: str,
    output_path: Optional[str] = None,
    requires_unzip: bool = False,
    requires_untar: bool = False,
    requires_gzip: bool = False,
    should_clean: bool = False,
    chunk_size: int = 1024,
) -> bool:
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
        logger.info(f"{final_destination.split('/')[-1]} already downloaded, skipping.")
        return False

    # Otherwise, download to the file in chunks.
    logger.info(f"Downloading {final_destination.split('/')[-1]}")
    r = requests.get(input_url, stream=True)
    with open(download_path, "wb") as fd:
        for chunk in r.iter_content(chunk_size=chunk_size):
            fd.write(chunk)
    # Optionally unzip the downloaded file.
    if requires_unzip:
        if output_path is None:
            raise ValueError("Unzipping requires an output_path destination.")
        with zipfile.ZipFile(download_path, "r") as zip_to_unzip:
            zip_to_unzip.extractall(output_path)
    # Optionally un-tar the downloaded file.
    elif requires_untar:
        if output_path is None:
            raise ValueError("Extracting a tar requires an output_path destination.")
        with tarfile.open(download_path) as tar:
            tar.extractall(output_path)
    elif requires_gzip:
        if output_path is None:
            raise ValueError("Extracting a gzip requires an output_path destination.")
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


def download_pudl_data(source: str = "aws"):
    """
    Downloads the pudl database. OGE currently supports two sources: zenodo and aws
    (i.e. nightly builds). For more information about data sources see:
    https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#data-access

    Zenodo provides stable, versioned data based on the output of the `main` branch of
    pudl but is updated less freqently.
    The most recent version can be found at:
    https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#zenodo-archives

    As of 12/2/2023, the most recent zenodo data was PUDL Data Release v2022.11.30.

    The `aws` source downloads data from the Catalyst's AWS Open Data Registry. This
    data is updated nightly based on the most recent `dev` branch of pudl so is less
    stable.

    Inputs:
        `source`: where to download pudl from, either "aws" or "zenodo"
    """
    os.makedirs(downloads_folder("pudl"), exist_ok=True)

    if source == "aws":
        # define the urls
        pudl_db_url = "https://s3.us-west-2.amazonaws.com/pudl.catalyst.coop/v2023.12.01/pudl.sqlite.gz"
        epacems_parquet_url = "https://s3.us-west-2.amazonaws.com/pudl.catalyst.coop/v2023.12.01/hourly_emissions_epacems.parquet"

        # download the pudl sqlite database
        if not os.path.exists(downloads_folder("pudl/pudl.sqlite")):
            output_filepath = downloads_folder("pudl/pudl.sqlite")
            download_helper(
                pudl_db_url,
                download_path=output_filepath + ".gz",
                output_path=output_filepath,
                requires_gzip=True,
                should_clean=True,
            )

            # add a version file
            with open(downloads_folder("pudl/pudl_sqlite_version.txt"), "w+") as v:
                v.write(f"{datetime.date.today()}")
        else:
            with open(downloads_folder("pudl/pudl_sqlite_version.txt"), "r") as f:
                existing_version = f.readlines()[0].replace("\n", "")
            logger.info(
                f"Using nightly build version of PUDL sqlite database downloaded {existing_version}"
            )

        if not os.path.exists(
            downloads_folder("pudl/hourly_emissions_epacems.parquet")
        ):
            # download the epacems parquet
            output_filepath = downloads_folder("pudl/hourly_emissions_epacems.parquet")
            download_helper(
                epacems_parquet_url,
                download_path=output_filepath,
            )

            # add a version file
            with open(downloads_folder("pudl/epacems_parquet_version.txt"), "w+") as v:
                v.write(f"{datetime.date.today()}")

        else:
            with open(downloads_folder("pudl/epacems_parquet_version.txt"), "r") as f:
                existing_version = f.readlines()[0].replace("\n", "")
            logger.info(
                f"Using nightly build version of PUDL epacems parquet file downloaded {existing_version}"
            )
    elif source == "zenodo":
        # NOTE: This is the most recent available version as of 12/2/2023
        zenodo_url = "https://zenodo.org/record/7472137/files/pudl-v2022.11.30.tgz"

        # get the version number
        pudl_version = zenodo_url.split("/")[-1].replace(".tgz", "")

        # if the pudl data already exists, do not re-download
        if os.path.exists(downloads_folder("pudl_zenodo")):
            pudl_version_file = downloads_folder("pudl_zenodo/pudl_version.txt")
            with open(pudl_version_file, "r") as f:
                existing_version = f.readlines()[0].replace("\n", "")
            if pudl_version == existing_version:
                logger.info("Most recent PUDL Zenodo archive already downloaded.")
                return
            else:
                logger.info("Downloading new version of pudl")
                shutil.rmtree(downloads_folder("pudl_zenodo"))

        download_pudl_from_zenodo(zenodo_url, pudl_version)
    else:
        raise ValueError(
            f"{source} is an invalid option for `source`. Must be 'aws' \
                         or 'zenodo'."
        )


def download_pudl_from_zenodo(zenodo_url, pudl_version):
    r = requests.get(zenodo_url, params={"download": "1"}, stream=True)
    # specify parameters for progress bar
    total_size_in_bytes = int(r.headers.get("content-length", 0))
    block_size = 1024 * 1024 * 10  # 10 MB
    downloaded = 0
    logger.info("Downloading PUDL data...")
    with open(downloads_folder("pudl.tgz"), "wb") as fd:
        for chunk in r.iter_content(chunk_size=block_size):
            print(
                f"Progress: {(round(downloaded/total_size_in_bytes*100,2))}%   \r",
                end="",
            )
            fd.write(chunk)
            downloaded += block_size
        print("Progress: 100.0%")

    # extract the tgz file
    logger.info("Extracting PUDL data...")
    with tarfile.open(downloads_folder("pudl.tgz")) as tar:
        tar.extractall(data_folder())

    # rename the extracted directory to pudl_zenodo
    os.rename(data_folder(pudl_version), downloads_folder("pudl_zenodo"))

    # add a version file
    with open(downloads_folder("pudl_zenodo/pudl_version.txt"), "w+") as v:
        v.write(pudl_version)

    # delete the downloaded tgz file
    os.remove(downloads_folder("pudl.tgz"))


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
        output_filepath = downloads_folder(f"eia930/chalendar/{output_filename}")
        download_helper(
            url,
            output_filepath + ".gz",
            output_filepath,
            requires_gzip=True,
            should_clean=True,
        )


def download_egrid_files():
    """
    Downloads the egrid excel files.

    Inputs:
        `urls_to_download`: a list of urls for the excel files that you want to download
    """
    os.makedirs(downloads_folder("egrid"), exist_ok=True)

    # the 2019 and 2020 data appear to be hosted on different urls
    egrid_urls = {
        2018: "https://www.epa.gov/sites/default/files/2020-03/egrid2018_data_v2.xlsx",
        2019: "https://www.epa.gov/sites/default/files/2021-02/egrid2019_data.xlsx",
        2020: "https://www.epa.gov/system/files/documents/2022-09/eGRID2020_Data_v2.xlsx",
        2021: "https://www.epa.gov/system/files/documents/2023-01/eGRID2021_data.xlsx",
        2022: "https://www.epa.gov/system/files/documents/2024-01/egrid2022_data.xlsx",
    }

    for year, url in egrid_urls.items():
        filepath = downloads_folder(f"egrid/egrid{year}_data.xlsx")
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
                download_filepath = downloads_folder(
                    f"eia930/EIA930_{description}_{year}_{months}.csv"
                )
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
        raise NotImplementedError(f"EIA-923 data is unavailable for '{year}'.")
    os.makedirs(downloads_folder("eia923"), exist_ok=True)
    url = f"https://www.eia.gov/electricity/data/eia923/xls/f923_{year}.zip"
    archive_url = (
        f"https://www.eia.gov/electricity/data/eia923/archive/xls/f923_{year}.zip"
    )
    filename = url.split("/")[-1].split(".")[0]
    download_filepath = downloads_folder(f"eia923/{filename}.zip")
    output_filepath = downloads_folder(f"eia923/{filename}")
    try:
        download_helper(
            url,
            download_filepath,
            output_filepath,
            requires_unzip=True,
            should_clean=True,
        )
    except Exception:
        download_helper(
            archive_url,
            download_filepath,
            output_filepath,
            requires_unzip=True,
            should_clean=True,
        )


def download_raw_eia_906_920(year):
    """
    For years before 2008, the EIA releases Form 906 and 920 instead of 923.

    Inputs:
        `year`: A four-digit year.
    """
    if year < 2005 or year > 2007:
        raise NotImplementedError(f"EIA-906/920 data is unavailable for '{year}'.")
    output_folder = f"f906920_{year}"
    download_helper(
        f"https://www.eia.gov/electricity/data/eia923/archive/xls/f906920_{year}.zip",
        downloads_folder(os.path.join("eia923", output_folder + ".zip")),
        downloads_folder(os.path.join("eia923", output_folder)),
        requires_unzip=True,
        should_clean=True,
    )


def download_raw_eia860(year):
    """
    Downloads raw EIA-860 data (zip files), and unzips them to the downloads folder.
    """
    if year < 2005:
        raise NotImplementedError(f"We haven't tested EIA-860 for '{year}'.")
    os.makedirs(downloads_folder("eia860"), exist_ok=True)
    url = f"https://www.eia.gov/electricity/data/eia860/xls/eia860{year}.zip"
    archive_url = (
        f"https://www.eia.gov/electricity/data/eia860/archive/xls/eia860{year}.zip"
    )
    filename = url.split("/")[-1].split(".")[0]
    zip_filepath = downloads_folder(f"eia860/{filename}.zip")
    output_filepath = downloads_folder(f"eia860/{filename}")
    try:
        download_helper(
            url, zip_filepath, output_filepath, requires_unzip=True, should_clean=True
        )
    except Exception:
        download_helper(
            archive_url,
            zip_filepath,
            output_filepath,
            requires_unzip=True,
            should_clean=True,
        )


def download_eia_electric_power_annual():
    """
    Downloads EIA Electric Power Annual uncontrolled emission factors.

    See: https://www.eia.gov/electricity/annual/
    """
    os.makedirs(downloads_folder("eia_electric_power_annual"), exist_ok=True)
    urls = [
        "https://www.eia.gov/electricity/annual/xls/epa_a_01.xlsx",
        "https://www.eia.gov/electricity/annual/xls/epa_a_02.xlsx",
        "https://www.eia.gov/electricity/annual/xls/epa_a_03.xlsx",
        "https://www.eia.gov/electricity/annual/xls/epa_a_04.xlsx",
    ]
    output_filenames = [
        "epa_a_01_so2_uncontrolled_efs.xlsx",
        "epa_a_02_nox_uncontrolled_efs.xlsx",
        "epa_a_03_co2_uncontrolled_efs.xlsx",
        "epa_a_04_nox_reduction_factors.xlsx",
    ]
    for url, output_filename in zip(urls, output_filenames):
        output_filepath = downloads_folder(
            f"eia_electric_power_annual/{output_filename}"
        )
        download_helper(url, output_filepath)
