from typing import Optional
import gzip
import os
import requests
import shutil
import tarfile
import zipfile

from filepaths import downloads_folder, data_folder


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
        print(f"    {final_destination.split('/')[-1]} already downloaded, skipping.")
        return False

    # Otherwise, download to the file in chunks.
    print(f"    Downloading {final_destination.split('/')[-1]}")
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
            existing_version = f.readlines()[0].replace("\n", "")
        if pudl_version == existing_version:
            print("    PUDL version already downloaded")
            return
        else:
            print("    Downloading new version of pudl")
            shutil.rmtree(downloads_folder("pudl"))

    download_pudl(zenodo_url, pudl_version)


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
    print("    Downloading PUDL. Progress: 100.0%")

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


def download_egrid_file(year):
    """
    Downloads the egrid excel files.
    """
    os.makedirs(downloads_folder("egrid"), exist_ok=True)

    egrid_urls = {
        2018: "https://www.epa.gov/sites/default/files/2020-03/egrid2018_data_v2.xlsx",
        2019: "https://www.epa.gov/sites/default/files/2021-02/egrid2019_data.xlsx",
        2020: "https://www.epa.gov/system/files/documents/2022-01/egrid2020_data.xlsx",
        2021: "https://www.epa.gov/system/files/documents/2023-01/eGRID2021_data.xlsx",
    }

    url = egrid_urls[year]
    filename = url.split("/")[-1]
    filepath = downloads_folder(f"egrid/{filename}")
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
        raise NotImplementedError(f"WARNING: We haven't tested EIA-860 for '{year}'.")
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
