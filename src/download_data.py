import gzip
import os
import requests
import shutil
import tarfile


def download_pudl_data(zenodo_url):
    """
    Downloads the archived PUDL data release. The most recent version can be found at https://catalystcoop-pudl.readthedocs.io/en/latest/data_access.html#zenodo-archives
    Inputs:
        zenodo_url: the url to the .tgz file hosted on zenodo
    """

    # get the version number
    pudl_version = zenodo_url.split("/")[-1].replace(".tgz", "")

    # if the pudl data already exists, do not re-download
    if os.path.exists("../data/downloads/pudl"):
        with open("../data/downloads/pudl/pudl_version.txt", "r") as f:
            existing_version = f.readlines()[0]
        if pudl_version == existing_version:
            print("PUDL data already downloaded")
        else:
            print("Downloading new version of pudl")
            shutil.rmtree("../data/downloads/pudl")
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
    with open("../data/downloads/pudl.tgz", "wb") as fd:
        for chunk in r.iter_content(chunk_size=block_size):
            print(
                f"Downloading PUDL. Progress: {(round(downloaded/total_size_in_bytes*100,2))}%   \r",
                end="",
            )
            fd.write(chunk)
            downloaded += block_size

    # extract the tgz file
    print("Extracting PUDL data...")
    with tarfile.open("../data/downloads/pudl.tgz") as tar:
        tar.extractall("../data/")

    # rename the extracted directory to pudl so that we don't have to update this for future versions
    os.rename(f"../data/{pudl_version}", "../data/downloads/pudl")

    # add a version file
    with open("../data/downloads/pudl/pudl_version.txt", "w+") as v:
        v.write(pudl_version)

    # delete the downloaded tgz file
    os.remove("../data/downloads/pudl.tgz")

    print("PUDL download complete")


def download_updated_pudl_database(download=True):
    """
    Downloaded the updated `pudl.sqlite` file from datasette.

    This is temporary until a new version of the data is published on zenodo
    """
    if download is True:
        print("Downloading updated pudl.sqlite from Datasette...")
        # remove the existing file from zenodo
        os.remove("../data/downloads/pudl/pudl_data/sqlite/pudl.sqlite")

        r = requests.get("https://data.catalyst.coop/pudl.db", stream=True)
        with open("../data/downloads/pudl/pudl_data/sqlite/pudl.sqlite", "wb") as fd:
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
    if not os.path.exists("../data/downloads/eia930"):
        os.mkdir("../data/downloads/eia930")
    # if there is not a directory for chalendar-formatted files, make it
    if not os.path.exists("../data/downloads/eia930/chalendar"):
        os.mkdir("../data/downloads/eia930/chalendar")

    # download the cleaned and raw files
    urls = [
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_elec.csv.gz",
        "https://gridemissions.s3.us-east-2.amazonaws.com/EBA_raw.csv.gz",
    ]
    for url in urls:
        filename = url.split("/")[-1].replace(".gz", "")
        # if the file already exists, do not re-download it
        if os.path.exists(f"../data/downloads/eia930/chalendar/{filename}"):
            print(f"{filename} already downloaded")
        else:
            r = requests.get(url, stream=True)

            with open(f"../data/downloads/eia930/chalendar/{filename}.gz", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)

            # Unzip
            with gzip.open(
                f"../data/downloads/eia930/chalendar/{filename}.gz", "rb"
            ) as f_in:
                with open(
                    f"../data/downloads/eia930/chalendar/{filename}", "wb"
                ) as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f"../data/downloads/eia930/chalendar/{filename}.gz")


def download_egrid_files(egrid_files_to_download):
    """
    Downloads the egrid excel files
    Inputs:
        egrid_files_to_download: a list of urls for the egrid excel files that you want to download
    """
    # if there is not yet a directory for egrid, make it
    if not os.path.exists("../data/downloads/egrid"):
        os.mkdir("../data/downloads/egrid")

    # download the egrid files
    for url in egrid_files_to_download:
        filename = url.split("/")[-1]
        # if the file already exists, do not re-download it
        if os.path.exists(f"../data/downloads/egrid/{filename}"):
            print(f"{filename} already downloaded")
        else:
            r = requests.get(url, stream=True)

            with open(f"../data/downloads/egrid/{filename}", "wb") as fd:
                for chunk in r.iter_content(chunk_size=1024):
                    fd.write(chunk)


def download_eia930_data(years_to_download):
    """
    Downloads the six month csv files from the EIA-930 website
    Inputs:
        years_to_download: list of four-digit year numbers to download from EIA-930
    """
    # if there is not yet a directory for EIA-930, make it
    if not os.path.exists("../data/downloads/eia930"):
        os.mkdir("../data/downloads/eia930")

    # download the egrid files
    for year in years_to_download:
        for description in ["BALANCE", "INTERCHANGE"]:
            for months in ["Jan_Jun", "Jul_Dec"]:
                if os.path.exists(
                    f"../data/downloads/eia930/EIA930_{description}_{year}_{months}.csv"
                ):
                    print(f"{description}_{year}_{months} data already downloaded")
                else:
                    print(f"downloading {description}_{year}_{months} data")
                    r = requests.get(
                        f"https://www.eia.gov/electricity/gridmonitor/sixMonthFiles/EIA930_{description}_{year}_{months}.csv",
                        stream=True,
                    )

                    with open(
                        f"../data/downloads/eia930/EIA930_{description}_{year}_{months}.csv",
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
    if not os.path.exists("../data/downloads/epa"):
        os.mkdir("../data/downloads/epa")

    filename = psdc_url.split("/")[-1]
    # if the file already exists, do not re-download it
    if os.path.exists(f"../data/downloads/epa/{filename}"):
        print(f"{filename} already downloaded")
    else:
        r = requests.get(psdc_url, stream=True)

        with open(f"../data/downloads/epa/{filename}", "wb") as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)

