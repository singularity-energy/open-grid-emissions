# --------------------------------------------------------------------------------------
# Options:
#  - Use `pytest -rP` to show print statements from PASSED tests after they finish
#  - Use `pytest -s` to direct print statements to the console
# 
# Run a specific test case with:
# pytest name_of_this_file.py -rP -k 'name_of_test_function'

import sys
sys.path.append('../')

import src.download_data as download_data


YEARS_TO_TEST = list(reversed(range(2005, 2021)))


def test_download_pudl_data():
    """
    Make sure that PUDL data download works.
    """
    download_data.download_pudl_data(
        zenodo_url="https://zenodo.org/record/6349861/files/pudl-v0.6.0-2022-03-12.tgz")
    print('DONE')


def test_download_egrid():
    """Make sure that eGRID download works for 2018-2020."""
    egrid_files_to_download = [
        "https://www.epa.gov/sites/default/files/2020-03/egrid2018_data_v2.xlsx",
        "https://www.epa.gov/sites/default/files/2021-02/egrid2019_data.xlsx",
        "https://www.epa.gov/system/files/documents/2022-01/egrid2020_data.xlsx",
    ]
    download_data.download_egrid_files(egrid_files_to_download)
    print('DONE')


def test_download_chalendar_files():
    """Make sure that Chalendar download works."""
    download_data.download_chalendar_files()
    print('DONE')


def test_download_eia930_data():
    """Test that EIA-930 data download works for all years."""
    print('Will test the following years:\n', YEARS_TO_TEST)
    for year in YEARS_TO_TEST:
        print(f'Testing EIA-930 download for {year}')
        download_data.download_eia930_data(years_to_download=[year])
    print('DONE')


def test_download_epa_psdc():
    """Test that the EPA Power Sector Data crosswalk download works."""
    download_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.2.1/epa_eia_crosswalk.csv")
    print('DONE')


def test_download_raw_eia860():
    """Test that EIA-860 data download works."""
    print('Will test the following years:\n', YEARS_TO_TEST)
    for year in YEARS_TO_TEST:
        print(f'Testing EIA-860 download for {year}')
        download_data.download_raw_eia860(year)


def test_format_raw_eia860():
    YEARS_WITH_DIFFERENT_ENVIRO_FILENAMES = range(2009, 2013)
    for year in YEARS_WITH_DIFFERENT_ENVIRO_FILENAMES:
        download_data.format_raw_eia860(year)
        download_data.check_required_files_raw_eia860(year)
