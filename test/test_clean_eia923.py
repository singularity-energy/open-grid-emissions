import sys
import pytest


@pytest.fixture
def data_cleaning():
    """Need to provide this import as a fixture to avoid complaints from the linter."""
    sys.path.append('../')
    import src.data_cleaning as data_cleaning
    return data_cleaning


def test_clean_eia923(data_cleaning):
    for year in list(reversed(range(2005, 2021))):
        print(f'--- Testing EIA-923 cleaning for {year}')
        data_cleaning.clean_eia923(year, False, add_subplant_id=False)


def test_clean_eia923_2015(data_cleaning):
    data_cleaning.clean_eia923(2015, False, add_subplant_id=False)
