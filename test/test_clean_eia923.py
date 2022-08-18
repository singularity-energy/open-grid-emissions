import sys
sys.path.append('../')

import src.data_cleaning as data_cleaning


def test_clean_eia923():
    for year in list(reversed(range(2009, 2020))):
        print(f'--- Testing EIA-923 cleaning for {year}')
        data_cleaning.clean_eia923(year, False, add_subplant_id=False)


def test_clean_eia923_2008():
    for year in list(reversed(range(2005, 2008))):
        data_cleaning.clean_eia923(year, False, add_subplant_id=False)

    # data_cleaning.clean_eia923(2008, False, add_subplant_id=False)
