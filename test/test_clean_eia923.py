import oge.data_cleaning as data_cleaning


def test_clean_eia923():
    for year in list(reversed(range(2005, 2021))):
        print(f"--- Testing EIA-923 cleaning for {year}")
        data_cleaning.clean_eia923(year)


def test_clean_eia923_2015():
    data_cleaning.clean_eia923(2015)
