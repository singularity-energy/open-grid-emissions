import oge.download_data as download_data


YEARS_TO_TEST = list(reversed(range(2005, 2021)))


def test_download_pudl_data():
    """Make sure that PUDL data download works."""
    download_data.download_pudl_data(source="aws")
    print("DONE")


def test_download_egrid():
    """Make sure that eGRID download works."""
    download_data.download_egrid_files()
    print("DONE")


def test_download_chalendar_files():
    """Make sure that Chalendar download works."""
    download_data.download_chalendar_files()
    print("DONE")


def test_download_eia_electric_power_annual():
    """Make sure that we can download electric power annual data from EIA."""
    download_data.download_eia_electric_power_annual()
    print("DONE")


def test_download_eia930():
    """Test that EIA-930 data download works for all years."""
    print("Will test the following years:\n", YEARS_TO_TEST)
    for year in YEARS_TO_TEST:
        print(f"Testing EIA-930 download for {year}")
        download_data.download_raw_eia930(years_to_download=[year])
    print("DONE")


def test_download_epa_psdc():
    """Test that the EPA Power Sector Data crosswalk download works."""
    download_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.3/epa_eia_crosswalk.csv"
    )
    print("DONE")


def test_download_raw_eia860():
    """Test that EIA-860 data download works."""
    print("Will test the following years:\n", YEARS_TO_TEST)
    for year in YEARS_TO_TEST:
        print(f"Testing EIA-860 download for {year}")
        download_data.download_raw_eia860(year)


def test_download_raw_eia923():
    """Test that EIA-923 data download works."""
    testable_years = range(2008, 2021)
    print("Will test the following years:\n", testable_years)
    for year in testable_years:
        print(f"Testing EIA-923 download for {year}")
        download_data.download_raw_eia923(year)


def test_download_raw_eia_906_920():
    testable_years = range(2005, 2008)
    for year in testable_years:
        print(f"Testing EIA-906/920 download for {year}")
        download_data.download_raw_eia_906_920(year)
