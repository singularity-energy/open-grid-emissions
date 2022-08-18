import sys
sys.path.append('../')

import src.load_data as load_data


def load_emissions_controls_helper(years):
    for year in years:
        print(f'-- Loading emission controls data from EIA-923 for {year}')
        df_emission_controls = load_data.load_emissions_controls_eia923(year)
        print('Columns:\n', df_emission_controls.columns)
        assert(len(df_emission_controls) > 0)
        print('Looks good!')


def load_boiler_nox_association_eia860_helper(years, expect_empty=False):
    for year in years:
        print(f'-- Loading boiler NOx information from EIA-923 for {year}')
        df_boiler_nox = load_data.load_boiler_nox_association_eia860(year)
        print('Columns:\n', df_boiler_nox.columns)
        if expect_empty:
            assert(len(df_boiler_nox) == 0)
        else:
            assert(len(df_boiler_nox) > 0)
        print('Looks good!')


def load_boiler_so2_association_eia860_helper(years, expect_empty=False):
    for year in years:
        print(f'-- Loading boiler SO2 information from EIA-923 for {year}')
        df_boiler_nox = load_data.load_boiler_so2_association_eia860(year)
        print('Columns:\n', df_boiler_nox.columns)
        if expect_empty:
            assert(len(df_boiler_nox) == 0)
        else:
            assert(len(df_boiler_nox) > 0)
        print('Looks good!')


def test_load_emissions_controls_eia923_post_2015():
    load_emissions_controls_helper(list(reversed(range(2016, 2021))))


def test_load_emissions_controls_eia923_2012_to_2015():
    load_emissions_controls_helper(list(reversed(range(2012, 2016))))


def test_load_emissions_controls_eia923_pre_2012():
    """Not available before 2008 because EIA publishes 906/920 instead of 923."""
    load_emissions_controls_helper(list(reversed(range(2008, 2012))))


def test_load_boiler_nox_association_eia860():
    load_boiler_nox_association_eia860_helper(list(reversed(range(2013, 2021))))
    load_boiler_nox_association_eia860_helper(list(reversed(range(2005, 2013))), expect_empty=True)


def test_load_boiler_so2_association_eia860():
    load_boiler_so2_association_eia860_helper(list(reversed(range(2013, 2021))))
    load_boiler_so2_association_eia860_helper(list(reversed(range(2005, 2013))), expect_empty=True)
