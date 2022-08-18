# --------------------------------------------------------------------------------------
# Options:
#  - Use `pytest -rP` to show print statements from PASSED tests after they finish
#  - Use `pytest -s` to direct print statements to the console
# 
# Run a specific test case with:
# pytest name_of_this_file.py -rP -k 'name_of_test_function'
#
# NOTE: The required input data must be downloaded first. See data_pipeline.py.

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


def load_boiler_nox_association_eia860_helper(years):
    for year in years:
        print(f'-- Loading boiler NOx information from EIA-923 for {year}')
        df_boiler_nox = load_data.load_emissions_controls_eia923(year)
        print('Columns:\n', df_boiler_nox.columns)
        assert(len(df_boiler_nox) > 0)
        print('Looks good!')


def test_load_emissions_controls_eia923_post_2015():
    load_emissions_controls_helper(list(reversed(range(2016, 2021))))


def test_load_emissions_controls_eia923_2012_to_2015():
    load_emissions_controls_helper(list(reversed(range(2012, 2016))))


def test_load_emissions_controls_eia923_pre_2012():
    load_emissions_controls_helper(list(reversed(range(2008, 2012))))


# def test_load_boiler_nox_association_eia860():
#     load_boiler_nox_association_eia860_helper(list(reversed(range(2008, 2012))))
