import sys
sys.path.append('../')

import pandas as pd

import src.load_data as load_data


def load_emissions_controls_helper(years):
    for year in years:
        print(f'-- Loading emission controls data from EIA-923 for {year}')
        df_emission_controls = load_data.load_emissions_controls_eia923(year)
        print('Columns:\n', df_emission_controls.columns)
        assert(len(df_emission_controls) > 0)
        print('Looks good!')


def test_load_emissions_controls_eia923_post_2015():
    load_emissions_controls_helper(list(reversed(range(2016, 2021))))


def test_load_emissions_controls_eia923_2012_to_2015():
    load_emissions_controls_helper(list(reversed(range(2012, 2016))))


def test_load_emissions_controls_eia923_pre_2012():
    load_emissions_controls_helper(list(reversed(range(2008, 2012))))
