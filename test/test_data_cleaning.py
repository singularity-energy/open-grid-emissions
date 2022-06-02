import sys
sys.path.append('../')

import pytest
import pandas as pd

import src.data_cleaning as data_cleaning


@pytest.fixture
def df_co2_ch4_no2():
    df = pd.DataFrame({
        'co2_mass_lb': [1.0, 2.0],
        'ch4_mass_lb': [3.0, 6.0],
        'n2o_mass_lb': [0.5, 0.2],
        'co2_mass_lb_adjusted': [1.0, 2.0],
        'ch4_mass_lb_adjusted': [3.0, 6.0],
        'n2o_mass_lb_adjusted': [0.5, 0.2],
        'co2_mass_lb_for_electricity': [1.0, 2.0],
        'ch4_mass_lb_for_electricity': [3.0, 6.0],
        'n2o_mass_lb_for_electricity': [0.5, 0.2]
    })
    return df


def test_co2_eq_AR5_20(df_co2_ch4_no2):
    """Compute CO2-eq using the AR5 20-year GWP factors."""
    df = data_cleaning.calculate_co2_eq_mass(
        df_co2_ch4_no2, ipcc_version='AR5', gwp_horizon=20, ar5_climate_carbon_feedback=False)

    assert('co2_eq_mass_lb' in df.columns and \
           'co2_eq_mass_lb_adjusted' in df.columns and \
           'co2_eq_mass_lb_for_electricity' in df.columns)

    # Check some values.
    assert(df['co2_eq_mass_lb'].iloc[0] == (1.0 + 84.0*3.0 + 264*0.5))
    assert(df['co2_eq_mass_lb_for_electricity'].iloc[0] == (1.0 + 84.0*3.0 + 264*0.5))


def test_co2_eq_AR5f_100(df_co2_ch4_no2):
    """Compute CO2-eq using the AR5 100-year GWP factors, and include CCF."""
    df = data_cleaning.calculate_co2_eq_mass(
        df_co2_ch4_no2, ipcc_version='AR5', gwp_horizon=100, ar5_climate_carbon_feedback=True)
    assert('co2_eq_mass_lb' in df.columns and \
           'co2_eq_mass_lb_adjusted' in df.columns and \
           'co2_eq_mass_lb_for_electricity' in df.columns)

    # Check some values.
    assert(df['co2_eq_mass_lb'].iloc[0] == (1.0 + 34.0*3.0 + 298*0.5))
    assert(df['co2_eq_mass_lb_for_electricity'].iloc[0] == (1.0 + 34.0*3.0 + 298*0.5))


def test_co2_eq_SAR_100(df_co2_ch4_no2):
    """Compute CO2-eq using the SAR 100-year GWP factors."""
    df = data_cleaning.calculate_co2_eq_mass(
        df_co2_ch4_no2, ipcc_version='SAR', gwp_horizon=100)

    assert('co2_eq_mass_lb' in df.columns and \
           'co2_eq_mass_lb_adjusted' in df.columns and \
           'co2_eq_mass_lb_for_electricity' in df.columns)

    # Check some values.
    assert(df['co2_eq_mass_lb'].iloc[0] == (1.0 + 21.0*3.0 + 310*0.5))
    assert(df['co2_eq_mass_lb_for_electricity'].iloc[0] == (1.0 + 21.0*3.0 + 310*0.5))


def test_co2_eq_exceptions(df_co2_ch4_no2):
    """Make sure that the input arguments are validated."""
    # Can't use CCF with any ipcc_version other than `AR5`.
    with pytest.raises(ValueError):
        df = data_cleaning.calculate_co2_eq_mass(
            df_co2_ch4_no2, ipcc_version='TAR', gwp_horizon=100, ar5_climate_carbon_feedback=True)

    # Only 20-year and 100-year factors.
    with pytest.raises(ValueError):
        df = data_cleaning.calculate_co2_eq_mass(
            df_co2_ch4_no2, ipcc_version='TAR', gwp_horizon=150)

    # Use an invalid IPCC version.
    with pytest.raises(ValueError):
        df = data_cleaning.calculate_co2_eq_mass(
            df_co2_ch4_no2, ipcc_version='SOME_OTHER_NAME')
    
    # Require CH4 and N2O columns.
    with pytest.raises(ValueError):
        df = df_co2_ch4_no2[['co2_mass_lb', 'co2_mass_lb_adjusted', 'co2_mass_lb_for_electricity']]
        df = data_cleaning.calculate_co2_eq_mass(df, ipcc_version='AR5')
