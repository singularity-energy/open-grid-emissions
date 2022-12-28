import sys

import pytest
import pandas as pd


@pytest.fixture
def emissions():
    """Need to provide this import as a fixture to avoid complaints from the linter."""
    sys.path.append("../")
    import src.emissions as emissions

    return emissions


@pytest.fixture
def df_co2_ch4_no2():
    df = pd.DataFrame(
        {
            "co2_mass_lb": [1.0, 2.0],
            "ch4_mass_lb": [3.0, 6.0],
            "n2o_mass_lb": [0.5, 0.2],
            "co2_mass_lb_adjusted": [1.0, 2.0],
            "ch4_mass_lb_adjusted": [3.0, 6.0],
            "n2o_mass_lb_adjusted": [0.5, 0.2],
            "co2_mass_lb_for_electricity": [1.0, 2.0],
            "ch4_mass_lb_for_electricity": [3.0, 6.0],
            "n2o_mass_lb_for_electricity": [0.5, 0.2],
        }
    )
    return df


# NOTE: Why are we testing these function specifically?
def test_co2e_AR6_100(df_co2_ch4_no2, emissions):
    """Compute CO2-eq using the AR6 20-year GWP factors."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2021, gwp_horizon=100, ar5_climate_carbon_feedback=False
    )

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_AR6_20(df_co2_ch4_no2, emissions):
    """Compute CO2-eq using the AR6 20-year GWP factors."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2021, gwp_horizon=20, ar5_climate_carbon_feedback=False
    )

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_AR5_20(df_co2_ch4_no2, emissions):
    """Compute CO2-eq using the AR5 20-year GWP factors."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2020, gwp_horizon=20, ar5_climate_carbon_feedback=False
    )

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_AR5f_100(df_co2_ch4_no2, emissions):
    """Compute CO2-eq using the AR5 100-year GWP factors, and include CCF."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2020, gwp_horizon=100, ar5_climate_carbon_feedback=True
    )
    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_SAR_100(df_co2_ch4_no2, emissions):
    """Compute CO2-eq using the SAR 100-year GWP factors."""
    df = emissions.calculate_co2e_mass(df_co2_ch4_no2, 2014, gwp_horizon=100)

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )
