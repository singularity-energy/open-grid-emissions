import numpy as np
import pytest
import pandas as pd

import oge.emissions as emissions


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


def test_co2e_AR6_100(df_co2_ch4_no2):
    """Compute CO2-eq using the AR6 20-year GWP factors."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2021, gwp_horizon=100, ar5_climate_carbon_feedback=False
    )

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_AR6_20(df_co2_ch4_no2):
    """Compute CO2-eq using the AR6 20-year GWP factors."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2021, gwp_horizon=20, ar5_climate_carbon_feedback=False
    )

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_AR5_20(df_co2_ch4_no2):
    """Compute CO2-eq using the AR5 20-year GWP factors."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2020, gwp_horizon=20, ar5_climate_carbon_feedback=False
    )

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_AR5f_100(df_co2_ch4_no2):
    """Compute CO2-eq using the AR5 100-year GWP factors, and include CCF."""
    df = emissions.calculate_co2e_mass(
        df_co2_ch4_no2, 2020, gwp_horizon=100, ar5_climate_carbon_feedback=True
    )
    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_co2e_SAR_100(df_co2_ch4_no2):
    """Compute CO2-eq using the SAR 100-year GWP factors."""
    df = emissions.calculate_co2e_mass(df_co2_ch4_no2, 2014, gwp_horizon=100)

    assert (
        "co2e_mass_lb" in df.columns
        and "co2e_mass_lb_adjusted" in df.columns
        and "co2e_mass_lb_for_electricity" in df.columns
    )


def test_adjust_emissions_for_biomass_cems():
    """Each hourly CEMS row is reduced by its subplant-month's biomass fraction, and
    rows without a fraction fall back to the primary fuel of the unit."""
    subplant_emission_factors = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 3],
            "subplant_id": [1, 1, 2, 1],
            "report_date": pd.to_datetime(
                ["2024-01-01", "2024-02-01", "2024-01-01", "2024-01-01"]
            ),
            "biomass_co2_fraction": [0.25, 0.0, 1.0, np.nan],
        }
    )
    cems = pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 1, 1, 3, 1, 2],
            "subplant_id": [1, 1, 2, 1, 1, 1, 2],
            "report_date": pd.to_datetime(
                [
                    "2024-01-01",
                    "2024-01-01",
                    "2024-01-01",
                    "2024-02-01",
                    "2024-01-01",
                    "2024-01-01",
                    "2024-01-01",
                ]
            ),
            "energy_source_code": ["NG", "WDS", "NG", "WDS", "MSB", "NG", "NG"],
            "co2_mass_lb": [100.0, 100.0, 100.0, 100.0, 100.0, 200.0, 100.0],
        },
        index=[10, 3, 7, 0, 25, 4, 12],
    )

    result = emissions.adjust_emissions_for_biomass(
        cems, subplant_emission_factors=subplant_emission_factors
    )

    expected = pd.Series(
        [
            75.0,  # fraction of 0.25
            0.0,  # no fraction, so the biomass primary fuel is used
            0.0,  # fraction of 1, even though the primary fuel is not biomass
            100.0,  # fraction of 0, even though the primary fuel is biomass
            0.0,  # undefined fraction, so the biomass primary fuel is used
            150.0,  # same subplant-month as the first row, in a different hour
            100.0,  # no fraction, so the fossil primary fuel is used
        ],
        index=[10, 3, 7, 0, 25, 4, 12],
        name="co2_mass_lb_adjusted",
    )
    pd.testing.assert_series_equal(result["co2_mass_lb_adjusted"], expected)


@pytest.fixture
def df_chp():
    """Subplant-month data modeled on the University of Washington Power Plant (54809)
    in 2023 (plant 1, subplant 1), including a startup fuel with no generation in a
    month when the subplant generated and a month with negative generation and some
    thermal output, plus a non-CHP subplant at the same plant
    (plant 1, subplant 2), a non-CHP plant (plant 2), and a plant with missing
    generation and reported thermal output (plant 3)."""
    df = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 1, 2, 3, 1, 2, 1, 1, 1],
            "subplant_id": [1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1],
            "report_date": pd.to_datetime(
                [
                    "2023-01-01",
                    "2023-02-01",
                    "2023-09-01",
                    "2023-11-01",
                    "2023-01-01",
                    "2023-01-01",
                    "2023-03-01",
                    "2023-02-01",
                    "2023-01-01",
                    "2023-11-01",
                    "2023-12-01",
                ]
            ),
            "energy_source_code": ["NG"] * 9 + ["DFO", "NG"],
            "fuel_consumed_mmbtu": [
                163542.0,
                150040.0,
                30768.0,
                56214.0,
                100.0,
                315.0,
                57819.0,
                100.0,
                100.0,
                50.0,
                100.0,
            ],
            "fuel_consumed_for_electricity_mmbtu": [
                163542.0,
                150040.0,
                30768.0,
                788.0,
                100.0,
                291.0,
                57819.0,
                100.0,
                100.0,
                50.0,
                90.0,
            ],
            "net_generation_mwh": [
                0.0,
                -5.0,
                2.0,
                478.0,
                0.0,
                None,
                None,
                None,
                0.0,
                0.0,
                -5.0,
            ],
        }
    )
    return emissions.calculate_electric_allocation_factor(df)


def test_eaf_chp_zero_generation(df_chp):
    """At a CHP subplant, zero or negative generation with no reported thermal output
    allocates nothing to electricity."""
    assert (df_chp.loc[[0, 1], "electric_allocation_factor"] == 0).all()


def test_eaf_chp_negative_generation_with_thermal_output(df_chp):
    """At a CHP subplant, negative generation allocates nothing to electricity even when
    some thermal output is reported, where the formula would otherwise give one."""
    assert df_chp.loc[10, "electric_allocation_factor"] == 0


def test_eaf_chp_positive_generation(df_chp):
    """At a CHP subplant, positive generation uses the eGRID formula."""
    assert df_chp.loc[2, "electric_allocation_factor"] == 1
    assert df_chp.loc[3, "electric_allocation_factor"] == pytest.approx(
        478 * 3.412 / (478 * 3.412 + 0.6 * (56214 - 788)), rel=1e-4
    )


def test_eaf_chp_zero_generation_on_one_fuel(df_chp):
    """At a CHP subplant, a fuel with zero generation in a month when another fuel
    generated electricity keeps all of its emissions for electricity."""
    assert df_chp.loc[9, "electric_allocation_factor"] == 1


def test_eaf_non_chp_zero_generation(df_chp):
    """At a subplant that is not CHP, zero or missing generation allocates everything to
    electricity, even if another subplant at the same plant is CHP."""
    assert (df_chp.loc[[4, 7, 8], "electric_allocation_factor"] == 1).all()


def test_eaf_missing_generation(df_chp):
    """When generation is missing, the reported fraction of fuel for electricity is
    used, including at a CHP subplant that reports no thermal output that month."""
    assert df_chp.loc[5, "electric_allocation_factor"] == pytest.approx(
        291 / 315, rel=1e-4
    )
    assert df_chp.loc[6, "electric_allocation_factor"] == 1
