import numpy as np
import pandas as pd

import oge.helpers as helpers


def test_flag_missing_hg():
    """Only data that consumed fuel but has no Hg is missing Hg."""
    df = pd.DataFrame(
        {
            "fuel_consumed_mmbtu": [10, 10, 0, np.nan],
            "hg_mass_lb": [0.1, np.nan, np.nan, np.nan],
        }
    )
    assert helpers.flag_missing_hg(df).tolist() == [False, True, False, False]


def test_remove_incomplete_hg():
    """Summed Hg is removed where any of the summed data was missing Hg, and the
    `missing_hg` column is dropped."""
    df = pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 3],
            "hg_mass_lb": [0.3, 0.1, np.nan],
            "hg_mass_lb_for_electricity": [0.2, 0.1, np.nan],
            "missing_hg": [0, 1, 0],
        }
    )

    result = helpers.remove_incomplete_hg(df)

    assert "missing_hg" not in result.columns
    np.testing.assert_allclose(result["hg_mass_lb"], [0.3, np.nan, np.nan])
    np.testing.assert_allclose(
        result["hg_mass_lb_for_electricity"], [0.2, np.nan, np.nan]
    )


def test_fill_hg_without_fuel():
    """Missing Hg is zero without fuel consumption, but only for groups that report
    Hg, and stays missing when fuel was consumed."""
    df = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 2],
            "fuel_consumed_mmbtu": [10, 0, 10, 0],
            "hg_mass_lb": [0.1, np.nan, np.nan, np.nan],
            "hg_mass_lb_for_electricity": [0.1, np.nan, np.nan, np.nan],
        }
    )

    result = helpers.fill_hg_without_fuel(df, ["plant_id_eia"])

    np.testing.assert_allclose(result["hg_mass_lb"], [0.1, 0, np.nan, np.nan])
    np.testing.assert_allclose(
        result["hg_mass_lb_for_electricity"], [0.1, 0, np.nan, np.nan]
    )


def test_combine_subplant_data_monthly_hg_requires_every_hour():
    """A subplant-month's Hg is missing if any hour that consumed fuel has no Hg."""
    jan, feb = pd.Timestamp("2024-01-01"), pd.Timestamp("2024-02-01")
    cems = pd.DataFrame(
        {
            "plant_id_eia": [1] * 5,
            "subplant_id": [1] * 5,
            "report_date": [jan, jan, jan, feb, feb],
            "fuel_consumed_mmbtu": [10, 10, 0, 10, 0],
            "hg_mass_lb": [0.1, np.nan, 0, 0.2, 0],
        }
    )
    empty = cems.iloc[0:0]

    result = helpers.combine_subplant_data(
        cems,
        empty,
        empty,
        empty,
        resolution="monthly",
        data_columns=["fuel_consumed_mmbtu", "hg_mass_lb"],
    )

    np.testing.assert_allclose(result["hg_mass_lb"], [np.nan, 0.2])
    assert result["fuel_consumed_mmbtu"].tolist() == [20, 10]
