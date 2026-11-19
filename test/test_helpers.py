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
    """Missing Hg should be set to zero where no fuel was consumed, but only for
    groups that report Hg in the hourly data.

    - Hg stays missing where fuel was consumed, since that is missing data.
    - Plant 2 doesn't report Hg, so it stays blank rather than zero.
    - Plant 3 reports Hg in the hourly data, but none of its Hg is left in the summed
      data (as when every month it operated has a gap), so it should still get zeros.
    """
    df = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 2, 3],
            "fuel_consumed_mmbtu": [10, 0, 10, 0, 0],
            "hg_mass_lb": [0.1, np.nan, np.nan, np.nan, np.nan],
            "hg_mass_lb_for_electricity": [0.1, np.nan, np.nan, np.nan, np.nan],
        }
    )
    hg_data = pd.DataFrame(
        {"plant_id_eia": [1, 2, 3], "hg_mass_lb": [0.1, np.nan, 0.2]}
    )

    result = helpers.fill_hg_without_fuel(df, ["plant_id_eia"], hg_data)

    np.testing.assert_allclose(result["hg_mass_lb"], [0.1, 0, np.nan, np.nan, 0])
    np.testing.assert_allclose(
        result["hg_mass_lb_for_electricity"], [0.1, 0, np.nan, np.nan, 0]
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


def test_add_plant_nerc_region(monkeypatch):
    plants = pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 3, 4, 5, 6],
            "state": ["TX", "CA", "NY", "AK", "PR", "OH"],
        }
    )
    yearly_plants = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 2, 2, 3, 3, 4],
            "report_date": pd.to_datetime(
                [
                    "2023-01-01",
                    "2024-01-01",
                    "2023-01-01",
                    "2024-01-01",
                    "2025-01-01",
                    "2026-01-01",
                    "2024-01-01",
                ]
            ),
            "nerc_region": ["SERC", "TRE", "WECC", pd.NA, "NPCC", "RFC", "WECC"],
        }
    ).astype({"nerc_region": "string"})
    monkeypatch.setattr(
        helpers.load_data,
        "load_pudl_table",
        lambda table, columns=None: yearly_plants,
    )
    result = helpers.add_plant_nerc_region(plants, 2024).set_index("plant_id_eia")
    # the current year's region is used where it is reported
    assert result.loc[1, "nerc_region"] == "TRE"
    # otherwise the most recent region before the current year
    assert result.loc[2, "nerc_region"] == "WECC"
    # otherwise the earliest region after the current year
    assert result.loc[3, "nerc_region"] == "NPCC"
    # Alaska and Puerto Rico are their own regions, as in eGRID
    assert result.loc[4, "nerc_region"] == "AK"
    assert result.loc[5, "nerc_region"] == "PR"
    # plants without a reported region are left blank
    assert pd.isna(result.loc[6, "nerc_region"])
