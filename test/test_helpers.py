import numpy as np
import pandas as pd
import pytest

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


@pytest.fixture
def yearly_plants_nerc(monkeypatch):
    """EIA-860 NERC regions of six plants over several years, returned in place of
    PUDL's `out_eia__yearly_plants` table.

    - Plant 1 reports SERC in 2023 and TRE in 2024.
    - Plant 2 reports WECC in 2023 and no region in 2024.
    - Plant 3 only reports after 2024 (NPCC in 2025, RFC in 2026).
    - Plant 4 is in Alaska and reports ASCC in 2019, and is unknown ("UNK") since 2020,
      as in EIA-860.
    - Plant 5 is in Puerto Rico and plant 6 is in Ohio, and neither reports a region.
    """
    yearly_plants = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 2, 2, 3, 3, 4, 4],
            "report_date": pd.to_datetime(
                [
                    "2023-01-01",
                    "2024-01-01",
                    "2023-01-01",
                    "2024-01-01",
                    "2025-01-01",
                    "2026-01-01",
                    "2019-01-01",
                    "2024-01-01",
                ]
            ),
            "nerc_region": [
                "SERC",
                "TRE",
                "WECC",
                pd.NA,
                "NPCC",
                "RFC",
                "ASCC",
                "UNK",
            ],
        }
    ).astype({"nerc_region": "string"})
    monkeypatch.setattr(
        helpers.load_data,
        "load_pudl_table",
        lambda table, columns=None: yearly_plants,
    )
    return pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 3, 4, 5, 6],
            "state": ["TX", "CA", "NY", "AK", "PR", "OH"],
        }
    )


def test_add_plant_nerc_region_uses_nearest_year(yearly_plants_nerc):
    """Each plant's NERC region should come from its EIA-860 record for the data year,
    or else its most recent earlier record, or else its earliest later record.

    This covers plants that don't report a region for the data year, such as retired or
    CEMS-only plants, the same way BA codes are assigned.
    """
    result = helpers.add_plant_nerc_region(yearly_plants_nerc, 2024).set_index(
        "plant_id_eia"
    )
    # the data year's region is used where it is reported
    assert result.loc[1, "nerc_region"] == "TRE"
    # otherwise the most recent region before the data year
    assert result.loc[2, "nerc_region"] == "WECC"
    # otherwise the earliest region after the data year
    assert result.loc[3, "nerc_region"] == "NPCC"
    # plants that never report a region are left blank
    assert pd.isna(result.loc[6, "nerc_region"])


@pytest.mark.parametrize(
    "year, alaska_region", [(2019, "ASCC"), (2020, "AK"), (2024, "AK")]
)
def test_add_plant_nerc_region_outside_nerc_regions(
    yearly_plants_nerc, year, alaska_region
):
    """Plants in Alaska, Hawaii and Puerto Rico should be labeled as in eGRID.

    EIA-860 reported Alaska plants in ASCC (and Hawaii plants in HICC) until 2019, and
    as unknown since 2020, when eGRID switched to state codes. So an Alaska plant keeps
    its EIA-860 region until 2019 and is "AK" from 2020, even though its last known
    EIA-860 region is ASCC. Puerto Rico plants have no EIA-860 region, and are always
    "PR".
    """
    result = helpers.add_plant_nerc_region(yearly_plants_nerc, year).set_index(
        "plant_id_eia"
    )
    assert result.loc[4, "nerc_region"] == alaska_region
    assert result.loc[5, "nerc_region"] == "PR"
