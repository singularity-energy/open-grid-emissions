from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

import oge.data_cleaning as data_cleaning
import oge.helpers as helpers


@pytest.fixture
def plant_timezone_lookup():
    """A fake `core_eia__entity_plants` table with a single Eastern-timezone plant."""
    return pd.DataFrame({"plant_id_eia": [99999], "timezone": ["America/New_York"]})


def test_complete_hourly_timeseries_reindexes_wrong_grid(plant_timezone_lookup):
    """A group with the correct row count but the wrong hours (e.g. shaped onto its
    BA's local year instead of its own) should be corrected, not skipped.

    This reproduces the "8761 hours" bug: a plant is shaped onto its BA's (US/Central)
    local year, but the plant's own timezone is America/New_York. Before the fix, the
    repair step was gated on row count, so a group that already had 8760 rows (just the
    wrong 8760) would be skipped and silently pass validation.
    """
    year = 2025
    central_grid = helpers.create_local_year_timestamps(year, "US/Central")
    df = pd.DataFrame(
        {
            "plant_id_eia": 99999,
            "datetime_utc": central_grid["datetime_utc"],
            "net_generation_mwh": 1.0,
        }
    )
    assert len(df) == 8760

    with patch.object(
        data_cleaning.load_data, "load_pudl_table", return_value=plant_timezone_lookup
    ):
        result = data_cleaning.complete_hourly_timeseries(
            df,
            year,
            group_cols=["plant_id_eia"],
            columns_to_fill_with_zero=["net_generation_mwh"],
        )

    eastern_grid = helpers.create_local_year_timestamps(year, "America/New_York")
    assert len(result) == 8760
    assert result["datetime_utc"].nunique() == 8760
    assert set(result["datetime_utc"]) == set(eastern_grid["datetime_utc"])
    assert result["net_generation_mwh"].isna().sum() == 0


def test_complete_hourly_timeseries_already_correct_grid_is_unchanged(
    plant_timezone_lookup,
):
    """A group already on its own correct local-year grid should come back unchanged,
    even though the function no longer skips groups based on row count."""
    year = 2025
    eastern_grid = helpers.create_local_year_timestamps(year, "America/New_York")
    df = pd.DataFrame(
        {
            "plant_id_eia": 99999,
            "datetime_utc": eastern_grid["datetime_utc"],
            "net_generation_mwh": 1.0,
        }
    )

    with patch.object(
        data_cleaning.load_data, "load_pudl_table", return_value=plant_timezone_lookup
    ):
        result = data_cleaning.complete_hourly_timeseries(
            df,
            year,
            group_cols=["plant_id_eia"],
            columns_to_fill_with_zero=["net_generation_mwh"],
        )

    assert len(result) == 8760
    assert set(result["datetime_utc"]) == set(eastern_grid["datetime_utc"])
    assert (result["net_generation_mwh"] == 1.0).all()


def test_filter_to_ba_local_year_drops_union_of_timezones():
    """Reproduces the power-sector "8761 hours" bug: a ba (e.g. SOCO) has real CEMS
    data from plants in two different physical timezones (Central and Eastern), each
    already completed onto its own local year (unchanged, via complete_hourly_
    timeseries). Summing them together to the ba level spans the union of both local
    years -- one extra hour -- unless the result is filtered down to the ba's own
    native local year.
    """
    year = 2025
    central_grid = helpers.create_local_year_timestamps(year, "US/Central")
    eastern_grid = helpers.create_local_year_timestamps(year, "America/New_York")

    # plant 1 (physically Central) and plant 2 (physically Eastern), both assigned
    # to SOCO, each already completed onto its own local year
    fleet_data = pd.concat(
        [
            pd.DataFrame(
                {
                    "ba_code": "SOCO",
                    "datetime_utc": central_grid["datetime_utc"],
                    "net_generation_mwh": 1.0,
                }
            ),
            pd.DataFrame(
                {
                    "ba_code": "SOCO",
                    "datetime_utc": eastern_grid["datetime_utc"],
                    "net_generation_mwh": 1.0,
                }
            ),
        ],
        ignore_index=True,
    )
    fleet_data = (
        fleet_data.groupby(["ba_code", "datetime_utc"], dropna=False)[
            "net_generation_mwh"
        ]
        .sum()
        .reset_index()
    )
    assert len(fleet_data) == 8761  # the union of both local years

    ba_timezone_lookup = pd.DataFrame(
        {"ba_code": ["SOCO"], "timezone_local": ["US/Central"]}
    )
    with patch.object(
        data_cleaning.load_data, "load_ba_reference", return_value=ba_timezone_lookup
    ):
        result = data_cleaning.filter_to_ba_local_year(fleet_data, year)

    assert len(result) == 8760
    assert set(result["datetime_utc"]) == set(central_grid["datetime_utc"])
    # the eastern-only edge hour is dropped, but every hour on soco's own local
    # year keeps its real (summed) value, including the eastern plant's overlap
    assert (result["net_generation_mwh"] > 0).all()


def test_fill_emissions_for_non_emitting_resources():
    df = pd.DataFrame(
        {
            "energy_source_code": ["SUN", "NG", "NG", "NG", "SUN"],
            "fuel_consumed_mmbtu": [100.0, 0.0, np.nan, 50.0, np.nan],
            "net_generation_mwh": [10.0, -1.0, np.nan, 5.0, np.nan],
            "nox_mass_lb": [np.nan, np.nan, np.nan, 2.0, np.nan],
            "so2_mass_lb": [np.nan, np.nan, np.nan, np.nan, np.nan],
        }
    )

    result = data_cleaning.fill_emissions_for_non_emitting_resources(df)

    # clean fuels and resources that did not consume fuel have zero emissions
    assert (result.loc[0, ["nox_mass_lb", "so2_mass_lb"]] == 0).all()
    assert (result.loc[1, ["nox_mass_lb", "so2_mass_lb"]] == 0).all()
    # emissions stay missing when fuel consumption data is missing, even for clean fuels
    assert result.loc[2, ["nox_mass_lb", "so2_mass_lb"]].isna().all()
    assert result.loc[4, ["nox_mass_lb", "so2_mass_lb"]].isna().all()
    # reported values are not changed, and other missing values are not filled
    assert result.loc[3, "nox_mass_lb"] == 2
    assert np.isnan(result.loc[3, "so2_mass_lb"])
    # data columns that are not emissions are not filled
    assert np.isnan(result.loc[2, "net_generation_mwh"])


def test_remove_generators_without_data():
    df = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 2, 2],
            "generator_id": ["A", "A", "B", "B"],
            "report_date": pd.to_datetime(["2024-01-01", "2024-02-01"] * 2),
            "net_generation_mwh": [10.0, np.nan, np.nan, np.nan],
            "fuel_consumed_mmbtu": [np.nan, np.nan, np.nan, np.nan],
            "fuel_consumed_for_electricity_mmbtu": [np.nan, np.nan, np.nan, np.nan],
        }
    )

    result = data_cleaning.remove_generators_without_data(df)

    # generator B does not report any data, so it is removed. Generator A reports data
    # in one month, so both of its months are kept
    assert result["generator_id"].unique().tolist() == ["A"]
    assert len(result) == 2
