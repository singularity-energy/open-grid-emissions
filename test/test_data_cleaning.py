from unittest.mock import patch

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
