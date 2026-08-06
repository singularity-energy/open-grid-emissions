import pandas as pd
import pytest

import oge.helpers as helpers
import oge.impute_hourly_profiles as impute_hourly_profiles


def test_create_local_year_timestamps_non_leap_year():
    df = helpers.create_local_year_timestamps(2025, "America/New_York")
    assert len(df) == 8760
    assert df["datetime_utc"].nunique() == 8760
    assert df["datetime_utc"].min() == pd.Timestamp("2025-01-01 05:00:00", tz="UTC")
    assert df["datetime_utc"].max() == pd.Timestamp("2026-01-01 04:00:00", tz="UTC")


def test_create_local_year_timestamps_leap_year():
    df = helpers.create_local_year_timestamps(2024, "America/New_York")
    assert len(df) == 8784


@pytest.fixture
def native_hourly_profiles():
    """A single BA's (SOCO, US/Central) absolute hourly profile for one fuel category,
    for a full non-leap year, with a constant profile value for simplicity."""
    year = 2025
    central_grid = helpers.create_local_year_timestamps(year, "US/Central")
    central_grid["ba_code"] = "SOCO"
    central_grid["fuel_category"] = "natural_gas"
    central_grid["profile"] = 1.0
    central_grid["flat_profile"] = 1.0
    central_grid["profile_method"] = "residual_profile"
    return central_grid[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "report_date",
            "profile",
            "flat_profile",
            "profile_method",
        ]
    ]


@pytest.fixture
def plant_attributes_mixed_timezones():
    """Plants in SOCO, one matching the BA's own timezone and one off-timezone."""
    return pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "ba_code": ["SOCO", "SOCO"],
            "timezone": ["US/Central", "America/New_York"],
        }
    )


def test_expand_hourly_profiles_to_plant_timezones(
    native_hourly_profiles, plant_attributes_mixed_timezones
):
    year = 2025
    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        native_hourly_profiles, plant_attributes_mixed_timezones, year
    )

    # only the off-native timezone should produce rows; the ba's own timezone is
    # already covered by the native profile and doesn't need reassembling
    assert set(expanded["timezone"].unique()) == {"America/New_York"}

    eastern_grid = helpers.create_local_year_timestamps(year, "America/New_York")
    assert len(expanded) == len(eastern_grid) == 8760
    assert expanded["datetime_utc"].nunique() == 8760
    assert set(expanded["datetime_utc"]) == set(eastern_grid["datetime_utc"])

    # every hour of the reassembled year (including the true year-edge hours that
    # required ffill/bfill) should have a real, non-null profile value
    assert (
        expanded[["profile", "flat_profile", "profile_method"]].isna().sum().sum() == 0
    )
    assert (expanded["ba_code"] == "SOCO").all()
    assert (expanded["fuel_category"] == "natural_gas").all()


def test_expand_hourly_profiles_to_plant_timezones_no_off_native_plants(
    native_hourly_profiles,
):
    """If every plant in a ba shares the ba's own timezone, no expansion is needed."""
    year = 2025
    plant_attributes = pd.DataFrame(
        {"plant_id_eia": [1], "ba_code": ["SOCO"], "timezone": ["US/Central"]}
    )

    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        native_hourly_profiles, plant_attributes, year
    )

    assert len(expanded) == 0


def test_expand_hourly_profiles_to_plant_timezones_warns_on_mid_year_nan(
    native_hourly_profiles, plant_attributes_mixed_timezones, caplog
):
    """A NaN in the ba's own native profile mid-year should be logged, not silently
    absorbed into the same fill as the true year-edge hours."""
    year = 2025
    mid_year_mask = native_hourly_profiles["report_date"] == "2025-07-01"
    mid_year_index = native_hourly_profiles.loc[mid_year_mask].index[0]
    native_hourly_profiles.loc[mid_year_index, "profile"] = None

    with caplog.at_level("WARNING"):
        expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
            native_hourly_profiles, plant_attributes_mixed_timezones, year
        )

    assert "missing profile values" in caplog.text
    assert "SOCO" in caplog.text
    # the gap should still get filled, just with a warning surfaced
    assert expanded["profile"].isna().sum() == 0


def test_expand_hourly_profiles_to_plant_timezones_warns_on_missing_hour(
    native_hourly_profiles, plant_attributes_mixed_timezones, caplog
):
    """A missing row (not just a missing value) in the ba's own native profile mid-year
    should also be logged, since the merge can't otherwise tell it apart from a true
    year-edge hour."""
    year = 2025
    mid_year_mask = native_hourly_profiles["report_date"] == "2025-07-01"
    mid_year_index = native_hourly_profiles.loc[mid_year_mask].index[0]
    native_hourly_profiles = native_hourly_profiles.drop(index=mid_year_index)

    with caplog.at_level("WARNING"):
        expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
            native_hourly_profiles, plant_attributes_mixed_timezones, year
        )

    assert "missing hours" in caplog.text
    assert "SOCO" in caplog.text
    assert expanded["profile"].isna().sum() == 0
