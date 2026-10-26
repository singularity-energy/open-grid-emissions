from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

import oge.helpers as helpers
import oge.impute_hourly_profiles as impute_hourly_profiles
import oge.validation as validation


def _offset(timezone, year=2025):
    """The standard-time UTC offset for a timezone name, matching what
    helpers.add_local_timezone_offset() would compute for the same (timezone, year)."""
    return pd.Timestamp(year, 1, 1, tz=timezone).utcoffset()


def test_calculate_national_average_wind_solar_profiles_groups_by_standard_time():
    """Averaging must group bas by local STANDARD time, not naive prevailing local
    time. During a DST month, a DST-observing ba's prevailing offset differs from its
    standard offset while a non-DST ba's (e.g. Arizona) does not, so grouping by
    prevailing time would misalign the two bas; grouping by standard time correctly
    aligns them onto the same clock slot. This is computed once, for the whole year
    and all bas at once -- not per (ba, fuel, month) -- which is also what avoids a
    real production bug: slicing this average by report_date (a ba-relative,
    prevailing-time tag) while grouping by a different (standard-time) offset let
    the same real hour get computed independently from two different report_dates'
    worth of data, producing duplicate datetime_utc values downstream.
    """
    year = 2025
    residual_profiles = pd.DataFrame(
        {
            "ba_code": ["GRIS", "AZPS"],
            "fuel_category": "wind",
            # 1 real utc hour apart, but the same local STANDARD time in both zones
            # (US/Central is UTC-6 standard; America/Phoenix is UTC-7 year-round,
            # since Arizona does not observe DST)
            "datetime_utc": pd.to_datetime(
                ["2025-07-01 07:00:00", "2025-07-01 08:00:00"], utc=True
            ),
            "eia930_profile": [2.0, 4.0],
        }
    )
    ba_timezone_lookup = pd.DataFrame(
        {
            "ba_code": ["GRIS", "AZPS"],
            "timezone_local": ["US/Central", "America/Phoenix"],
        }
    )

    with patch.object(
        impute_hourly_profiles.load_data,
        "load_ba_reference",
        return_value=ba_timezone_lookup,
    ):
        result = impute_hourly_profiles.calculate_national_average_wind_solar_profiles(
            residual_profiles, year
        )

    # both rows land on local standard time 2025-07-01 01:00:00, so they're
    # averaged into a single row rather than kept separate
    matching = result[
        (result["fuel_category"] == "wind")
        & (result["datetime_local_standard"] == pd.Timestamp("2025-07-01 01:00:00"))
    ]
    assert len(matching) == 1
    assert matching["national_average_profile"].iloc[0] == 3.0


def test_get_national_average_profile_for_fuel_month_uses_targets_own_report_date():
    """get_national_average_profile_for_fuel_month must assign report_date from the
    TARGET ba's own reconstructed local calendar, not from a source ba's
    report_date -- since the precomputed national-average table (see
    calculate_national_average_wind_solar_profiles) is no longer sliced by
    report_date at all, source attribution is structurally impossible. This
    reproduces the ba-relative month boundary that caused the original bug: two
    naive-standard-time entries that straddle GRIS's own March/April boundary must
    land in GRIS's own correct, distinct months.
    """
    year = 2025
    # GRIS (US/Central) standard offset is -6h, but by April 1 its prevailing (DST)
    # offset is already -5h (CDT started March 9). naive-standard 2025-03-31 22:00
    # -> utc 2025-04-01 04:00 -> GRIS's own prevailing local 2025-03-31 23:00
    # (still March); naive-standard 2025-03-31 23:00 -> utc 2025-04-01 05:00 ->
    # GRIS's own prevailing local 2025-04-01 00:00 (April)
    national_average_profiles = pd.DataFrame(
        {
            "fuel_category": ["wind", "wind"],
            "datetime_local_standard": pd.to_datetime(
                ["2025-03-31 22:00:00", "2025-03-31 23:00:00"]
            ),
            "national_average_profile": [10.0, 20.0],
        }
    )
    ba_timezone_lookup = pd.DataFrame(
        {"ba_code": ["GRIS"], "timezone_local": ["US/Central"]}
    )

    with patch.object(
        impute_hourly_profiles.load_data,
        "load_ba_reference",
        return_value=ba_timezone_lookup,
    ):
        march_result = (
            impute_hourly_profiles.get_national_average_profile_for_fuel_month(
                national_average_profiles,
                "GRIS",
                "wind",
                pd.Timestamp("2025-03-01"),
                year,
            )
        )
        april_result = (
            impute_hourly_profiles.get_national_average_profile_for_fuel_month(
                national_average_profiles,
                "GRIS",
                "wind",
                pd.Timestamp("2025-04-01"),
                year,
            )
        )

    march_matches = march_result[march_result["imputed_profile"].notna()]
    april_matches = april_result[april_result["imputed_profile"].notna()]
    assert march_matches["imputed_profile"].tolist() == [10.0]
    assert april_matches["imputed_profile"].tolist() == [20.0]
    # the hour belongs to exactly one of the two adjacent months, never both
    combined_matches = pd.concat([march_matches, april_matches])
    assert not combined_matches["datetime_utc"].duplicated().any()


def test_get_national_average_profile_for_fuel_month_no_dst_ambiguity_or_gap():
    """Standard-time arithmetic is a fixed, bijective offset, so unlike naive
    prevailing local time, reconstructing a target ba's own timeline from it never
    produces an ambiguous ("fall back") or nonexistent ("spring forward") hour --
    even across GRIS's own DST fall-back transition (US/Central falls back on
    2025-11-02).
    """
    year = 2025
    hours = pd.date_range(
        "2025-11-01 22:00:00", "2025-11-02 10:00:00", freq="h", tz="UTC"
    )
    residual_profiles = pd.DataFrame(
        {
            "ba_code": "GRIS",
            "fuel_category": "wind",
            "datetime_utc": hours,
            "eia930_profile": 1.0,
        }
    )
    ba_timezone_lookup = pd.DataFrame(
        {"ba_code": ["GRIS"], "timezone_local": ["US/Central"]}
    )

    with patch.object(
        impute_hourly_profiles.load_data,
        "load_ba_reference",
        return_value=ba_timezone_lookup,
    ):
        national_average_profiles = (
            impute_hourly_profiles.calculate_national_average_wind_solar_profiles(
                residual_profiles, year
            )
        )
        result = impute_hourly_profiles.get_national_average_profile_for_fuel_month(
            national_average_profiles,
            "GRIS",
            "wind",
            pd.Timestamp("2025-11-01"),
            year,
        )

    # result spans all of GRIS's own November, not just the hours seeded above;
    # narrow to those before checking for duplicates/gaps
    result = result[result["datetime_utc"].isin(hours)]
    assert len(result) == len(hours)
    assert not result["datetime_utc"].duplicated().any()
    assert set(result["datetime_utc"]) == set(hours)


def test_create_local_year_timestamps_non_leap_year():
    df = helpers.create_local_year_timestamps(2025, "America/New_York")
    assert len(df) == 8760
    assert df["datetime_utc"].nunique() == 8760
    assert df["datetime_utc"].min() == pd.Timestamp("2025-01-01 05:00:00", tz="UTC")
    assert df["datetime_utc"].max() == pd.Timestamp("2026-01-01 04:00:00", tz="UTC")


def test_create_local_year_timestamps_leap_year():
    df = helpers.create_local_year_timestamps(2024, "America/New_York")
    assert len(df) == 8784


def test_add_local_timezone_offset_collapses_aliases():
    """A legacy alias (e.g. ba_reference.csv's "US/Central") and its canonical IANA
    name (e.g. "America/Chicago") must produce the same local_timezone_offset, while
    a genuinely different timezone (e.g. "America/New_York") must not."""
    year = 2025
    df = pd.DataFrame(
        {"timezone": ["US/Central", "America/Chicago", "America/New_York"]}
    )
    result = helpers.add_local_timezone_offset(df, year)

    assert (
        result.loc[0, "local_timezone_offset"] == result.loc[1, "local_timezone_offset"]
    )
    assert (
        result.loc[0, "local_timezone_offset"] != result.loc[2, "local_timezone_offset"]
    )


def test_convert_profile_to_percent_month_with_mixed_profile_methods():
    """Reproduces the FPL bug: a month's hours split across two profile_method
    values must still be normalized together, summing to 1.0 for the month. If
    "profile_method" is included in group_keys, each method's hours get normalized
    separately, and the month's percentages sum to one per distinct method used
    (2.0 here) instead of 1.0 -- silently doubling the shaped monthly total.
    """
    hours = pd.date_range("2025-01-01", periods=744, freq="h", tz="UTC")
    df = pd.DataFrame(
        {
            "ba_code": "FPL",
            "fuel_category": "natural_gas",
            "report_date": pd.Timestamp("2025-01-01"),
            "datetime_utc": hours,
            # most of the month used one method; a handful of hours fell back to another
            "profile_method": ["residual_profile"] * 700
            + ["shifted_residual_profile"] * 44,
            "profile": 1.0,
        }
    )

    correctly_grouped = impute_hourly_profiles.convert_profile_to_percent(
        df.copy(),
        group_keys=["ba_code", "fuel_category"],
        columns_to_convert=["profile"],
    )
    assert correctly_grouped["profile"].sum() == pytest.approx(1.0)

    incorrectly_grouped = impute_hourly_profiles.convert_profile_to_percent(
        df.copy(),
        group_keys=["ba_code", "fuel_category", "profile_method"],
        columns_to_convert=["profile"],
    )
    assert incorrectly_grouped["profile"].sum() == pytest.approx(2.0)


def test_native_plant_join_uses_consistent_local_timezone_offset():
    """Reproduces a real production bug (AECI): a plant on its ba's own native
    timezone, but tagged with a different raw timezone name than ba_reference.csv
    uses for that ba (the canonical "America/Chicago" vs. the legacy alias
    "US/Central"), must still successfully match its ba's native hourly_profiles
    slice in the actual shaping merge. It's not enough for
    expand_hourly_profiles_to_plant_timezones to correctly judge the plant "native"
    in isolation -- native hourly_profiles rows and the plant's own
    local_timezone_offset must match for the join in shape_monthly_eia_data_as_hourly
    to find anything at all.
    """
    year = 2025
    native_grid = helpers.create_local_year_timestamps(year, "US/Central")
    native_grid["ba_code"] = "AECI"
    native_grid["fuel_category"] = "natural_gas"
    native_grid["profile"] = 1.0  # absolute; convert_profile_to_percent normalizes it
    native_grid["flat_profile"] = 1.0
    native_grid["profile_method"] = "residual_profile"
    # tagged with the ba's native alias -- exactly as data_pipeline.py tags
    # hourly_profiles's native rows
    native_grid["local_timezone_offset"] = _offset("US/Central", year)

    hourly_profiles = native_grid[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "report_date",
            "profile",
            "flat_profile",
            "profile_method",
            "local_timezone_offset",
        ]
    ].copy()
    # normalize per (ba_code, fuel_category, local_timezone_offset, report_date),
    # exactly as data_pipeline.py does, so each month's percentages sum to 1.0
    hourly_profiles = impute_hourly_profiles.convert_profile_to_percent(
        hourly_profiles,
        group_keys=["ba_code", "fuel_category", "local_timezone_offset"],
        columns_to_convert=["profile", "flat_profile"],
    )

    monthly_eia_data_to_shape = pd.DataFrame(
        {
            "plant_id_eia": 1127,
            "subplant_id": 1,
            "ba_code": "AECI",
            "fuel_category": "natural_gas",
            "report_date": pd.date_range("2025-01-01", periods=12, freq="MS"),
            "net_generation_mwh": 100.0,
            "fuel_consumed_mmbtu": 50.0,
            # the plant's own canonical PUDL timezone name -- exactly as
            # combine_and_export_hourly_plant_data tags monthly_eia_data_to_shape
            "local_timezone_offset": _offset("America/Chicago", year),
        }
    )

    shaped = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        monthly_eia_data_to_shape, hourly_profiles, year
    )

    assert shaped["datetime_utc"].isna().sum() == 0
    assert shaped["net_generation_mwh"].sum() == pytest.approx(
        monthly_eia_data_to_shape["net_generation_mwh"].sum()
    )


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
    central_grid["local_timezone_offset"] = _offset("US/Central", year)
    return central_grid[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "report_date",
            "profile",
            "flat_profile",
            "profile_method",
            "local_timezone_offset",
        ]
    ]


@pytest.fixture
def plant_attributes_mixed_timezones():
    """Plants in SOCO, one matching the BA's own timezone and one off-timezone.

    The native plant is tagged with the canonical IANA name ("America/Chicago"),
    matching how real plant_attributes timezones are populated from
    core_eia__entity_plants, rather than the legacy alias ("US/Central") that
    ba_reference.csv uses for the ba's own native timezone. Both names refer to the
    same physical timezone, and so must produce the same local_timezone_offset.
    """
    df = pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "ba_code": ["SOCO", "SOCO"],
            "timezone": ["America/Chicago", "America/New_York"],
        }
    )
    return helpers.add_local_timezone_offset(df, 2025)


def test_expand_hourly_profiles_to_plant_timezones(
    native_hourly_profiles, plant_attributes_mixed_timezones
):
    year = 2025
    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        native_hourly_profiles, plant_attributes_mixed_timezones, year
    )

    # only the off-native timezone should produce rows; the ba's own timezone is
    # already covered by the native profile and doesn't need reassembling
    assert set(expanded["local_timezone_offset"].unique()) == {
        _offset("America/New_York", year)
    }

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


def test_expand_hourly_profiles_fills_hours_from_months_with_another_method(
    native_hourly_profiles, plant_attributes_mixed_timezones
):
    """An Eastern plant in SOCO (Central) starts each month an hour before SOCO. If
    SOCO's previous month uses another method, that hour's value isn't on the same
    scale, so it takes the nearest value of the plant's own month instead. If both
    months use the same method, the hour keeps SOCO's value.

    SOCO's January and April are flat (1), February and March are residuals in MW
    (200 and 300).
    """
    year = 2025
    profiles = native_hourly_profiles.copy()
    month = profiles["report_date"].dt.month
    profiles["profile"] = month.map({2: 200.0, 3: 300.0}).fillna(1.0)
    profiles["profile_method"] = month.map(
        {2: "residual_profile", 3: "residual_profile"}
    ).fillna("assumed_flat")

    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        profiles, plant_attributes_mixed_timezones, year
    ).set_index("datetime_utc")

    # the first hour of the plant's month, which is the last hour of SOCO's month
    def first_hour(day):
        return expanded.loc[pd.Timestamp(day, tz="America/New_York").tz_convert("UTC")]

    # flat January to February in MW: the hour takes February's value
    assert first_hour("2025-02-01")["profile"] == 200.0
    assert first_hour("2025-02-01")["profile_method"] == "residual_profile"
    # February to March, both in MW: the hour keeps SOCO's February value
    assert first_hour("2025-03-01")["profile"] == 200.0
    # March in MW to flat April: the hour takes April's flat value
    assert first_hour("2025-04-01")["profile"] == 1.0
    assert first_hour("2025-04-01")["profile_method"] == "assumed_flat"
    assert expanded["profile"].notna().all()


def test_expand_hourly_profiles_to_plant_timezones_carries_storage_profiles(
    native_hourly_profiles, plant_attributes_mixed_timezones
):
    """Storage profiles should be reassembled onto plant timezones like the other
    profiles, so that storage plants in a different timezone than their BA can be
    shaped with them.
    """
    year = 2025
    profiles = native_hourly_profiles.assign(
        storage_discharge_profile=2.0,
        storage_charge_profile=3.0,
        storage_charge_profile_method="net_load_proxy",
    )
    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        profiles, plant_attributes_mixed_timezones, year
    )
    # every reassembled hour, including the filled year-edge hours, has the storage
    # profiles
    assert (expanded["storage_discharge_profile"] == 2.0).all()
    assert (expanded["storage_charge_profile"] == 3.0).all()
    assert (expanded["storage_charge_profile_method"] == "net_load_proxy").all()


def test_remove_storage_from_eia930_data():
    """The EIA-930 storage categories should never be used as residual profiles, and
    "other" should be removed in the BA-months when it is used for battery storage.
    """
    report_dates = pd.to_datetime(["2025-01-01", "2025-02-01"])
    eia930_data = pd.DataFrame(
        {
            "ba_code": "CISO",
            "fuel_category_eia930": ["other", "other", "battery_storage", "solar"],
            "report_date": report_dates[[0, 1, 1, 1]],
        }
    )
    other_storage_months = pd.DataFrame(
        {"ba_code": ["CISO"], "report_date": report_dates[:1]}
    )

    result = impute_hourly_profiles.remove_storage_from_eia930_data(
        eia930_data, other_storage_months
    )

    # January's "other" and the battery category are removed
    assert result[["fuel_category_eia930", "report_date"]].values.tolist() == [
        ["other", report_dates[1]],
        ["solar", report_dates[1]],
    ]


def test_expand_hourly_profiles_to_plant_timezones_no_off_native_plants(
    native_hourly_profiles,
):
    """If every plant in a ba shares the ba's own timezone, no expansion is needed."""
    year = 2025
    plant_attributes = helpers.add_local_timezone_offset(
        pd.DataFrame(
            {"plant_id_eia": [1], "ba_code": ["SOCO"], "timezone": ["America/Chicago"]}
        ),
        year,
    )

    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        native_hourly_profiles, plant_attributes, year
    )

    assert len(expanded) == 0


def test_expand_hourly_profiles_to_plant_timezones_recognizes_timezone_aliases(
    native_hourly_profiles,
):
    """A plant's timezone should be treated as native even when it doesn't match the
    ba's own timezone string exactly, as long as they refer to the same physical
    timezone. ba_reference.csv uses legacy aliases (e.g. "US/Central") for a ba's
    native timezone, while plant timezones use canonical IANA names (e.g.
    "America/Chicago"); these should not be treated as different timezones just
    because the strings differ.
    """
    year = 2025
    plant_attributes = helpers.add_local_timezone_offset(
        pd.DataFrame(
            {
                "plant_id_eia": [1, 2, 3],
                "ba_code": ["SOCO", "SOCO", "SOCO"],
                # all three of these strings refer to the same physical timezone as
                # ba_reference.csv's "US/Central" entry for SOCO
                "timezone": ["America/Chicago", "US/Central", "Etc/GMT+6"],
            }
        ),
        year,
    )

    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        native_hourly_profiles, plant_attributes, year
    )

    assert len(expanded) == 0


def test_expand_hourly_profiles_to_plant_timezones_warns_only_once_per_ba(
    native_hourly_profiles, caplog
):
    """A ba with more than one non-native timezone should still only be warned about
    once, since the gap being warned about is a property of the ba's own native
    profile, not of any particular target timezone."""
    year = 2025
    mid_year_mask = native_hourly_profiles["report_date"] == "2025-07-01"
    mid_year_index = native_hourly_profiles.loc[mid_year_mask].index[0]
    native_hourly_profiles.loc[mid_year_index, "profile"] = None

    plant_attributes = helpers.add_local_timezone_offset(
        pd.DataFrame(
            {
                "plant_id_eia": [1, 2],
                "ba_code": ["SOCO", "SOCO"],
                # two distinct non-native timezones for the same ba
                "timezone": ["America/New_York", "America/Denver"],
            }
        ),
        year,
    )

    with caplog.at_level("WARNING"):
        expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
            native_hourly_profiles, plant_attributes, year
        )

    matching_warnings = [
        r for r in caplog.records if "missing profile values" in r.message
    ]
    assert len(matching_warnings) == 1
    assert set(expanded["local_timezone_offset"].unique()) == {
        _offset("America/New_York", year),
        _offset("America/Denver", year),
    }
    # the mid-year gap is not a true year-edge run, so it stays missing in each of
    # the two reassembled timezones rather than being filled
    assert expanded["profile"].isna().sum() == 2


def test_expand_hourly_profiles_to_plant_timezones_warns_on_mid_year_nan(
    native_hourly_profiles, plant_attributes_mixed_timezones, caplog
):
    """A NaN in the ba's own native profile mid-year should be logged, and should
    stay missing in the reassembled profile rather than being silently smoothed
    into the same fill as the true year-edge hours."""
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
    # the mid-year gap is not a true year-edge run, so it should remain missing
    assert expanded["profile"].isna().sum() == 1


def test_expand_hourly_profiles_to_plant_timezones_warns_on_missing_hour(
    native_hourly_profiles, plant_attributes_mixed_timezones, caplog
):
    """A missing row (not just a missing value) in the ba's own native profile mid-year
    should also be logged, since the merge can't otherwise tell it apart from a true
    year-edge hour, and should also stay missing rather than being filled."""
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
    assert expanded["profile"].isna().sum() == 1


def test_expand_hourly_profiles_to_plant_timezones_only_fills_true_edges(
    native_hourly_profiles, plant_attributes_mixed_timezones
):
    """A mid-year gap in the ba's own native profile should remain missing in the
    reassembled output, while the genuine start/end-of-year edge hours -- which are
    also missing prior to filling, for the same structural reason -- are still
    filled. This distinguishes an unrelated data gap from the expected edge gap
    that reassembling onto a different timezone always introduces.
    """
    year = 2025
    july_index = native_hourly_profiles.loc[
        native_hourly_profiles["report_date"] == "2025-07-01"
    ].index
    mid_year_index = july_index[300]  # an hour in the middle of july, not an edge
    native_hourly_profiles.loc[mid_year_index, "profile"] = None
    gap_hour = native_hourly_profiles.loc[mid_year_index, "datetime_utc"]

    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        native_hourly_profiles, plant_attributes_mixed_timezones, year
    )

    eastern = expanded[
        expanded["local_timezone_offset"] == _offset("America/New_York", year)
    ].sort_values("datetime_utc")
    # the mid-year gap carries through as missing, not filled
    assert eastern.loc[eastern["datetime_utc"] == gap_hour, "profile"].isna().all()
    # but the true year-edge hours (first/last rows, once reassembled onto a
    # different timezone) are still filled
    assert eastern["profile"].iloc[0] is not None and not pd.isna(
        eastern["profile"].iloc[0]
    )
    assert not pd.isna(eastern["profile"].iloc[-1])
    assert eastern["profile"].isna().sum() == 1


def test_expand_hourly_profiles_to_plant_timezones_caps_fill_at_edge_window(
    native_hourly_profiles, plant_attributes_mixed_timezones
):
    """A leading/trailing gap wider than EDGE_WINDOW_HOURS should only be filled up
    to that window, not filled in its entirety. This bounds how much a single
    unexpected gap (e.g. a larger-than-expected timezone offset, or a real gap in
    the ba's native data that happens to sit at the edge) can be silently smoothed
    over, regardless of how wide it is.
    """
    year = 2025
    # remove the first 6 hours of SOCO's own native year -- combined with the 1 hour
    # America/New_York's reassembled year always needs to borrow from the prior
    # year (which doesn't exist), this produces a 7-hour missing run at the start
    # of the reassembled Eastern year, wider than the 4-hour edge window
    sorted_native = native_hourly_profiles.sort_values("datetime_utc")
    rows_to_drop = sorted_native.index[:6]
    native_hourly_profiles = native_hourly_profiles.drop(index=rows_to_drop)

    expanded = impute_hourly_profiles.expand_hourly_profiles_to_plant_timezones(
        native_hourly_profiles, plant_attributes_mixed_timezones, year
    )

    eastern = expanded[
        expanded["local_timezone_offset"] == _offset("America/New_York", year)
    ].sort_values("datetime_utc")
    # only the first 4 hours of the 7-hour gap should be filled
    assert eastern["profile"].iloc[:4].isna().sum() == 0
    # the remaining 3 hours of the gap should still be missing, not filled
    assert eastern["profile"].iloc[4:7].isna().sum() == 3
    assert eastern["profile"].iloc[7:].isna().sum() == 0


@pytest.fixture
def zero_january_profiles():
    """Hourly profiles for one 2025 fleet, with a profile of zero in every hour of
    January, already converted to a percent of the monthly total.
    """
    year = 2025
    grid = helpers.create_local_year_timestamps(year, "US/Central")
    grid["ba_code"] = "AECI"
    grid["fuel_category"] = "natural_gas"
    grid["profile"] = (grid["report_date"].dt.month != 1).astype(float)
    grid["flat_profile"] = 1.0
    grid["profile_method"] = "residual_profile"
    grid["local_timezone_offset"] = _offset("US/Central", year)
    return impute_hourly_profiles.convert_profile_to_percent(
        grid[
            [
                "ba_code",
                "fuel_category",
                "datetime_utc",
                "report_date",
                "profile",
                "flat_profile",
                "profile_method",
                "local_timezone_offset",
            ]
        ],
        group_keys=["ba_code", "fuel_category", "local_timezone_offset"],
        columns_to_convert=["profile", "flat_profile"],
    )


@pytest.fixture
def monthly_fleet_data():
    """Monthly EIA-923 data for the fleet in `zero_january_profiles`."""
    return pd.DataFrame(
        {
            "ba_code": "AECI",
            "fuel_category": "natural_gas",
            "report_date": pd.date_range("2025-01-01", periods=12, freq="MS"),
            "net_generation_mwh": 100.0,
            "fuel_consumed_mmbtu": 50.0,
            "local_timezone_offset": _offset("US/Central"),
        }
    )


def test_shape_monthly_eia_data_as_hourly_flat_profile_for_zero_month(
    zero_january_profiles, monthly_fleet_data
):
    """A month whose profile is zero in every hour should be shaped with a flat
    profile, rather than dropped.

    A profile of zeros gives blank shares of the monthly total (0 / 0), which used to
    make the month's generation disappear from the hourly data.
    """
    shaped = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        monthly_fleet_data, zero_january_profiles, 2025
    )

    # no generation is lost
    assert shaped["net_generation_mwh"].sum() == pytest.approx(1200.0)
    january = shaped[shaped["report_date"].dt.month == 1]
    # January is spread evenly over its hours, and flagged as such
    assert january["net_generation_mwh"].sum() == pytest.approx(100.0)
    assert january["net_generation_mwh"].nunique() == 1
    assert (january["profile_method"] == "flat_zero_profile").all()
    # the other months keep their own profile
    assert (
        shaped.loc[shaped["report_date"].dt.month != 1, "profile_method"]
        == "residual_profile"
    ).all()


def test_validate_shaped_totals_raises_on_lost_months(
    zero_january_profiles, monthly_fleet_data
):
    """Shaped data that lose part of the monthly totals should fail validation.

    The check used to raise only if shaping added to the totals, so months lost in
    shaping passed.
    """
    shaped = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        monthly_fleet_data, zero_january_profiles, 2025
    )
    group_keys = ["ba_code", "fuel_category"]
    # complete data pass
    validation.validate_shaped_totals(shaped, monthly_fleet_data, 2025, group_keys)
    # a month with no shaped data fails
    with pytest.raises(UserWarning):
        validation.validate_shaped_totals(
            shaped[shaped["report_date"].dt.month != 1],
            monthly_fleet_data,
            2025,
            group_keys,
        )
    # so does a month whose shaped generation is too low
    lower = shaped.copy()
    lower.loc[lower["report_date"].dt.month == 2, "net_generation_mwh"] *= 0.5
    with pytest.raises(UserWarning):
        validation.validate_shaped_totals(lower, monthly_fleet_data, 2025, group_keys)


def _hourly_series(ba, timezone, months, values, columns, extra=None):
    """Hourly rows for `ba` in the given months of 2025, in its local time.

    `values` maps each column in `columns` to a function of the local hour of the day
    that returns the value for that hour.
    """
    data = helpers.create_local_year_timestamps(2025, timezone)
    data = data[data["report_date"].dt.month.isin(months)].copy()
    local_time = data["datetime_utc"].dt.tz_convert(timezone)
    data["ba_code"] = ba
    data["datetime_local"] = local_time.astype(str)
    data["report_date"] = data["report_date"].astype("datetime64[s]")
    for column in columns:
        data[column] = local_time.dt.hour.map(values[column]).astype(float)
    for column, value in (extra or {}).items():
        data[column] = value
    return data.reset_index(drop=True)


def test_assign_storage_shaping_category():
    """Storage subplants should get the storage shaping category of their prime mover
    and plant, and other subplants should keep their fuel category.

    - pumped storage gets "pumped_storage"
    - a battery at a co-located solar plant gets "solar_w_integrated_battery_storage"
    - a battery at a hybrid wind plant gets "wind_w_integrated_battery_storage"
    - a standalone battery gets "battery_storage"
    """
    primary_fuel_table = pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 2, 3, 3, 4],
            "subplant_id": [1, 1, 2, 1, 2, 1],
            "energy_source_code": ["WAT", "MWH", "SUN", "MWH", "WND", "MWH"],
            "subplant_primary_prime_mover_code": ["PS", "BA", "PV", "BA", "WT", "BA"],
            "plant_storage_category": [
                pd.NA,
                "co_located",
                "co_located",
                "hybrid",
                "hybrid",
                "standalone",
            ],
        }
    )
    subplant_data = primary_fuel_table[["plant_id_eia", "subplant_id"]].assign(
        fuel_category=["storage", "storage", "solar", "storage", "wind", "storage"]
    )

    result = impute_hourly_profiles.assign_storage_shaping_category(
        subplant_data, primary_fuel_table
    )

    assert result["fuel_category"].to_list() == [
        "pumped_storage",
        "solar_w_integrated_battery_storage",
        "solar",
        "wind_w_integrated_battery_storage",
        "wind",
        "battery_storage",
    ]


def test_calculate_eia930_storage_series_splits_each_energy_source():
    """Discharging and charging should be split for each energy source before they're
    added up, so that a BA that reports charging separately keeps both directions.

    ERCO reports battery discharging as battery storage and charging as unknown energy
    storage. In the first hour, adding the two first would give a net discharge of 20
    MWh instead of 105 MWh of discharge (including other energy storage) and 80 MWh of
    charge. Batteries still reported in "other" (before the BA's cutover) are left out.
    """
    hours = pd.to_datetime(["2025-07-01 00:00", "2025-07-01 01:00"], utc=True)
    eia930_generation = pd.DataFrame(
        {
            "ba_code": "ERCO",
            "generation_energy_source": [
                "battery_storage",
                "unknown_energy_storage",
                "other_energy_storage",
                "battery_storage",
                "unknown_energy_storage",
                "pumped_storage",
                "battery_storage",
            ],
            "fuel_category_eia930": [
                "battery_storage",
                "battery_storage",
                "other_energy_storage",
                "battery_storage",
                "battery_storage",
                "pumped_storage",
                "other",
            ],
            "datetime_utc": hours[[0, 0, 0, 1, 1, 1, 1]],
            "datetime_local": "",
            "report_date": pd.Timestamp("2025-06-01"),
            "net_generation_mwh_930": [100.0, -80.0, 5.0, -10.0, -50.0, 30.0, 999.0],
        }
    )

    result = impute_hourly_profiles.calculate_eia930_storage_series(
        eia930_generation
    ).set_index(["fuel_category", "datetime_utc"])

    battery = result.loc["battery_storage"]
    assert battery["storage_discharge_profile"].to_list() == [105.0, 0.0]
    assert battery["storage_charge_profile"].to_list() == [80.0, 60.0]
    pumped = result.loc["pumped_storage"]
    assert pumped["storage_discharge_profile"].to_list() == [30.0]
    assert pumped["storage_charge_profile"].to_list() == [0.0]
    assert (result["method"] == "eia930_profile").all()


def test_calculate_integrated_storage_series_backs_out_battery():
    """The battery's discharging and charging should be recovered from solar with
    integrated battery storage, by subtracting the plain solar profile scaled to the
    solar part's monthly energy.

    In January, the combined series is half of the plain solar profile plus a battery
    that charges 30 MW at 10-11h and discharges 25 MW at 19-20h. The battery loses 10
    MWh a day, which EIA-923 reports, so the scaling factor works out to exactly 0.5
    and the battery is recovered exactly. In February, the combined series adds up to
    less than zero, so no battery series is backed out.
    """
    timezone = "US/Pacific"

    def solar(hour):
        return 100.0 if 8 <= hour <= 16 else 0.0

    def battery(hour):
        return {10: -30.0, 11: -30.0, 19: 25.0, 20: 25.0}.get(hour, 0.0)

    plain = _hourly_series(
        "LDWP",
        timezone,
        [1, 2],
        {"net_generation_mwh_930": solar},
        ["net_generation_mwh_930"],
        extra={"generation_energy_source": "solar", "fuel_category_eia930": "solar"},
    )
    combined = _hourly_series(
        "LDWP",
        timezone,
        [1, 2],
        {"net_generation_mwh_930": lambda hour: 0.5 * solar(hour) + battery(hour)},
        ["net_generation_mwh_930"],
        extra={
            "generation_energy_source": "solar_w_integrated_battery_storage",
            "fuel_category_eia930": "solar_w_integrated_battery_storage",
        },
    )
    # February's combined series is negative
    combined.loc[combined["report_date"].dt.month == 2, "net_generation_mwh_930"] = -1.0
    monthly_storage_by_ba = pd.DataFrame(
        {
            "ba_code": "LDWP",
            "fuel_category": "solar_w_integrated_battery_storage",
            "report_date": pd.to_datetime(["2025-01-01"]).astype("datetime64[s]"),
            "storage_charge_mwh": [60.0 * 31],
            "storage_discharge_mwh": [50.0 * 31],
        }
    )

    result = impute_hourly_profiles.calculate_integrated_storage_series(
        pd.concat([plain, combined], ignore_index=True), monthly_storage_by_ba
    )

    # only January is backed out
    assert set(result["report_date"].dt.month) == {1}
    local_hour = result["datetime_utc"].dt.tz_convert(timezone).dt.hour
    expected = local_hour.map(battery)
    assert result["storage_discharge_profile"].to_numpy() == pytest.approx(
        expected.clip(lower=0).to_numpy()
    )
    assert result["storage_charge_profile"].to_numpy() == pytest.approx(
        (-expected).clip(lower=0).to_numpy()
    )
    assert (result["fuel_category"] == "solar_w_integrated_battery_storage").all()
    assert (result["method"] == "integrated_residual_profile").all()


def test_calculate_storage_series_from_other_sizes_with_eia923_ratio():
    """For a BA whose "other" is mostly battery storage, the batteries should be backed
    out of "other" by removing a flat offset sized so that the ratio of discharging to
    charging matches EIA-923.

    CISO's "other" is a flat 700 MW of other resources plus batteries that charge 400 MW
    at 10-11h and discharge 340 MW at 18-19h, a ratio of 0.85, which EIA-923 also reports. So the batteries are
    recovered exactly in January. In February, CISO reports battery storage
    separately, so "other" isn't used. SCEG's "other" isn't used at all, since only
    CISO's "other" is used for battery storage.
    """
    timezone = "US/Pacific"

    def battery(hour):
        return {10: -400.0, 11: -400.0, 18: 340.0, 19: 340.0}.get(hour, 0.0)

    other = pd.concat(
        [
            _hourly_series(
                ba,
                timezone,
                [1, 2],
                {"net_generation_mwh_930": lambda hour: 700.0 + battery(hour)},
                ["net_generation_mwh_930"],
                extra={
                    "generation_energy_source": "other",
                    "fuel_category_eia930": "other",
                },
            )
            for ba in ["CISO", "SCEG"]
        ],
        ignore_index=True,
    )
    february = other[
        (other["ba_code"] == "CISO") & (other["report_date"].dt.month == 2)
    ]
    eia930_storage_series = february[
        ["ba_code", "datetime_utc", "datetime_local", "report_date"]
    ].assign(
        fuel_category="battery_storage",
        storage_discharge_profile=1.0,
        storage_charge_profile=1.0,
        method="eia930_profile",
    )
    report_dates = pd.to_datetime(["2025-01-01", "2025-02-01"]).astype("datetime64[s]")
    monthly_storage_by_ba = pd.DataFrame(
        {
            "ba_code": ["CISO", "CISO", "SCEG", "SCEG"],
            "fuel_category": "battery_storage",
            "report_date": report_dates[[0, 1, 0, 1]],
            "storage_charge_mwh": 1000.0,
            "storage_discharge_mwh": 850.0,
        }
    )
    result = impute_hourly_profiles.calculate_storage_series_from_other(
        other, eia930_storage_series, monthly_storage_by_ba
    )

    # only CISO in January
    assert set(result["ba_code"]) == {"CISO"}
    assert set(result["report_date"].dt.month) == {1}
    local_hour = result["datetime_utc"].dt.tz_convert(timezone).dt.hour
    expected = local_hour.map(battery)
    assert result["storage_discharge_profile"].to_numpy() == pytest.approx(
        expected.clip(lower=0).to_numpy(), abs=1e-6
    )
    assert result["storage_charge_profile"].to_numpy() == pytest.approx(
        (-expected).clip(lower=0).to_numpy(), abs=1e-6
    )
    assert (result["method"] == "other_residual_profile").all()


def test_identify_usable_storage_series():
    """Only months of storage series that look like storage and represent enough of
    the storage they shape should be used.

    In January 2025:
    - ERCO's batteries charge 100 MW at noon and discharge 85 MW at 19h every day,
      like its EIA-923 batteries: used in both directions
    - ISNE's pumped storage only reports discharging, like its EIA-923 discharge: its
      discharging is used, and its charging is blank
    - CPLE's batteries cycle like ERCO's, but at a tenth of their EIA-923 volume:
      not used
    - PSCO's batteries discharge much more than they charge: not used
    - FPC's batteries charge on only two days: not used
    - PNM's batteries discharge only half of what they charge (e.g. because some
      batteries charge while others discharge): used
    """
    timezone = "US/Eastern"
    shapes = {
        "ERCO": (
            lambda hour: 85.0 if hour == 19 else 0.0,
            lambda hour: 100.0 if hour == 12 else 0.0,
        ),
        "ISNE": (lambda hour: 85.0 if hour == 19 else 0.0, lambda hour: 0.0),
        "CPLE": (
            lambda hour: 8.5 if hour == 19 else 0.0,
            lambda hour: 10.0 if hour == 12 else 0.0,
        ),
        "PSCO": (
            lambda hour: 85.0 if hour == 19 else 0.0,
            lambda hour: 10.0 if hour == 12 else 0.0,
        ),
        "FPC": (lambda hour: 85.0 if hour == 19 else 0.0, lambda hour: 0.0),
        "PNM": (
            lambda hour: 50.0 if hour == 19 else 0.0,
            lambda hour: 100.0 if hour == 12 else 0.0,
        ),
    }
    storage_series = []
    for ba, (discharge, charge) in shapes.items():
        series = _hourly_series(
            ba,
            timezone,
            [1],
            {"storage_discharge_profile": discharge, "storage_charge_profile": charge},
            impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
            extra={
                "fuel_category": "pumped_storage"
                if ba == "ISNE"
                else "battery_storage",
                "method": "eia930_profile",
            },
        )
        storage_series.append(series)
    storage_series = pd.concat(storage_series, ignore_index=True)
    # FPC charges 1,240 MWh in total, but only on two days
    fpc_charge = (
        (storage_series["ba_code"] == "FPC")
        & storage_series["datetime_local"].str[:10].isin(["2025-01-05", "2025-01-06"])
        & storage_series["datetime_local"]
        .str[11:13]
        .isin(["10", "11", "12", "13", "14"])
    )
    storage_series.loc[fpc_charge, "storage_charge_profile"] = 124.0
    monthly_storage_to_shape = pd.DataFrame(
        {
            "ba_code": ["ERCO", "ISNE", "CPLE", "PSCO", "FPC", "PNM"],
            "fuel_category": [
                "battery_storage",
                "pumped_storage",
                "battery_storage",
                "battery_storage",
                "battery_storage",
                "battery_storage",
            ],
            "report_date": pd.Timestamp("2025-01-01"),
            "storage_charge_mwh": [3100.0, 3100.0, 3100.0, 3100.0, 1240.0, 3100.0],
            "storage_discharge_mwh": [2635.0, 2635.0, 2635.0, 2635.0, 2635.0, 2635.0],
        }
    )
    monthly_storage_to_shape["report_date"] = monthly_storage_to_shape[
        "report_date"
    ].astype("datetime64[s]")

    result = impute_hourly_profiles.identify_usable_storage_series(
        storage_series, monthly_storage_to_shape
    )

    used = result.groupby("ba_code")[
        impute_hourly_profiles.STORAGE_PROFILE_COLUMNS
    ].apply(lambda month: month.notna().all())
    assert used.to_dict(orient="index") == {
        "CPLE": {"storage_discharge_profile": False, "storage_charge_profile": False},
        "ERCO": {"storage_discharge_profile": True, "storage_charge_profile": True},
        "FPC": {"storage_discharge_profile": False, "storage_charge_profile": False},
        "ISNE": {"storage_discharge_profile": True, "storage_charge_profile": False},
        "PNM": {"storage_discharge_profile": True, "storage_charge_profile": True},
        "PSCO": {"storage_discharge_profile": False, "storage_charge_profile": False},
    }


def test_select_storage_profiles_hierarchy():
    """Each direction should use the best series with data in that month: the BA's
    own, then (for integrated storage) the BA's battery series, then the net-load
    proxy, then the national average, then flat.

    In January 2025, AZPS's batteries discharge at 19h but report no charging, and
    SRP's pumped storage reports both directions. SRP and WALC have net-load proxies
    (discharging at 20h and charging at 12h), and TEPC has none (e.g. no demand
    data). So:
    - SRP pumped storage uses its own series for both directions
    - AZPS integrated storage uses AZPS's battery discharging, and charging from its
      net-load proxy (none here), the national average (no BA reports battery
      charging) and then flat
    - SRP and WALC batteries use their net-load proxies in both directions
    - TEPC batteries use the national average discharging, which is AZPS's, and flat
      charging
    """
    timezone = "US/Arizona"
    keys = ["ba_code", "fuel_category"]
    azps_battery = _hourly_series(
        "AZPS",
        timezone,
        [1],
        {
            "storage_discharge_profile": lambda hour: 10.0 if hour == 19 else 0.0,
            "storage_charge_profile": lambda hour: 0.0,
        },
        impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
        extra={"fuel_category": "battery_storage", "method": "eia930_profile"},
    )
    srp_pumped = _hourly_series(
        "SRP",
        timezone,
        [1],
        {
            "storage_discharge_profile": lambda hour: 5.0 if hour == 18 else 0.0,
            "storage_charge_profile": lambda hour: 7.0 if hour == 3 else 0.0,
        },
        impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
        extra={"fuel_category": "pumped_storage", "method": "eia930_profile"},
    )
    storage_series = pd.concat([azps_battery, srp_pumped], ignore_index=True)
    net_load_proxy = pd.concat(
        [
            _hourly_series(
                ba,
                timezone,
                [1],
                {
                    "storage_discharge_profile": lambda hour: (
                        3.0 if hour == 20 else 0.0
                    ),
                    "storage_charge_profile": lambda hour: 4.0 if hour == 12 else 0.0,
                },
                impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
            )
            for ba in ["SRP", "WALC"]
        ],
        ignore_index=True,
    )
    storage_fleet_months = pd.DataFrame(
        {
            "ba_code": ["SRP", "AZPS", "SRP", "TEPC", "WALC"],
            "fuel_category": [
                "pumped_storage",
                "solar_w_integrated_battery_storage",
                "battery_storage",
                "battery_storage",
                "battery_storage",
            ],
            "report_date": pd.Timestamp("2025-01-01"),
        }
    )
    storage_fleet_months["report_date"] = storage_fleet_months["report_date"].astype(
        "datetime64[s]"
    )

    result = impute_hourly_profiles.select_storage_profiles(
        storage_fleet_months, storage_series, net_load_proxy, 2025
    )

    methods = (
        result.groupby(keys)[["profile_method", "storage_charge_profile_method"]]
        .first()
        .apply(tuple, axis=1)
        .to_dict()
    )
    assert methods == {
        ("AZPS", "solar_w_integrated_battery_storage"): (
            "eia930_profile",
            "assumed_flat",
        ),
        ("SRP", "battery_storage"): ("net_load_proxy", "net_load_proxy"),
        ("SRP", "pumped_storage"): ("eia930_profile", "eia930_profile"),
        ("TEPC", "battery_storage"): ("national_average", "assumed_flat"),
        ("WALC", "battery_storage"): ("net_load_proxy", "net_load_proxy"),
    }
    # every fleet has every hour of January
    assert (result.groupby(keys).size() == 744).all()
    local_hour = result["datetime_utc"].dt.tz_convert(timezone).dt.hour
    # AZPS and TEPC discharge at 19h (AZPS's series), SRP and WALC batteries at 20h
    # (their proxies)
    for bas, hour, value in [(["AZPS", "TEPC"], 19, 10.0), (["SRP", "WALC"], 20, 3.0)]:
        x = result[
            result["ba_code"].isin(bas) & (result["fuel_category"] != "pumped_storage")
        ]
        assert (
            x["storage_discharge_profile"]
            == local_hour[x.index].map(lambda h: value if h == hour else 0.0)
        ).all()
    # the discharging profile also shapes the other data columns
    assert (result["profile"] == result["storage_discharge_profile"]).all()


def test_calculate_net_load_storage_proxy():
    """Storage should discharge in proportion to how far each hour's net load (demand
    minus solar and wind) is above the local day's average, and charge in proportion
    to how far below.

    CISO's demand is 100 MW in every hour, with 60 MW of solar from 10h to 15h (local
    time), and LDWP has demand but no solar or wind. On each day, CISO's net load
    averages 85 MW, so it discharges 15 MW in the 18 hours without solar and charges
    45 MW in the 6 hours with solar. LDWP's net load is flat, so it neither charges
    nor discharges. BANC has no demand data, so it has no proxy.
    """
    timezone = "US/Pacific"
    generation = pd.concat(
        [
            _hourly_series(
                ba,
                timezone,
                [1],
                {"net_generation_mwh_930": solar},
                ["net_generation_mwh_930"],
                extra={"generation_energy_source": source},
            )
            for ba, source, solar in [
                ("CISO", "solar", lambda hour: 60.0 if 10 <= hour <= 15 else 0.0),
                ("LDWP", "natural_gas", lambda hour: 100.0),
                ("BANC", "solar", lambda hour: 10.0),
            ]
        ],
        ignore_index=True,
    )
    operations = (
        generation[generation["ba_code"] != "BANC"][["ba_code", "datetime_utc"]]
        .rename(columns={"ba_code": "balancing_authority_code_eia"})
        .assign(demand_mwh=100.0)
    )

    result = impute_hourly_profiles.calculate_net_load_storage_proxy(
        operations, generation
    )

    assert set(result["ba_code"]) == {"CISO", "LDWP"}
    ciso = result[result["ba_code"] == "CISO"]
    local_hour = ciso["datetime_utc"].dt.tz_convert(timezone).dt.hour
    is_solar_hour = local_hour.between(10, 15)
    assert ciso.loc[~is_solar_hour, "storage_discharge_profile"].to_numpy() == (
        pytest.approx(15.0)
    )
    assert ciso.loc[is_solar_hour, "storage_charge_profile"].to_numpy() == (
        pytest.approx(45.0)
    )
    assert (ciso.loc[is_solar_hour, "storage_discharge_profile"] == 0).all()
    ldwp = result[result["ba_code"] == "LDWP"]
    assert (ldwp[impute_hourly_profiles.STORAGE_PROFILE_COLUMNS] == 0).all(axis=None)


@pytest.fixture
def storage_profiles_for_shaping():
    """Hourly profiles (as shares of the month) for an AECI battery fleet in January and
    February 2025, which discharges at 19h and charges at 3h."""
    profiles = _hourly_series(
        "AECI",
        "US/Central",
        [1, 2],
        {
            "storage_discharge_profile": lambda hour: 1.0 if hour == 19 else 0.0,
            "storage_charge_profile": lambda hour: 1.0 if hour == 3 else 0.0,
        },
        impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
        extra={
            "fuel_category": "battery_storage",
            "flat_profile": 1.0,
            "profile_method": "eia930_profile",
            "local_timezone_offset": _offset("US/Central"),
        },
    )
    profiles["profile"] = profiles["storage_discharge_profile"]
    return impute_hourly_profiles.convert_profile_to_percent(
        profiles,
        group_keys=["ba_code", "fuel_category", "local_timezone_offset"],
        columns_to_convert=["profile", "flat_profile"]
        + impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
    )


def test_shape_monthly_eia_data_as_hourly_storage_keeps_monthly_net(
    storage_profiles_for_shaping,
):
    """Storage should be shaped as the shaped discharge minus the shaped charge, which
    keeps the monthly net generation.

    In January, the net generation (-100 MWh) is lower than the discharge (310 MWh)
    minus the charge (372 MWh), so the 38 MWh difference is added to the charge, and
    shaped with the charging profile (at 3h). In February, there is no charging or
    discharging data, so the net generation is shaped as before (flat, since it's
    negative).
    """
    monthly_data = pd.DataFrame(
        {
            "ba_code": "AECI",
            "fuel_category": "battery_storage",
            "report_date": pd.to_datetime(["2025-01-01", "2025-02-01"]).astype(
                "datetime64[s]"
            ),
            "net_generation_mwh": [-100.0, -56.0],
            "fuel_consumed_mmbtu": 0.0,
            "storage_charge_mwh": [372.0, np.nan],
            "storage_discharge_mwh": [310.0, np.nan],
            "local_timezone_offset": _offset("US/Central"),
        }
    )

    shaped = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        monthly_data, storage_profiles_for_shaping, 2025
    )

    monthly_totals = shaped.groupby("report_date")["net_generation_mwh"].sum()
    assert monthly_totals.to_list() == pytest.approx([-100.0, -56.0])
    january = shaped[shaped["report_date"].dt.month == 1]
    local_hour = january["datetime_utc"].dt.tz_convert("US/Central").dt.hour
    expected = local_hour.map({19: 310.0 / 31, 3: -(372.0 + 38.0) / 31}).fillna(0)
    assert january["net_generation_mwh"].to_numpy() == pytest.approx(
        expected.to_numpy()
    )
    assert (january["profile_method"] == "eia930_profile").all()
    february = shaped[shaped["report_date"].dt.month == 2]
    assert february["net_generation_mwh"].nunique() == 1
    assert (february["profile_method"] == "flat_negative_generation").all()


def test_residual_cems_switches_category_in_ba_local_month(monkeypatch):
    """CEMS generation should switch to a new EIA-930 energy source in the BA's local
    month, as EIA-930 does, not the plant's.

    A battery in Georgia (US/Eastern) is in SOCO (US/Central), which started to report
    battery storage on 2025-05-15, so from SOCO's June onwards. At 2025-06-01 04:00 UTC
    the plant is in June but SOCO is still in May, so that hour stays in "other".
    """
    hours = pd.date_range("2025-06-01 03:00", periods=3, freq="h", tz="UTC")
    cems = pd.DataFrame(
        {
            "plant_id_eia": 1,
            "subplant_id": 1,
            "datetime_utc": hours,
            # the plant's local month
            "report_date": pd.to_datetime(["2025-05-01", "2025-06-01", "2025-06-01"]),
            "net_generation_mwh": [1.0, 2.0, 4.0],
        }
    )

    def assign_fleet(subplant_data, fuel_category_col, **_):
        return subplant_data.assign(
            ba_code="SOCO",
            distribution_flag=False,
            **{
                fuel_category_col: "battery_storage"
                if fuel_category_col == "fuel_category_eia930_expanded"
                else "other"
            },
        )

    monkeypatch.setattr(
        impute_hourly_profiles, "assign_fleet_to_subplant_data", assign_fleet
    )
    cutover_dates = pd.DataFrame(
        {
            "ba_code": ["SOCO"],
            "generation_energy_source": ["battery_storage"],
            "cutover_date": [pd.Timestamp("2025-05-15")],
        }
    )

    cems_agg = impute_hourly_profiles.aggregate_cems_to_fleet_for_residual_calc(
        cems,
        cems.iloc[0:0],
        cems.iloc[0:0],
        plant_attributes=None,
        primary_fuel_table=None,
        year=2025,
        cutover_dates=cutover_dates,
    )

    assert cems_agg.set_index("datetime_utc")["fuel_category"].to_dict() == {
        # 2025-05-31 22:00 in SOCO
        hours[0]: "other",
        # 2025-05-31 23:00 in SOCO, but already June at the plant
        hours[1]: "other",
        # 2025-06-01 00:00 in SOCO
        hours[2]: "battery_storage",
    }
    assert cems_agg.set_index("datetime_utc")["cems_profile"].to_dict() == {
        hours[0]: 1.0,
        hours[1]: 2.0,
        hours[2]: 4.0,
    }


def test_add_storage_data_to_monthly_eia_data_fills_pumped_storage_discharge():
    """The blank discharge of pumped storage with inflow should be filled with the
    discharge implied by the net generation and charge (but not less than zero), and
    subplants without storage data should stay blank.
    """
    report_date = pd.Timestamp("2025-01-01")
    monthly_eia_data_to_shape = pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 3, 4],
            "subplant_id": 1,
            "report_date": report_date,
            "net_generation_mwh": [-20.0, -150.0, 50.0, 10.0],
        }
    )
    monthly_storage_data = pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 3],
            "subplant_id": 1,
            "report_date": report_date,
            "storage_charge_mwh": [100.0, 100.0, 40.0],
            "storage_discharge_mwh": [np.nan, np.nan, 90.0],
        }
    )

    result = impute_hourly_profiles.add_storage_data_to_monthly_eia_data(
        monthly_eia_data_to_shape, monthly_storage_data
    )

    # filled for pumped storage, clipped at zero, kept where reported, blank otherwise
    assert result["storage_discharge_mwh"].to_list()[:3] == [80.0, 0.0, 90.0]
    assert result[["storage_charge_mwh", "storage_discharge_mwh"]].iloc[3].isna().all()


def test_select_storage_profiles_ignores_national_average_at_year_edges():
    """A month in which no BA has its own series should get a flat profile, rather
    than the few year-edge hours that the national average fills from the next month.

    Only AZPS's February data is available, so the national average only has values
    in the first hours of January (filled from February at the edge of the year).
    Using them would put all of TEPC's January charging and discharging into those
    hours.
    """
    timezone = "US/Arizona"
    azps_battery = _hourly_series(
        "AZPS",
        timezone,
        [2],
        {
            "storage_discharge_profile": lambda hour: 10.0 if hour == 19 else 0.0,
            "storage_charge_profile": lambda hour: 10.0 if hour == 12 else 0.0,
        },
        impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
        extra={"fuel_category": "battery_storage", "method": "eia930_profile"},
    )
    storage_fleet_months = pd.DataFrame(
        {
            "ba_code": ["TEPC"],
            "fuel_category": ["battery_storage"],
            "report_date": pd.to_datetime(["2025-01-01"]).astype("datetime64[s]"),
        }
    )
    result = impute_hourly_profiles.select_storage_profiles(
        storage_fleet_months,
        azps_battery,
        pd.DataFrame(
            columns=["ba_code", "datetime_utc"]
            + impute_hourly_profiles.STORAGE_PROFILE_COLUMNS
        ),
        2025,
    )

    assert (result["profile_method"] == "assumed_flat").all()
    assert (result["storage_charge_profile_method"] == "assumed_flat").all()
    assert (result[impute_hourly_profiles.STORAGE_PROFILE_COLUMNS] == 1.0).all(
        axis=None
    )


def test_storage_fleet_shaping_matches_sum_of_subplants(
    storage_profiles_for_shaping, monkeypatch
):
    """Storage fleets should be shaped with the sum of their subplants' shaping
    volumes, so that the hourly fleet data is the sum of the hourly subplant data.

    Two batteries each discharge 100 MWh and charge 120 MWh. The first's net
    generation (-30 MWh) is 10 MWh lower than discharge minus charge, and the second's
    (-10 MWh) is 10 MWh higher. Each is shaped with 100/130 and 110/120 MWh, so the
    fleet should be shaped with 210/250 MWh, not with 200/240 MWh, where the
    differences cancel.
    """
    monthly_data = pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "subplant_id": 1,
            "ba_code": "AECI",
            "fuel_category": "battery_storage",
            "report_date": pd.to_datetime(["2025-01-01"] * 2).astype("datetime64[s]"),
            "net_generation_mwh": [-30.0, -10.0],
            "fuel_consumed_mmbtu": 0.0,
            "storage_charge_mwh": 120.0,
            "storage_discharge_mwh": 100.0,
            "local_timezone_offset": _offset("US/Central"),
        }
    )
    monkeypatch.setattr(
        impute_hourly_profiles,
        "assign_fleet_to_subplant_data",
        lambda data, *_, fuel_category_col, **__: data.assign(
            **{fuel_category_col: "battery_storage"}
        ),
    )
    monkeypatch.setattr(
        impute_hourly_profiles,
        "assign_storage_shaping_category",
        lambda data, *_, **__: data,
    )
    monkeypatch.setattr(
        impute_hourly_profiles,
        "get_shaped_plant_id_from_ba_fuel",
        lambda data: data.assign(shaped_plant_id=9001),
    )

    subplants = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        monthly_data, storage_profiles_for_shaping, 2025
    )
    fleet = impute_hourly_profiles.aggregate_eia_data_to_fleet(
        monthly_data.drop(columns="local_timezone_offset"), None, None, 2025
    ).assign(local_timezone_offset=_offset("US/Central"))
    fleet = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        fleet,
        storage_profiles_for_shaping,
        2025,
        fuel_category_col_for_shaping="fuel_category_for_shaping",
    )

    subplant_total = subplants.groupby("datetime_utc")["net_generation_mwh"].sum()
    fleet_total = fleet.set_index("datetime_utc")["net_generation_mwh"]
    assert fleet_total.to_numpy() == pytest.approx(
        subplant_total.loc[fleet_total.index].to_numpy()
    )
    local_hour = fleet_total.index.tz_convert("US/Central").hour
    assert fleet_total[local_hour == 19].sum() == pytest.approx(210.0)
    assert fleet_total[local_hour == 3].sum() == pytest.approx(-250.0)


def test_calculate_storage_profiles_without_storage_series(monkeypatch):
    """Without any EIA-930 storage series (e.g. before any BA reports storage), the
    profiles should come from the net-load proxy instead of failing."""
    timezone = "US/Arizona"
    generation = _hourly_series(
        "SRP",
        timezone,
        [1],
        {"net_generation_mwh_930": lambda hour: 5.0},
        ["net_generation_mwh_930"],
        extra={"fuel_category_eia930": "solar", "generation_energy_source": "solar"},
    )
    monthly_storage = pd.DataFrame(
        {
            "ba_code": ["SRP"],
            "fuel_category": ["battery_storage"],
            "report_date": pd.to_datetime(["2025-01-01"]).astype("datetime64[s]"),
            "storage_charge_mwh": [100.0],
            "storage_discharge_mwh": [90.0],
        }
    )
    net_load_proxy = _hourly_series(
        "SRP",
        timezone,
        range(1, 13),
        {
            "storage_discharge_profile": lambda hour: 3.0 if hour == 20 else 0.0,
            "storage_charge_profile": lambda hour: 4.0 if hour == 12 else 0.0,
        },
        impute_hourly_profiles.STORAGE_PROFILE_COLUMNS,
    )
    monkeypatch.setattr(
        impute_hourly_profiles,
        "assign_local_time_and_category_to_eia930_data",
        lambda *_: generation,
    )
    monkeypatch.setattr(
        impute_hourly_profiles,
        "aggregate_storage_data_to_ba",
        lambda *_: monthly_storage,
    )
    monkeypatch.setattr(
        impute_hourly_profiles,
        "calculate_net_load_storage_proxy",
        lambda *_: net_load_proxy,
    )

    profiles, other_storage_months = impute_hourly_profiles.calculate_storage_profiles(
        None, None, None, None, None, None, 2025
    )

    assert len(other_storage_months) == 0
    january = profiles[profiles["report_date"].dt.month == 1]
    assert (january["profile_method"] == "net_load_proxy").all()
    assert (january["storage_charge_profile_method"] == "net_load_proxy").all()
