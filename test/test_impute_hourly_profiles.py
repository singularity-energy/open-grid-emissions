from unittest.mock import patch

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


@pytest.fixture
def eia_only_storage_fleet():
    """Monthly EIA-923 data, plant attributes and primary fuels for one battery
    subplant in AECI that only reports to EIA, for January and February 2025.

    January has positive net generation and February negative net generation (more
    charging than discharging), as is common for storage. The residual profiles only
    have an AECI natural gas fleet, since EIA-930 has no storage category.
    """
    monthly_data = pd.DataFrame(
        {
            "plant_id_eia": 1,
            "subplant_id": 1,
            "report_date": pd.to_datetime(["2025-01-01", "2025-02-01"]).astype(
                "datetime64[s]"
            ),
            "net_generation_mwh": [744.0, -672.0],
            "fuel_consumed_mmbtu": 0.0,
            "fuel_consumed_for_electricity_mmbtu": 0.0,
        }
    )
    plant_attributes = pd.DataFrame({"plant_id_eia": [1], "ba_code": ["AECI"]})
    primary_fuel_table = pd.DataFrame(
        {
            "plant_id_eia": [1],
            "subplant_id": [1],
            "subplant_primary_fuel": ["MWH"],
            "subplant_primary_fuel_from_capacity_mw": ["MWH"],
            "subplant_primary_prime_mover_code": ["BA"],
        }
    )
    residual_profiles = helpers.create_local_year_timestamps(2025, "US/Central")
    residual_profiles["ba_code"] = "AECI"
    residual_profiles["datetime_local"] = (
        residual_profiles["datetime_utc"].dt.tz_convert("US/Central").astype(str)
    )
    residual_profiles["fuel_category"] = "natural_gas"
    residual_profiles["eia930_profile"] = 100.0
    residual_profiles["residual_profile"] = 100.0
    return monthly_data, plant_attributes, primary_fuel_table, residual_profiles


def test_impute_missing_hourly_profiles_flat_profile_for_eia_only_storage(
    eia_only_storage_fleet,
):
    """A storage fleet that only reports to EIA should get a flat imputed profile, and
    shaping it should keep its monthly totals, including a month of negative net
    generation.

    EIA-930 has no storage category, so every EIA-only storage fleet-month needs an
    imputed profile. There used to be no imputation method for the "storage" fuel
    category, which raised "Fuel category storage not recognized" in step 14.
    """
    monthly_data, plant_attributes, primary_fuel_table, residual_profiles = (
        eia_only_storage_fleet
    )
    with patch(
        "oge.impute_hourly_profiles.load_data.load_diba_data",
        return_value=pd.DataFrame(
            columns=[
                "ba_code",
                "ba_region",
                "diba_code",
                "diba_region",
                "timezone_local",
                "timezone_local_diba",
            ]
        ),
    ):
        hourly_profiles = impute_hourly_profiles.impute_missing_hourly_profiles(
            monthly_data, residual_profiles, plant_attributes, primary_fuel_table, 2025
        )

    storage_profiles = hourly_profiles[hourly_profiles["fuel_category"] == "storage"]
    # a flat profile for each month of storage data
    assert set(storage_profiles["report_date"].dt.month) == {1, 2}
    assert (storage_profiles["imputation_method"] == "assumed_flat").all()
    assert (storage_profiles["imputed_profile"] == 1.0).all()

    # shape the monthly data as in data_pipeline step 14 and
    # shape_monthly_eia_data_as_hourly
    storage_profiles = storage_profiles.assign(
        profile=storage_profiles["imputed_profile"],
        flat_profile=1.0,
        profile_method="assumed_flat",
        local_timezone_offset=_offset("US/Central"),
    )
    storage_profiles = impute_hourly_profiles.convert_profile_to_percent(
        storage_profiles,
        group_keys=["ba_code", "fuel_category", "local_timezone_offset"],
        columns_to_convert=["profile", "flat_profile"],
    )
    monthly_fleet_data = monthly_data.assign(
        ba_code="AECI",
        fuel_category="storage",
        local_timezone_offset=_offset("US/Central"),
    ).drop(columns=["plant_id_eia", "subplant_id"])
    shaped = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        monthly_fleet_data, storage_profiles, 2025
    )
    monthly_totals = shaped.groupby("report_date", dropna=False)[
        "net_generation_mwh"
    ].sum()
    # the monthly totals are kept
    assert monthly_totals.to_list() == pytest.approx([744.0, -672.0])
    # each month is spread evenly over its hours
    assert (shaped.groupby("report_date")["net_generation_mwh"].nunique() == 1).all()
    # the negative month is flagged as such
    february = shaped[shaped["report_date"].dt.month == 2]
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
