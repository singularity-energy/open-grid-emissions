import os

import numpy as np
import pandas as pd
import pytest

import oge.physics_reconciliation as physics_reconciliation
from oge.physics_reconciliation import (
    DEMAND,
    DEMAND_LOWER_BOUND,
    GENERATION_BY_ENERGY_SOURCE_START_DATE,
    INTERCHANGE,
    NET_GENERATION,
    NET_GENERATION_BY_ENERGY_SOURCE,
    SERIES_LEVELS,
    TOTAL_INTERCHANGE,
    EIA930Tables,
    apply_basic_cleaning,
    apply_rolling_window_cleaning,
    eia930_tables_exist,
    find_zero_months,
    melt_series_to_eia930_tables,
    pivot_eia930_tables_to_series,
    read_eia930_tables,
    reconcile_eia930,
    reconcile_energy_balance,
    check_bas_start_or_stop_reporting,
    write_eia930_tables,
    zero_inactive_bas,
)

FIXTURE_FOLDER = os.path.join(os.path.dirname(__file__), "data")

# each fixture is two weeks of real EIA-930 data for every BA, identified by its dates
# (see test/data/README.md)
FIXTURE_DATES = ["2024-11-01_2024-11-14", "2018-06-24_2018-07-07"]

# stages of `reconcile_eia930` saved in the reconciled fixtures
RECONCILED_STAGES = ["reconciled_without_energy_sources", "reconciled"]

# reconciled values are saved rounded to this many decimals, which is well within
# the solver's precision
RECONCILED_DECIMALS = 3

# largest constraint violation, and largest distance from the exact solution, allowed
# for the reconciled data
CONSTRAINT_TOLERANCE = 1e-6  # MWh
SOLUTION_TOLERANCE = 0.01  # MWh


def melt_series_to_long(series: pd.DataFrame) -> pd.DataFrame:
    """Melts a series data frame to one row per series and hour.

    Args:
        series (pd.DataFrame): series data frame.

    Returns:
        pd.DataFrame: table with columns "datetime_utc", "field",
            "balancing_authority_code_eia", "other" and "value_mwh".
    """
    long = (
        series.stack(SERIES_LEVELS, future_stack=True).rename("value_mwh").reset_index()
    )
    return long[["datetime_utc", *SERIES_LEVELS, "value_mwh"]]


def pivot_long_to_series(long: pd.DataFrame) -> pd.DataFrame:
    """Pivots a table written by `melt_series_to_long` back to a series data frame.

    Args:
        long (pd.DataFrame): table with columns "datetime_utc", "field",
            "balancing_authority_code_eia", "other" and "value_mwh".

    Returns:
        pd.DataFrame: series data frame, with columns in sorted order.
    """
    series = long.pivot(index="datetime_utc", columns=SERIES_LEVELS, values="value_mwh")
    return series.sort_index(axis="columns")


def read_raw_fixture(dates: str) -> EIA930Tables:
    """Reads a raw EIA-930 fixture in `test/data` as EIA-930 tables.

    Args:
        dates (str): dates of the fixture, e.g. "2024-11-01_2024-11-14".

    Returns:
        EIA930Tables: EIA-930 data.
    """
    long = pd.read_parquet(os.path.join(FIXTURE_FOLDER, f"eia930_raw_{dates}.parquet"))
    return melt_series_to_eia930_tables(pivot_long_to_series(long))


def reconcile_fixture(dates: str, output_folder: str) -> pd.DataFrame:
    """Reconciles a raw fixture and returns the reconciled data of each stage.

    Args:
        dates (str): dates of the fixture, e.g. "2024-11-01_2024-11-14".
        output_folder (str): folder for the files written by `reconcile_eia930`.

    Returns:
        pd.DataFrame: table with a "stage" column and the columns of
            `melt_series_to_long`, with values rounded to `RECONCILED_DECIMALS`.
    """
    reconcile_eia930(read_raw_fixture(dates), output_folder)
    stages = []
    for stage in RECONCILED_STAGES:
        if not eia930_tables_exist(output_folder, stage):
            continue
        series = pivot_eia930_tables_to_series(read_eia930_tables(output_folder, stage))
        stages.append(melt_series_to_long(series).assign(stage=stage))
    reconciled = pd.concat(stages, ignore_index=True)
    reconciled["value_mwh"] = reconciled["value_mwh"].round(RECONCILED_DECIMALS)
    return reconciled[["stage", "datetime_utc", *SERIES_LEVELS, "value_mwh"]]


def write_reconciled_fixture(dates: str, output_folder: str):
    """Writes the reconciled fixture for a raw fixture to `test/data`.

    Run this to update the reconciled fixtures after an intended change to the
    reconciliation, and review the differences before committing them.

    Args:
        dates (str): dates of the fixture, e.g. "2024-11-01_2024-11-14".
        output_folder (str): folder for the files written by `reconcile_eia930`.
    """
    reconcile_fixture(dates, output_folder).to_parquet(
        os.path.join(FIXTURE_FOLDER, f"eia930_reconciled_{dates}.parquet"),
        index=False,
        compression="zstd",
    )


@pytest.mark.parametrize("dates", FIXTURE_DATES)
def test_reconcile_eia930_matches_reconciled_fixture(dates, tmp_path):
    """The reconciled data should not change unless the reconciliation is changed on
    purpose.

    Each raw fixture is two weeks of real EIA-930 data for every BA. Every BA is
    included because the optimization couples BAs through interchange, so dropping a BA
    would change the results for the others.

    - The 2024 fixture covers the first weeks of storage reporting in several BAs
      (e.g. ISNE, ERCO, NEVP), when negative storage values started to appear. It
      includes the energy sources that EIA-930 introduced in 2024 (e.g. battery
      storage, pumped storage) as their own series.
    - The 2018 fixture crosses 2018-07-01, so it exercises both the optimization
      without generation by energy source (before) and with it (after).

    The reconciled fixtures were written by `write_reconciled_fixture`. The test checks
    that the same series and hours are reconciled, to values within
    `SOLUTION_TOLERANCE`.
    """
    result = reconcile_fixture(dates, str(tmp_path))
    expected = pd.read_parquet(
        os.path.join(FIXTURE_FOLDER, f"eia930_reconciled_{dates}.parquet")
    )
    for stage in RECONCILED_STAGES:
        stage_result = result[result["stage"] == stage].drop(columns="stage")
        stage_expected = expected[expected["stage"] == stage].drop(columns="stage")
        pd.testing.assert_frame_equal(
            pivot_long_to_series(stage_result),
            pivot_long_to_series(stage_expected),
            check_exact=False,
            rtol=0,
            atol=SOLUTION_TOLERANCE,
            obj=stage,
        )


@pytest.mark.parametrize("dates", FIXTURE_DATES)
def test_reconcile_energy_balance_solves_the_problem(dates):
    """The reconciled data should satisfy every constraint and be the exact solution.

    This protects the formulation in `HourlyReconciliationProblem`. For each fixture, after basic and rolling window cleaning, the
    hours are reconciled as in `reconcile_eia930`: without generation by energy source
    before 2018-07-01 (only in the 2018 fixture) and with it from that date, and:
        - every constraint of the problem holds to `CONSTRAINT_TOLERANCE`
        - the reconciled data are within `SOLUTION_TOLERANCE` of a solve with much
          tighter tolerances, which shows the default solver settings are precise
          enough
    """
    series = pivot_eia930_tables_to_series(read_raw_fixture(dates))
    series, weights = apply_rolling_window_cleaning(apply_basic_cleaning(series))

    periods = {
        False: slice(None, GENERATION_BY_ENERGY_SOURCE_START_DATE),
        True: slice(GENERATION_BY_ENERGY_SOURCE_START_DATE, None),
    }
    for include_energy_sources, hours in periods.items():
        if len(series.loc[hours]) == 0:
            continue
        reconciled, _ = reconcile_energy_balance(
            series.loc[hours], weights.loc[hours], include_energy_sources
        )
        check_reconciliation_constraints(
            reconciled, series.loc[hours], include_energy_sources
        )

        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(
                physics_reconciliation,
                "CLARABEL_SETTINGS",
                {"tol_gap_abs": 1e-12, "tol_gap_rel": 1e-12, "tol_feas": 1e-12},
            )
            reference, _ = reconcile_energy_balance(
                series.loc[hours], weights.loc[hours], include_energy_sources
            )
        # every reconciled value is close to the exact solution
        pd.testing.assert_frame_equal(
            reconciled, reference, check_exact=False, rtol=0, atol=SOLUTION_TOLERANCE
        )


@pytest.fixture
def series_with_retired_ba():
    """Two months of series for GRIF, which retires on 2023-11-01 (07:00 UTC) and
    reports nothing after, and its neighbor WALC."""
    hours = pd.date_range("2023-10-01", periods=60 * 24, freq="h", tz="UTC")
    grif_keys = [
        (DEMAND, "GRIF", ""),
        (NET_GENERATION, "GRIF", ""),
        (TOTAL_INTERCHANGE, "GRIF", ""),
        (NET_GENERATION_BY_ENERGY_SOURCE, "GRIF", "gas"),
        (INTERCHANGE, "GRIF", "WALC"),
        (INTERCHANGE, "WALC", "GRIF"),
    ]
    values = {key: 500.0 for key in grif_keys}
    values[(DEMAND, "GRIF", "")] = 1.0
    values[(INTERCHANGE, "WALC", "GRIF")] = -500.0
    values[(NET_GENERATION, "WALC", "")] = 3000.0
    series = pd.DataFrame(values, index=hours)
    series.columns = pd.MultiIndex.from_tuples(series.columns, names=SERIES_LEVELS)
    series.loc[hours >= pd.Timestamp("2023-11-01 07:00", tz="UTC"), grif_keys] = np.nan
    return series, grif_keys


def test_zero_inactive_bas_stops_gap_filling_after_retirement(series_with_retired_ba):
    """A retired BA's series, and its neighbors' interchange with it, should be zero
    from its retirement, rather than filled with its last reported days.

    E.g. GRIF retired on 2023-11-01 and reports nothing after, but the data for 2024
    start in October 2023, and the gap filling carried its October values through all
    of 2024.
    """
    series, grif_keys = series_with_retired_ba
    retirement_date = pd.Timestamp("2023-11-01 07:00", tz="UTC")
    retired = series.index >= retirement_date
    walc_key = (NET_GENERATION, "WALC", "")
    # a neighbor still reports a small value for a few hours after the retirement
    series.loc[retirement_date, (INTERCHANGE, "WALC", "GRIF")] = -9.0

    zeroed = zero_inactive_bas(series, {"GRIF": (pd.NaT, retirement_date)})

    assert (zeroed.loc[retired, grif_keys] == 0).all(axis=None)
    pd.testing.assert_frame_equal(zeroed.loc[~retired], series.loc[~retired])
    pd.testing.assert_series_equal(zeroed[walc_key], series[walc_key])
    # without the zeros, the gap filling carries the October values forward
    filled, _ = apply_rolling_window_cleaning(series)
    assert (filled.loc[retired, (NET_GENERATION, "GRIF", "")] == 500).all()
    cleaned, _ = apply_rolling_window_cleaning(zeroed)
    assert (cleaned.loc[retired, grif_keys] == 0).all(axis=None)
    # BAs that are active throughout the data, or that aren't in the data, are
    # unchanged
    pd.testing.assert_frame_equal(
        zero_inactive_bas(
            series,
            {
                "GRIF": (pd.Timestamp("2004-01-01", tz="UTC"), pd.NaT),
                "WALC": (pd.NaT, series.index[-1] + pd.Timedelta("1h")),
                "GRMA": (pd.NaT, series.index[0]),
            },
        ),
        series,
    )


def test_zero_inactive_bas_stops_gap_filling_before_activation(series_with_retired_ba):
    """A BA's series should be zero before its activation, rather than filled with
    its first reported days, e.g. SIKE, activated on 2025-06-01."""
    series, grif_keys = series_with_retired_ba
    # reverse the hours, so that GRIF starts reporting a month into the data
    series = pd.DataFrame(
        series.to_numpy()[::-1], index=series.index, columns=series.columns
    )
    activation_date = (
        series.index[-1]
        - pd.Timedelta(pd.Timestamp("2023-11-01 07:00", tz="UTC") - series.index[0])
        + pd.Timedelta("1h")
    )
    inactive = series.index < activation_date
    assert series.loc[inactive, grif_keys].isna().all(axis=None)
    assert series.loc[~inactive, grif_keys].notna().all(axis=None)

    zeroed = zero_inactive_bas(series, {"GRIF": (activation_date, pd.NaT)})

    assert (zeroed.loc[inactive, grif_keys] == 0).all(axis=None)
    pd.testing.assert_frame_equal(zeroed.loc[~inactive], series.loc[~inactive])


def test_check_bas_start_or_stop_reporting(series_with_retired_ba, caplog):
    """A BA that stops reporting for the rest of the data should be reported, unless
    its retirement date explains it. The data only have 28 days after GRIF stops, so
    the shortest gap is lowered to 20 days."""
    series, _ = series_with_retired_ba
    with caplog.at_level("WARNING"):
        check_bas_start_or_stop_reporting(series, {}, min_gap=pd.Timedelta(days=20))
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert len(warnings) == 1
    assert warnings[0].startswith("GRIF stops reporting")

    caplog.clear()
    with caplog.at_level("WARNING"):
        check_bas_start_or_stop_reporting(
            series,
            {"GRIF": (pd.NaT, pd.Timestamp("2023-11-01 07:00", tz="UTC"))},
            min_gap=pd.Timedelta(days=20),
        )
    assert not [r for r in caplog.records if r.levelname == "WARNING"]


def test_find_zero_months_uses_local_months():
    """Zero months should be the months of each BA's local time.

    The series is zero for all of February in US/Pacific time, but not in UTC, since
    UTC February starts with 8 hours of local January. Hourly profiles are calculated
    for local months, so the months must match.
    """
    hours = pd.date_range("2025-01-31 00:00", "2025-03-01 23:00", freq="h", tz="UTC")
    local_month = hours.tz_convert("US/Pacific").month
    key = (NET_GENERATION_BY_ENERGY_SOURCE, "CISO", "oil")
    series = pd.DataFrame(
        {key: np.where(local_month == 2, 0.0, 5.0)},
        index=hours,
    )
    series.columns = pd.MultiIndex.from_tuples([key], names=SERIES_LEVELS)

    local = find_zero_months(series, {"CISO": "US/Pacific"})[key]
    # every hour of local February, and no other hour
    assert (local == (local_month == 2)).all()
    # in UTC, no month is zero in every hour
    assert not find_zero_months(series, {})[key].any()


def test_find_zero_months_skips_bas_without_energy_sources():
    """A BA that reports net generation but no generation by energy source for a
    whole month isn't reporting zero generation, so its energy sources shouldn't be
    kept at zero. A BA with no generation at all should.

    E.g. SRP in July 2018, before it reported generation by energy source. Keeping its
    energy sources at zero would force its net generation to zero.
    """
    hours = pd.date_range("2018-07-01", periods=48, freq="h", tz="UTC")
    keys = [
        (NET_GENERATION, "SRP", ""),
        (NET_GENERATION_BY_ENERGY_SOURCE, "SRP", "gas"),
        (NET_GENERATION_BY_ENERGY_SOURCE, "SRP", "solar"),
        (NET_GENERATION, "HST", ""),
        (NET_GENERATION_BY_ENERGY_SOURCE, "HST", "solar"),
    ]
    series = pd.DataFrame(
        [[6500.0, 0.0, 0.0, 0.0, 0.0]] * len(hours),
        index=hours,
        columns=pd.MultiIndex.from_tuples(keys, names=SERIES_LEVELS),
    )

    zero_months = find_zero_months(series, {})

    assert not zero_months.loc[:, (NET_GENERATION_BY_ENERGY_SOURCE, "SRP")].any(
        axis=None
    )
    assert zero_months[(NET_GENERATION_BY_ENERGY_SOURCE, "HST", "solar")].all()


def test_reconcile_energy_balance_keeps_zero_months_at_zero():
    """Generation by energy source that is zero for a whole month should stay exactly
    zero after reconciliation, while every other constraint still holds.

    Without this, the reconciliation adds small amounts of generation to energy
    sources a BA doesn't have, which are then used as hourly profiles.
    """
    series = pivot_eia930_tables_to_series(read_raw_fixture(FIXTURE_DATES[0]))
    series, weights = apply_rolling_window_cleaning(apply_basic_cleaning(series))
    zero_months = find_zero_months(series, {})
    # the fixture has energy sources that are zero for a whole month
    assert zero_months.any(axis=None)

    reconciled, _ = reconcile_energy_balance(
        series, weights, include_energy_sources=True, zero_months=zero_months
    )

    locked = zero_months.reindex(columns=reconciled.columns, fill_value=False)
    assert (reconciled.where(locked) == 0).sum().sum() == locked.sum().sum()
    check_reconciliation_constraints(reconciled, series, include_energy_sources=True)


def check_reconciliation_constraints(
    reconciled: pd.DataFrame, series: pd.DataFrame, include_energy_sources: bool
):
    """Asserts that reconciled data satisfy every constraint of the reconciliation.

    The constraints are checked directly on the series, independently of how
    `HourlyReconciliationProblem` builds them.

    Args:
        reconciled (pd.DataFrame): reconciled series data frame.
        series (pd.DataFrame): series data frame before reconciliation.
        include_energy_sources (bool): whether generation by energy source was
            reconciled.
    """
    fields = reconciled.columns.get_level_values("field")

    def by_ba(field: str) -> pd.DataFrame:
        values = reconciled.loc[:, fields == field]
        return values.T.groupby(level="balancing_authority_code_eia").sum().T

    demand = by_ba(DEMAND)
    net_generation = by_ba(NET_GENERATION)
    total_interchange = by_ba(TOTAL_INTERCHANGE)
    interchange = reconciled.loc[:, fields == INTERCHANGE]

    # demand + total interchange == net generation
    assert (
        (demand + total_interchange - net_generation).abs() < CONSTRAINT_TOLERANCE
    ).all(axis=None)
    # total interchange == sum of interchange with each adjacent BA
    interchange_sum = by_ba(INTERCHANGE).reindex(
        columns=total_interchange.columns, fill_value=0
    )
    assert ((total_interchange - interchange_sum).abs() < CONSTRAINT_TOLERANCE).all(
        axis=None
    )
    # interchange from A to B == -(interchange from B to A)
    opposite = reconciled.loc[
        :,
        [(INTERCHANGE, adjacent_ba, ba) for _, ba, adjacent_ba in interchange.columns],
    ]
    assert (
        abs(interchange.to_numpy() + opposite.to_numpy()) < CONSTRAINT_TOLERANCE
    ).all()
    # demand is at least the lower bound
    assert (demand >= DEMAND_LOWER_BOUND - CONSTRAINT_TOLERANCE).all(axis=None)
    # generation is at least the lower of zero and its value before reconciliation
    is_generation = fields.isin([NET_GENERATION, NET_GENERATION_BY_ENERGY_SOURCE])
    generation = reconciled.loc[:, is_generation]
    lower_bound = np.minimum(series.reindex(columns=generation.columns).fillna(0), 0)
    assert (generation >= lower_bound - CONSTRAINT_TOLERANCE).all(axis=None)

    # generation by energy source is only in the problem if it is reconciled
    assert (fields == NET_GENERATION_BY_ENERGY_SOURCE).any() == include_energy_sources
    if include_energy_sources:
        # net generation == sum of generation by energy source
        generation_sum = by_ba(NET_GENERATION_BY_ENERGY_SOURCE)
        assert ((net_generation - generation_sum).abs() < CONSTRAINT_TOLERANCE).all(
            axis=None
        )


@pytest.fixture
def tables_without_energy_sources():
    """EIA-930 tables with an empty generation by energy source table.

    Two hours before 2018-07-01 for one BA ("CISO") with interchange to a second BA
    ("BANC"), shaped like the "reconciled_without_energy_sources" tables written by
    `reconcile_eia930`.
    """
    hours = pd.to_datetime(["2018-06-30 00:00", "2018-06-30 01:00"], utc=True)
    return EIA930Tables(
        operations=pd.DataFrame(
            {
                "datetime_utc": hours,
                "balancing_authority_code_eia": "CISO",
                "demand_mwh": [100.0, 110.0],
                "net_generation_mwh": [90.0, 95.0],
                "interchange_mwh": [-10.0, -15.0],
            }
        ),
        net_generation_by_energy_source=pd.DataFrame(
            columns=[
                "datetime_utc",
                "balancing_authority_code_eia",
                "generation_energy_source",
                "net_generation_mwh",
            ]
        ),
        interchange=pd.DataFrame(
            {
                "datetime_utc": hours,
                "balancing_authority_code_eia": "CISO",
                "balancing_authority_code_adjacent_eia": "BANC",
                "interchange_mwh": [-10.0, -15.0],
            }
        ),
    )


def test_pivot_eia930_tables_with_empty_energy_sources(tables_without_energy_sources):
    """Tables with no generation by energy source should pivot and melt back.

    EIA-930 has no generation by energy source before 2018-07-01, so the
    "reconciled_without_energy_sources" tables have an empty generation table. Pivoting
    used to fail because the column levels couldn't be inferred from no columns.
    """
    series = pivot_eia930_tables_to_series(tables_without_energy_sources)
    # the column levels are kept, and there are no generation by energy source series
    assert list(series.columns.names) == SERIES_LEVELS
    assert (
        series.columns.get_level_values("field") == NET_GENERATION_BY_ENERGY_SOURCE
    ).sum() == 0
    # melting gives back the same tables
    result = melt_series_to_eia930_tables(series)
    assert len(result.net_generation_by_energy_source) == 0
    pd.testing.assert_frame_equal(
        result.operations, tables_without_energy_sources.operations, check_like=True
    )
    pd.testing.assert_frame_equal(
        result.interchange, tables_without_energy_sources.interchange, check_like=True
    )


def test_eia930_tables_exist_requires_every_table(
    tables_without_energy_sources, tmp_path
):
    """Cached tables should only count as written if every table's file exists.

    The pipeline reuses reconciled data when it exists. If a run was interrupted after
    writing some of the tables, the next run must reconcile again rather than fail
    when reading the missing table.
    """
    folder = str(tmp_path)
    assert not eia930_tables_exist(folder, "reconciled")
    write_eia930_tables(tables_without_energy_sources, folder, "reconciled")
    assert eia930_tables_exist(folder, "reconciled")
    # a run interrupted before the interchange table was written
    os.remove(os.path.join(folder, "eia930_reconciled_interchange.parquet"))
    assert not eia930_tables_exist(folder, "reconciled")
