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
    melt_series_to_eia930_tables,
    pivot_eia930_tables_to_series,
    read_eia930_tables,
    reconcile_eia930,
    reconcile_energy_balance,
    write_eia930_tables,
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
      (e.g. ISNE, ERCO, NEVP), when negative storage values started to appear.
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
