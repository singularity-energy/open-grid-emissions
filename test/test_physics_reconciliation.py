import gzip
import io
import os
import re

import numpy as np
import pandas as pd
import pytest

import oge.physics_reconciliation as physics_reconciliation
from oge.physics_reconciliation import (
    DEMAND,
    GENERATION_BY_ENERGY_SOURCE_START_DATE,
    INTERCHANGE,
    NET_GENERATION,
    NET_GENERATION_BY_ENERGY_SOURCE,
    OPTIMIZATION_LOWER_BOUND,
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

# energy source codes used by gridemissions, mapped to the PUDL names used by OGE
ENERGY_SOURCE_CODES = {
    "COL": "coal",
    "NG": "gas",
    "NUC": "nuclear",
    "OIL": "oil",
    "OTH": "other",
    "SUN": "solar",
    "UNK": "unknown",
    "WAT": "hydro",
    "WND": "wind",
}

# gridemissions output file, mapped to the matching OGE stage
OUTPUT_STAGES = {
    "eia930_basic.csv": "basic",
    "eia930_rolling.csv": "rolling",
    "eia930_elec.csv": "reconciled",
}
PRE_2018_07_OUTPUT_STAGES = {
    "eia930_opt_no_src.csv": "reconciled_without_energy_sources",
}

# OGE solves the reconciliation with Clarabel at tight tolerances, while gridemissions
# uses OSQP at its default tolerances, which leaves errors of up to about 1.1 MWh (on
# the 2018 fixture) compared to the exact solution. The reconciled data are compared
# within these tolerances, which cover gridemissions' solver error but not a change to
# the problem itself.
RECONCILED_STAGES = ["reconciled", "reconciled_without_energy_sources"]
RECONCILED_TOLERANCE = 2.0  # MWh
OBJECTIVE_RELATIVE_TOLERANCE = 1e-4

# largest constraint violation, and largest distance from the exact solution, allowed
# for the reconciled data
CONSTRAINT_TOLERANCE = 1e-6  # MWh
SOLUTION_TOLERANCE = 0.01  # MWh


@pytest.fixture(scope="module")
def gridemissions_make_dataset(tmp_path_factory):
    """`make_dataset` from the original `gridemissions` package.

    `gridemissions` is a dev-only dependency, kept so that `oge.physics_reconciliation`
    can be checked against the code it was ported from. Importing it writes config and
    logging files, so its config folder is pointed at a temporary folder first. The
    test is skipped if the package isn't installed.
    """
    config_folder = tmp_path_factory.mktemp("gridemissions_config")
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GRIDEMISSIONS_CONFIG_DIR_PATH", str(config_folder))
        mp.setenv("GRIDEMISSIONS_CONFIG_FILE_PATH", str(config_folder / "config.json"))
        workflows = pytest.importorskip("gridemissions.workflows")
    return workflows.make_dataset


def read_raw_fixture(fixture_name: str) -> str:
    """Returns the text of a gzipped raw EIA-930 fixture in `test/data`.

    Args:
        fixture_name (str): file name of the fixture.

    Returns:
        str: CSV text, in the wide format used by gridemissions.
    """
    with gzip.open(os.path.join(FIXTURE_FOLDER, fixture_name), "rt") as f:
        return f.read()


def convert_wide_to_tables(wide: pd.DataFrame, value_suffix: str = "_mwh"):
    """Converts data in the wide gridemissions format to EIA-930 tables.

    This parses the gridemissions column names independently of
    `oge.physics_reconciliation`, so the test also checks how OGE builds its series.

    Args:
        wide (pd.DataFrame): one column per series, named like "EBA.{BA}-ALL.D.H".
        value_suffix (str, optional): suffix of the value columns. Defaults to "_mwh".

    Returns:
        EIA930Tables: EIA-930 data.
    """
    long = wide.rename_axis("datetime_utc").reset_index().melt(id_vars="datetime_utc")
    parts = long["variable"].str.extract(
        r"^EBA\.(?P<ba>[^-]+)-(?P<other>[^.]+)\.(?P<field>[A-Z]+)(?:\.(?P<src>[A-Z]+))?\.H$"
    )
    long = pd.concat([long, parts], axis="columns")

    operations_fields = {
        "D": f"demand{value_suffix}",
        "NG": f"net_generation{value_suffix}",
        "TI": f"interchange{value_suffix}",
    }
    is_operations = (long["other"] == "ALL") & long["src"].isna()
    operations = (
        long[is_operations]
        .assign(field=lambda df: df["field"].map(operations_fields))
        .pivot(index=["datetime_utc", "ba"], columns="field", values="value")
        .reset_index()
        .rename(columns={"ba": "balancing_authority_code_eia"})
    )
    operations.columns.name = None

    net_generation_by_energy_source = (
        long[long["src"].notna()]
        .assign(generation_energy_source=lambda df: df["src"].map(ENERGY_SOURCE_CODES))
        .rename(
            columns={
                "ba": "balancing_authority_code_eia",
                "value": f"net_generation{value_suffix}",
            }
        )[
            [
                "datetime_utc",
                "balancing_authority_code_eia",
                "generation_energy_source",
                f"net_generation{value_suffix}",
            ]
        ]
    )

    interchange = long[long["field"] == "ID"].rename(
        columns={
            "ba": "balancing_authority_code_eia",
            "other": "balancing_authority_code_adjacent_eia",
            "value": f"interchange{value_suffix}",
        }
    )[
        [
            "datetime_utc",
            "balancing_authority_code_eia",
            "balancing_authority_code_adjacent_eia",
            f"interchange{value_suffix}",
        ]
    ]
    return EIA930Tables(operations, net_generation_by_energy_source, interchange)


def convert_tables_to_wide(tables: EIA930Tables, value_suffix: str = "_mwh"):
    """Converts EIA-930 tables to the wide gridemissions format.

    Args:
        tables (EIA930Tables): EIA-930 data.
        value_suffix (str, optional): suffix of the value columns. Defaults to "_mwh".

    Returns:
        pd.DataFrame: one column per series, named like "EBA.{BA}-ALL.D.H", indexed
            by "datetime_utc", with columns in alphabetical order.
    """
    codes = {name: code for code, name in ENERGY_SOURCE_CODES.items()}
    frames = []
    for field, code in [
        ("demand", "D"),
        ("net_generation", "NG"),
        ("interchange", "TI"),
    ]:
        frame = tables.operations.pivot(
            index="datetime_utc",
            columns="balancing_authority_code_eia",
            values=f"{field}{value_suffix}",
        )
        frame.columns = [f"EBA.{ba}-ALL.{code}.H" for ba in frame.columns]
        frames.append(frame)
    if len(tables.net_generation_by_energy_source) > 0:
        frame = tables.net_generation_by_energy_source.pivot(
            index="datetime_utc",
            columns=["balancing_authority_code_eia", "generation_energy_source"],
            values=f"net_generation{value_suffix}",
        )
        frame.columns = [f"EBA.{ba}-ALL.NG.{codes[src]}.H" for ba, src in frame.columns]
        frames.append(frame)
    frame = tables.interchange.pivot(
        index="datetime_utc",
        columns=[
            "balancing_authority_code_eia",
            "balancing_authority_code_adjacent_eia",
        ],
        values=f"interchange{value_suffix}",
    )
    frame.columns = [f"EBA.{ba}-{adjacent}.ID.H" for ba, adjacent in frame.columns]
    frames.append(frame)

    wide = pd.concat(frames, axis="columns")
    wide.index = wide.index.astype("datetime64[ns, UTC]")
    return wide[sorted(wide.columns)].astype(float)


def read_gridemissions_output(path: str) -> pd.DataFrame:
    """Reads a gridemissions output CSV, with columns in alphabetical order.

    Uses the round-trip float parser, since pandas' default parser can be off by one
    unit in the last place for some values.

    Args:
        path (str): path to the CSV file.

    Returns:
        pd.DataFrame: one column per series, indexed by "datetime_utc".
    """
    wide = pd.read_csv(
        path, index_col=0, parse_dates=True, float_precision="round_trip"
    )
    wide.index = wide.index.astype("datetime64[ns, UTC]")
    wide.index.name = "datetime_utc"
    return wide[sorted(wide.columns)].astype(float)


@pytest.mark.parametrize(
    "fixture_name, output_stages",
    [
        (
            "eia930_raw_2024-11-01_2024-11-14.csv.gz",
            OUTPUT_STAGES,
        ),
        (
            "eia930_raw_2018-06-24_2018-07-07.csv.gz",
            OUTPUT_STAGES | PRE_2018_07_OUTPUT_STAGES,
        ),
    ],
)
def test_reconcile_eia930_matches_gridemissions(
    fixture_name, output_stages, gridemissions_make_dataset, tmp_path
):
    """The data at every step of the reconciliation should match the output of the
    original `gridemissions` package.

    `oge.physics_reconciliation` is a port of the `gridemissions` cleaning code to
    PUDL-shaped tables, and must not change its results beyond the precision of
    gridemissions' solver. Each fixture is two weeks of
    real EIA-930 data for every BA, in the wide format `gridemissions` reads. Every BA
    is included because the optimization couples BAs through interchange, so dropping
    a BA would change the results for the others.

    - The 2024 fixture covers the first weeks of storage reporting in several BAs
      (e.g. ISNE, ERCO, NEVP), when negative storage values started to appear.
    - The 2018 fixture crosses 2018-07-01, so it exercises both the optimization
      without generation by energy source (before) and with it (after).

    For each step, the test checks that both versions have the same series and hours,
    and:
        - exactly the same values after basic and rolling window cleaning, and for the
          weights
        - values within `RECONCILED_TOLERANCE` after reconciliation, and objective
          values within `OBJECTIVE_RELATIVE_TOLERANCE`, since OGE uses a more precise
          solver
    """
    raw_csv = read_raw_fixture(fixture_name)
    # both versions parse the raw fixture with pandas' default float parser, as
    # gridemissions does

    gridemissions_folder = tmp_path / "gridemissions"
    os.makedirs(gridemissions_folder)
    (gridemissions_folder / "eia930_raw.csv").write_text(raw_csv)
    gridemissions_make_dataset(
        start="",
        end="",
        file_name="eia930",
        tmp_folder=str(gridemissions_folder),
        folder_hist=str(gridemissions_folder),
        scrape=False,
        add_ca_fuels=False,
        calc_consumed=False,
    )

    oge_folder = str(tmp_path / "oge")
    raw = pd.read_csv(io.StringIO(raw_csv), index_col=0, parse_dates=True)
    reconcile_eia930(convert_wide_to_tables(raw), oge_folder)

    for file_name, stage in output_stages.items():
        expected = read_gridemissions_output(gridemissions_folder / file_name)
        result = convert_tables_to_wide(read_eia930_tables(oge_folder, stage))
        if stage in RECONCILED_STAGES:
            # same series and hours, and values within the solver's precision
            pd.testing.assert_frame_equal(
                result,
                expected,
                check_exact=False,
                rtol=0,
                atol=RECONCILED_TOLERANCE,
                check_names=False,
                obj=stage,
            )
        else:
            # same series, same hours and exactly the same values
            pd.testing.assert_frame_equal(
                result, expected, check_exact=True, check_names=False, obj=stage
            )

    # weights are only written for every series in the rolling step
    expected = read_gridemissions_output(gridemissions_folder / "eia930_weights.csv")
    result = convert_tables_to_wide(
        read_eia930_tables(oge_folder, "weights"), value_suffix="_weight"
    )
    pd.testing.assert_frame_equal(
        result, expected, check_exact=True, check_names=False, obj="weights"
    )

    # the objective value of each hour
    objective_files = {"eia930_objective.csv": "eia930_objective.parquet"}
    if "eia930_opt_no_src.csv" in output_stages:
        objective_files["eia930_objective_no_src.csv"] = (
            "eia930_objective_without_energy_sources.parquet"
        )
    for gridemissions_file, oge_file in objective_files.items():
        expected = pd.read_csv(
            gridemissions_folder / gridemissions_file,
            index_col=0,
            parse_dates=True,
            float_precision="round_trip",
        )["CleaningObjective"]
        result = pd.read_parquet(os.path.join(oge_folder, oge_file))["objective"]
        np.testing.assert_allclose(
            result.to_numpy(),
            expected.to_numpy(),
            rtol=OBJECTIVE_RELATIVE_TOLERANCE,
            err_msg=oge_file,
        )


@pytest.mark.parametrize(
    "fixture_name",
    [
        "eia930_raw_2024-11-01_2024-11-14.csv.gz",
        "eia930_raw_2018-06-24_2018-07-07.csv.gz",
    ],
)
def test_reconcile_energy_balance_solves_the_problem(fixture_name):
    """The reconciled data should satisfy every constraint and be the exact solution.

    This protects the formulation in `HourlyReconciliationProblem`, independently of
    `gridemissions`. For each fixture, after basic and rolling window cleaning, the
    hours are reconciled as in `reconcile_eia930`: without generation by energy source
    before 2018-07-01 (only in the 2018 fixture) and with it from that date, and:
        - every constraint of the problem holds to `CONSTRAINT_TOLERANCE`
        - the reconciled data are within `SOLUTION_TOLERANCE` of a solve with much
          tighter tolerances, which shows the default solver settings are precise
          enough
    """
    series = pivot_eia930_tables_to_series(
        convert_wide_to_tables(
            pd.read_csv(
                io.StringIO(read_raw_fixture(fixture_name)),
                index_col=0,
                parse_dates=True,
            )
        )
    )
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
        check_reconciliation_constraints(reconciled, include_energy_sources)

        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(
                physics_reconciliation,
                "CLARABEL_SETTINGS",
                {"tol_gap_abs": 1e-13, "tol_gap_rel": 1e-13, "tol_feas": 1e-13},
            )
            reference, _ = reconcile_energy_balance(
                series.loc[hours], weights.loc[hours], include_energy_sources
            )
        # every reconciled value is close to the exact solution
        pd.testing.assert_frame_equal(
            reconciled, reference, check_exact=False, rtol=0, atol=SOLUTION_TOLERANCE
        )


def check_reconciliation_constraints(
    reconciled: pd.DataFrame, include_energy_sources: bool
):
    """Asserts that reconciled data satisfy every constraint of the reconciliation.

    The constraints are checked directly on the series, independently of how
    `HourlyReconciliationProblem` builds them.

    Args:
        reconciled (pd.DataFrame): reconciled series data frame.
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
    # demand and net generation are at least the lower bound
    lower_bound = OPTIMIZATION_LOWER_BOUND - CONSTRAINT_TOLERANCE
    assert (demand >= lower_bound).all(axis=None)
    assert (net_generation >= lower_bound).all(axis=None)

    is_generation = fields == NET_GENERATION_BY_ENERGY_SOURCE
    # generation by energy source is only in the problem if it is reconciled
    assert is_generation.any() == include_energy_sources
    if include_energy_sources:
        # net generation == sum of generation by energy source
        generation_sum = by_ba(NET_GENERATION_BY_ENERGY_SOURCE)
        assert ((net_generation - generation_sum).abs() < CONSTRAINT_TOLERANCE).all(
            axis=None
        )
        # generation by energy source is at least the lower bound
        assert (reconciled.loc[:, is_generation] >= lower_bound).all(axis=None)


def test_convert_wide_to_tables_round_trips():
    """The test's own conversion helpers should round trip the raw fixture exactly.

    Otherwise a difference in the parity test could come from the helpers rather than
    from `oge.physics_reconciliation`.
    """
    raw = pd.read_csv(
        io.StringIO(read_raw_fixture("eia930_raw_2024-11-01_2024-11-14.csv.gz")),
        index_col=0,
        parse_dates=True,
    )
    raw.index = raw.index.astype("datetime64[ns, UTC]")
    raw.index.name = "datetime_utc"
    expected = raw[sorted(raw.columns)].astype(float)
    result = convert_tables_to_wide(convert_wide_to_tables(raw))
    # every column name in the fixture is parsed
    assert all(
        re.match(r"^EBA\.[^-]+-[^.]+\.[A-Z]+(\.[A-Z]+)?\.H$", c) for c in raw.columns
    )
    pd.testing.assert_frame_equal(result, expected, check_exact=True)


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
