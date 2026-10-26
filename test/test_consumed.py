import sys

import numpy as np
import pandas as pd

import oge.consumed as consumed

import pytest

REGIONS = ["A", "B", "C", "D", "E"]
HOURS = pd.date_range("2025-01-01", periods=6, freq="h", tz="UTC", name="datetime_utc")


@pytest.fixture
def hourly_consumed() -> consumed.HourlyConsumed:
    """A small `HourlyConsumed` with synthetic data, built without reading any files.

    Regions A, B and C are U.S. regions, D is import-only (generation from EIA-930) and
    E never generates or trades, so its row and column of the system matrix are zero.
    Hour 1 has a missing interchange value, in hour 2 region C has zero generation but
    exports, and hour 5 has no data except one inconsistent interchange value, so its
    linear system is singular and can't be solved.
    """
    rng = np.random.default_rng(0)
    hc = object.__new__(consumed.HourlyConsumed)
    hc.regions = REGIONS
    hc.import_regions = ["D"]

    # generation for U.S. regions, from the power sector data
    generation = pd.DataFrame(
        rng.uniform(100, 1000, size=(len(HOURS), 4)),
        index=HOURS,
        columns=["A", "B", "C", "E"],
    )
    generation["E"] = 0.0
    generation.loc[HOURS[2], "C"] = 0.0
    generation.loc[HOURS[5]] = 0.0
    hc.generation = generation

    # generation for import-only regions, from EIA-930
    net_generation = pd.DataFrame(
        rng.uniform(100, 1000, size=(len(HOURS), 1)), index=HOURS, columns=["D"]
    )
    net_generation.loc[HOURS[5]] = 0.0
    hc.net_generation = net_generation

    # interchange in both directions for each adjacent pair (positive is exports),
    # plus a pair with a region that isn't part of the linear system
    interchange = {}
    for ri, rj in [("A", "B"), ("B", "C"), ("C", "A"), ("D", "A"), ("A", "Z")]:
        flow = rng.uniform(-100, 100, size=len(HOURS))
        interchange[(ri, rj)] = flow
        interchange[(rj, ri)] = -flow
    interchange = pd.DataFrame(interchange, index=HOURS)
    interchange.columns = pd.MultiIndex.from_tuples(interchange.columns)
    interchange.loc[HOURS[1], ("A", "B")] = np.nan
    interchange.loc[HOURS[2], ("C", "A")] = 50.0
    interchange.loc[HOURS[2], ("A", "C")] = -50.0
    interchange.loc[HOURS[2], ("B", "C")] = 20.0
    interchange.loc[HOURS[2], ("C", "B")] = -20.0
    # in the hour with no data, A reports importing from C but C reports nothing, so
    # C's row of the system matrix is zero but its column isn't
    interchange.loc[HOURS[5]] = 0.0
    interchange.loc[HOURS[5], ("A", "C")] = -30.0
    hc.interchange = interchange

    # mass emissions, with a different emission rate per region and series
    all_generation = pd.concat([generation, net_generation], axis=1)[REGIONS]
    hc.pollutant_mass_emissions = {
        (adj, pol): all_generation * rng.uniform(0, 2000, size=len(REGIONS))
        for pol, adj in consumed.POLLUTANT_ADJUSTMENTS
    }

    hc.results = hc._build_results()
    return hc


def reference_consumption_emissions(F, P, ID):
    """Per-series consumption_emissions from before it accepted several series."""
    Imp = (-ID).clip(min=0)
    I_tot = Imp.sum(axis=1)
    A = np.diag(P + I_tot) - Imp
    b = F
    if np.linalg.cond(A) > (1.0 / sys.float_info.epsilon):
        for i in range(len(A)):
            if (np.abs(A[:, i]).sum() == 0.0) & (np.abs(A[i, :]).sum() == 0.0):
                A[i, i] = 1.0
                b[i] = 0.0
    return np.linalg.solve(A, b)


def reference_build_matrices(hc, pol, adj, date):
    """Per-hour, per-series build_matrices from before it was vectorized."""
    ID = np.zeros((len(hc.regions), len(hc.regions)))
    for i, ri in enumerate(hc.regions):
        for j, rj in enumerate(hc.regions):
            if (ri, rj) in hc.interchange.columns:
                ID[i][j] = hc.interchange.loc[date, (ri, rj)]
    E = hc.pollutant_mass_emissions[(adj, pol)].loc[date, hc.regions].to_numpy()
    G = np.zeros(len(hc.regions))
    for i, r in enumerate(hc.regions):
        if r in hc.import_regions:
            G[i] = hc.net_generation.loc[date, r]
        else:
            G[i] = hc.generation.loc[date, r]
    E = np.nan_to_num(E)
    G = np.nan_to_num(G)
    ID = np.nan_to_num(ID)
    to_fix = (ID.sum(axis=1) > 0) & (G == 0)
    ID[:, to_fix] = 0
    ID[to_fix, :] = 0
    return E, G, ID


def reference_run(hc):
    """Per-hour, per-series run from before it was vectorized.

    Returns:
        dict[str, pd.DataFrame]: consumed emission rates by region.
    """
    results = {r: pd.DataFrame(index=hc.generation.index) for r in hc.regions}
    for pol, adj in consumed.POLLUTANT_ADJUSTMENTS:
        col = consumed.get_pollutant_rate_column_name(pol, adj, generated=False)
        for date in hc.generation.index:
            E, G, ID = reference_build_matrices(hc, pol, adj, date)
            try:
                X = reference_consumption_emissions(E, G, ID)
            except np.linalg.LinAlgError:
                X = np.full(len(hc.regions), np.nan)
            for i, r in enumerate(hc.regions):
                results[r].loc[date, col] = X[i]
    return results


def test_build_matrices_matches_per_hour(hourly_consumed):
    """The vectorized arrays equal the per-hour, per-series matrices."""
    E, G, ID = hourly_consumed.build_matrices()
    for t, date in enumerate(HOURS):
        for s, (pol, adj) in enumerate(consumed.POLLUTANT_ADJUSTMENTS):
            ref_E, ref_G, ref_ID = reference_build_matrices(
                hourly_consumed, pol, adj, date
            )
            np.testing.assert_array_equal(E[t, :, s], ref_E)
            np.testing.assert_array_equal(G[t], ref_G)
            np.testing.assert_array_equal(ID[t], ref_ID)

    # region C's row and column are dropped in the hour it exports without generating
    assert (ID[2, 2, :] == 0).all() and (ID[2, :, 2] == 0).all()
    assert (ID[3, 2, :] != 0).any()


def test_run_matches_per_hour(hourly_consumed):
    """Consumed rates equal the per-hour, per-series solves, NaN in the same hours."""
    expected = reference_run(hourly_consumed)
    hourly_consumed.run()
    for r in REGIONS:
        for col in expected[r].columns:
            np.testing.assert_allclose(
                hourly_consumed.results[r][col].to_numpy(),
                expected[r][col].to_numpy(),
                rtol=1e-12,
                atol=0,
                err_msg=f"{r} {col}",
            )

    # only the hour with no data fails to solve, and the isolated region is zero
    rates = pd.concat(hourly_consumed.results, axis=1).filter(regex="consumed_.*_rate_")
    assert rates.isna().all(axis=1).tolist() == [False] * 5 + [True]
    assert (rates.loc[HOURS[:5], "E"] == 0).all().all()


def test_consumption_emissions_with_several_series():
    """Solving several series at once matches solving them one at a time."""
    rng = np.random.default_rng(1)
    P = np.array([500.0, 300.0, 0.0, 200.0])
    ID = np.array(
        [
            [0.0, 40.0, 0.0, -30.0],
            [-40.0, 0.0, 0.0, 10.0],
            [0.0, 0.0, 0.0, 0.0],
            [30.0, -10.0, 0.0, 0.0],
        ]
    )
    F = rng.uniform(0, 1e6, size=(4, 3))
    F[2] = 5.0  # an isolated region's emissions are forced to zero
    F_before = F.copy()

    X, n_isolated = consumed.consumption_emissions(F, P, ID)

    assert n_isolated == 1
    np.testing.assert_array_equal(F, F_before)
    for s in range(F.shape[1]):
        np.testing.assert_allclose(
            X[:, s],
            reference_consumption_emissions(F[:, s].copy(), P, ID.copy()),
            rtol=1e-12,
        )
    assert (X[2] == 0).all()
