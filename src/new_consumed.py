import numpy as np
import pandas as pd
import os
import sys

from gridemissions.load import BaData

from src.output_data import (
    GENERATED_EMISSION_RATE_COLS,
    output_to_results,
    TIME_RESOLUTIONS,
)

# Regions outside the US that are import-only.
IMPORT_REGIONS = [
    "IESO",
    "HQT",
    "AESO",
    "NBSO",
    "BCHA",
    "MHEB",
    # "CFE",
]

# Defined in output_data, written to each BA file
EMISSION_COLS = [
    "co2_mass_lb_for_electricity",
    "ch4_mass_lb_for_electricity",
    "n2o_mass_lb_for_electricity",
    "co2e_mass_lb_for_electricity",
    "nox_mass_lb_for_electricity",
    "so2_mass_lb_for_electricity",
    "co2_mass_lb_for_electricity_adjusted",
    "ch4_mass_lb_for_electricity_adjusted",
    "n2o_mass_lb_for_electricity_adjusted",
    "co2e_mass_lb_for_electricity_adjusted",
    "nox_mass_lb_for_electricity_adjusted",
    "so2_mass_lb_for_electricity_adjusted",
]

CONSUMED_EMISSION_RATE_COLS = [
    "consumed_co2_rate_lb_per_mwh_for_electricity",
    "consumed_ch4_rate_lb_per_mwh_for_electricity",
    "consumed_n2o_rate_lb_per_mwh_for_electricity",
    "consumed_co2e_rate_lb_per_mwh_for_electricity",
    "consumed_nox_rate_lb_per_mwh_for_electricity",
    "consumed_so2_rate_lb_per_mwh_for_electricity",
    "consumed_co2_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_ch4_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_n2o_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_co2e_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_nox_rate_lb_per_mwh_for_electricity_adjusted",
    "consumed_so2_rate_lb_per_mwh_for_electricity_adjusted",
]

FUEL_TYPE_MAP = {
    "COL": "coal",
    "NG": "natural_gas",
    "NUC": "nuclear",
    "OIL": "petroleum",
    "OTH": "total",
    "SUN": "solar",
    "UNK": "total",
    "WAT": "hydro",
    "WND": "wind",
    "GEO": "total",  # TODO double-check that geo is unused
    "BIO": "biomass",
}

POLLS = ["CO2", "CH4", "N2O", "CO2E", "NOX", "SO2"]

ADJUSTMENTS = ["for_electricity", "for_electricity_adjusted"]

SRC = ["COL", "NG", "NUC", "OIL", "OTH", "SUN", "UNK", "WAT", "WND", "GEO", "BIO"]


def get_column(poll: str, adjustment: str, ba: str = ""):
    """
    Return output file column name for a pollutant and adjustment type
    Returns mass columns, not rate columns
    """
    assert poll in POLLS
    assert adjustment in ADJUSTMENTS
    if ba == "":  # no BA, looking for output file column
        column = poll.lower() + "_mass_lb_" + adjustment
        assert column in EMISSION_COLS
    else:
        column = ba + "_" + poll.lower() + "_mass_lb_" + adjustment
    return column


def get_rate_column(poll: str, adjustment: str, generated: bool = True, ba: str = ""):
    """
    Return either generated or consumed output file rate column
    for pollutant `poll` and adjustment `adjustment`
    """
    assert poll in POLLS
    assert adjustment in ADJUSTMENTS
    if generated:
        column = "generated_" + poll.lower() + "_rate_lb_per_mwh_" + adjustment
        assert column in GENERATED_EMISSION_RATE_COLS
    else:
        column = "consumed_" + poll.lower() + "_rate_lb_per_mwh_" + adjustment
        assert column in CONSUMED_EMISSION_RATE_COLS
    if ba != "":  # For internal column use, add ba
        column = ba + "_" + column
    return column


def get_average_emission_factors(prefix: str = "2020/", year: int = 2020):
    """
    Locate per-fuel, per-adjustment, per-poll emission factors.
    Used to fill in emissions from BAs outside of US, where we have generation by
    fuel (from gridemissions) but no hourly-egrid data

    Structure: EMISSIONS_FACTORS[poll][adjustment][fuel]
    """
    genavg = pd.read_csv(
        f"../data/outputs/{prefix}annual_generation_averages_by_fuel_{year}.csv",
        index_col="fuel_category",
    )
    efs = {}
    for pol in POLLS:
        efs[pol] = {}
        for adjustment in ADJUSTMENTS:
            efs[pol][adjustment] = {}
            for fuel in SRC:
                column = get_rate_column(pol, adjustment, generated=True)
                if FUEL_TYPE_MAP[fuel] not in genavg.index:
                    print(
                        f"WARNING: fuel {FUEL_TYPE_MAP[fuel]} not found in file annual_generation_averages_by_fuel_{year}.csv, using average"
                    )
                    efs[pol][adjustment][fuel] = genavg.loc["total", column]
                else:
                    efs[pol][adjustment][fuel] = genavg.loc[FUEL_TYPE_MAP[fuel], column]
    return efs


def consumption_emissions(F, P, ID):
    """
    FROM GRIDEMISSIONS: https://github.com/jdechalendar/gridemissions

    Form and solve linear system to compute consumption emissions

    Parameters
    ----------
    F: np.array
        emissions
    P: np.array
        production
    ID: np.array
        exchanges

    Notes
    -----
    Create linear system to calculate consumption emissions
    - Create import matrix
    - Create linear system and solve:
    f_i^c*(d_i+sum_j t_{ji}) - sum_j t_{ij}*f_j^c = F_i^p
    where:
        f_i^c: consumption emissions at node i
        d_i: demand at node i
        t: trade matrix - t_{ij} is from node i to j
        F_i^p: emissions produced at node i
    Note: np version must be high enough, otherwise np.linalg.cond fails
    on a matrix with only zeros.
    """
    from distutils.version import LooseVersion

    assert LooseVersion(np.__version__) >= LooseVersion("1.15.1")

    # Create and solve linear system
    Imp = (-ID).clip(min=0)  # trade matrix reports exports - we want imports
    I_tot = Imp.sum(axis=1)  # sum over columns
    A = np.diag(P + I_tot) - Imp
    b = F

    perturbed = []
    if np.linalg.cond(A) > (1.0 / sys.float_info.epsilon):
        # matrix is ill-conditioned
        for i in range(len(A)):
            if (np.abs(A[:, i]).sum() == 0.0) & (np.abs(A[i, :]).sum() == 0.0):
                A[i, i] = 1.0  # slightly perturb that element
                perturbed += [i]
                # force this to be zero so the linear system makes sense
                b[i] = 0.0

    X = np.linalg.solve(A, b)

    for j in perturbed:
        if X[j] != 0.0:
            print(b[j])
            print(np.abs(A[j, :]).sum())
            print(np.abs(A[:, j]).sum())
            raise ValueError("X[%d] is %.2f instead of 0" % (j, X[j]))

    return X, len(perturbed)
