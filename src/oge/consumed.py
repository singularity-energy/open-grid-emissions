import numpy as np
import pandas as pd
import os
import sys

from packaging.version import Version

from oge.physics_reconciliation import ENERGY_SOURCES, read_eia930_tables
from oge.filepaths import reference_table_folder, results_folder
from oge.logging_util import get_logger
from oge.constants import TIME_RESOLUTIONS
from oge.output_data import (
    GENERATED_EMISSION_RATE_COLS,
    CONSUMED_EMISSION_RATE_COLS,
    output_to_results,
)

logger = get_logger(__name__)


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

# maps EIA-930 energy sources to the fuel categories of the US fleet average emission
# rates used for BAs outside of the US. Energy storage uses the rate of the US storage
# fleet, which has no direct emissions, so its rate is zero. Emissions from the
# electricity used to charge storage are not yet tracked through storage.
FUEL_TYPE_MAP = {
    "coal": "coal",
    "gas": "natural_gas",
    "nuclear": "nuclear",
    "oil": "petroleum",
    "other": "total",
    "solar": "solar",
    "unknown": "total",
    "hydro": "hydro",
    "wind": "wind",
    "pumped_storage": "storage",
    "solar_w_integrated_battery_storage": "solar",
    "wind_w_integrated_battery_storage": "wind",
    "battery_storage": "storage",
    "unknown_energy_storage": "storage",
    "other_energy_storage": "storage",
    "geothermal": "geothermal",
}

POLLUTANTS = ["CO2", "CH4", "N2O", "CO2E", "NOX", "SO2"]

ADJUSTMENTS = ["for_electricity", "for_electricity_adjusted"]

# every (pollutant, adjustment) series solved for consumed emission rates, in the
# order of the last axis of the emission array built by `HourlyConsumed.build_matrices`
POLLUTANT_ADJUSTMENTS = [(pol, adj) for pol in POLLUTANTS for adj in ADJUSTMENTS]


def get_pollutant_mass_column_name(poll: str, adjustment: str, ba: str = "") -> str:
    """Get name of pollutant mass column.

    Args:
        poll (str): pollutant name (e.g., "CO2")
        adjustment (str): adjustment type (e.g., "for_electricity")
        ba (str, optional): Balancing Authority code. Defaults to "".

    Returns:
        str: column name for the pollutant mass
    """
    assert poll in POLLUTANTS
    assert adjustment in ADJUSTMENTS
    if ba == "":  # no BA, looking for output file column
        column = poll.lower() + "_mass_lb_" + adjustment
        assert column in EMISSION_COLS
    else:
        column = ba + "_" + poll.lower() + "_mass_lb_" + adjustment
    return column


def get_pollutant_rate_column_name(
    poll: str, adjustment: str, generated: bool = True, ba: str = ""
) -> str:
    """Get name of pollutant rate column.

    Args:
        poll (str): pollutant name (e.g., "CO2")
        adjustment (str): adjustment type (e.g., "for_electricity")
        generated (bool, optional): whether the rate is for generated emissions. If
            False, the rate is for consumed emissions. Defaults to True.
        ba (str, optional): Balancing Authority code. Defaults to "".

    Returns:
        str: column name for the pollutant rate
    """
    assert poll in POLLUTANTS
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


def get_average_emission_factors(prefix: str) -> dict[str, dict[str, dict[str, float]]]:
    """Get U.S. average emission factors by-pollutant, by-adjustment and by-fuel type.

    Used to fill in emissions from BAs outside of US, where we have generation by
    fuel (from gridemissions) but no open-grid-emissions data

    We use `gridemissions` assumptions for fuel mix for non-US BAs, which are simple
    and not time-varying

    Args:
        prefix (str): path prefix to results folder

    Returns:
        dict[str, dict[str, dict[str, float]]]: emission factors.
    """
    genavg = pd.read_csv(
        results_folder(f"{prefix}/power_sector_data/annual/us_units/US.csv"),
        index_col="fuel_category",
    )
    efs = {}
    for pol in POLLUTANTS:
        efs[pol] = {}
        for adjustment in ADJUSTMENTS:
            efs[pol][adjustment] = {}
            for fuel in ENERGY_SOURCES:
                column = get_pollutant_rate_column_name(pol, adjustment, generated=True)
                # guard: storage has no direct emissions, so use zero rather than the
                # total average if the US fleet data has no storage row
                if FUEL_TYPE_MAP[fuel] == "storage" and "storage" not in genavg.index:
                    efs[pol][adjustment][fuel] = 0.0
                elif FUEL_TYPE_MAP[fuel] not in genavg.index:
                    logger.warning(
                        f"fuel {FUEL_TYPE_MAP[fuel]} not found in US fleet average "
                        "data, using total average"
                    )
                    efs[pol][adjustment][fuel] = genavg.loc["total", column]
                else:
                    efs[pol][adjustment][fuel] = genavg.loc[FUEL_TYPE_MAP[fuel], column]
    return efs


def consumption_emissions(
    F: np.ndarray, P: np.ndarray, ID: np.ndarray
) -> tuple[np.ndarray, int]:
    """Form and solve linear system to compute consumption emissions

    From GRIDEMISSIONS: https://github.com/jdechalendar/gridemissions

    The system matrix only depends on generation and interchange, so several emission
    series for the same hour can be solved at once by passing them as the columns of
    `F`.

    Args:
        F (np.ndarray): emission vector (n,), or emission matrix (n, k) with one column
            per emission series.
        P (np.ndarray): generation vector (n,).
        ID (np.ndarray): interchange matrix (n, n).

    Returns:
        tuple[np.ndarray, int]: intensity consumption emissions, with the same shape as
            `F`, and number of perturbed nodes

    Notes:
    Create linear system to calculate consumption emissions
    - See https://www.pnas.org/doi/full/10.1073/pnas.1912950116 (Equations 1 to 4)
    - Start from Equation 1: x_i * d_i = f_i + sum_j x_j u_{ij}) - sum_k x_i u_{ki}
    where:
        x_i: intensity of electricity consumed at node i
        x_j: intensity of electricity consumed at node j
        d_i: electricity consumed at node i
        u_{ij}: electricity imported from j to i
        u_{ki}: electricity exported from i to k
        f_i: pollutant production at node i
    """

    # numpy version must be high enough, otherwise np.linalg.cond fails on a matrix with
    # only zeros.
    assert Version(np.__version__) >= Version("1.15.1")

    # Create and solve linear system
    Imp = (-ID).clip(min=0)  # trade matrix reports exports - we want imports
    I_tot = Imp.sum(axis=1)  # sum over columns
    A = np.diag(P + I_tot) - Imp
    # copy so that zeroing perturbed nodes below doesn't modify the caller's array
    b = np.array(F, dtype=np.float64)

    perturbed = []
    if np.linalg.cond(A) > (1.0 / sys.float_info.epsilon):
        # matrix is ill-conditioned
        for i in range(len(A)):
            if (np.abs(A[:, i]).sum() == 0.0) & (np.abs(A[i, :]).sum() == 0.0):
                A[i, i] = 1.0  # slightly perturb that element
                perturbed += [i]
                # force this to be zero (in every emission series) so the linear
                # system makes sense
                b[i] = 0.0

    X = np.linalg.solve(A, b)

    for j in perturbed:
        if np.any(X[j] != 0.0):
            logger.warning(f"b[{j}]: {b[j]}")
            logger.warning(f"sum of abs(A[{j}, :]): {np.abs(A[j, :]).sum()}")
            logger.warning(f"sum of abs(A[:, {j}]): {np.abs(A[:, j]).sum()}")
            raise ValueError(f"X[{j}] is {X[j]} instead of 0")

    return X, len(perturbed)


class HourlyConsumed:
    """Class to load data, calculate consumed rates, and write results to files."""

    def __init__(
        self,
        eia930_folder: str,
        prefix: str,
        year: int,
        skip_outputs: bool = False,
    ):
        self.prefix = prefix
        self.year = year
        self.skip_outputs = skip_outputs

        # reconciled 930 data
        self._load_eia930_data(eia930_folder)

        # Emission factors for non-US bas
        self.default_factors = get_average_emission_factors(prefix)

        # Look up lists of BAs with specific requirements
        self.import_regions, self.generation_regions = self._get_special_regions()

        # Load generated pollutant mass emissions and generation for import-only BAs
        self.pollutant_mass_emissions, self.generation = (
            self._load_pollutant_mass_emissions_and_generation_for_import_only_regions()
        )

        # Identify shared BAs
        regions = set(self.eia930_regions)
        regions = regions.intersection(set(self.generation.columns))
        # Add back import-only regions
        regions = regions.union(set(self.import_regions))
        self.regions = list(regions)

        # Build result df
        self.results = self._build_results()

    def _load_eia930_data(self, eia930_folder: str):
        """Loads reconciled EIA-930 data into one wide table per series type.

        Sets `demand` and `net_generation` (one column per BA), `interchange` (one
        column per (BA, adjacent BA)) and `net_generation_by_energy_source` (one column
        per (BA, energy source)), all indexed by "datetime_utc", and `eia930_regions`,
        the BAs with any EIA-930 data.

        Args:
            eia930_folder (str): folder with the reconciled EIA-930 parquet files.
        """
        tables = read_eia930_tables(eia930_folder, "reconciled")

        def pivot(table, columns, values):
            wide = table.pivot(index="datetime_utc", columns=columns, values=values)
            wide.index = wide.index.astype("datetime64[ns, UTC]")
            return wide

        self.demand = pivot(
            tables.operations, "balancing_authority_code_eia", "demand_mwh"
        )
        self.net_generation = pivot(
            tables.operations, "balancing_authority_code_eia", "net_generation_mwh"
        )
        self.interchange = pivot(
            tables.interchange,
            ["balancing_authority_code_eia", "balancing_authority_code_adjacent_eia"],
            "interchange_mwh",
        )
        self.net_generation_by_energy_source = pivot(
            tables.net_generation_by_energy_source,
            ["balancing_authority_code_eia", "generation_energy_source"],
            "net_generation_mwh",
        )
        self.eia930_regions = sorted(
            set(tables.operations["balancing_authority_code_eia"])
            | set(tables.interchange["balancing_authority_code_eia"])
            | set(
                tables.net_generation_by_energy_source["balancing_authority_code_eia"]
            )
        )

    def _get_special_regions(self) -> tuple[list[str], list[str]]:
        """Get import-only and generation-only regions.

        Returns:
            tuple[list[str], list[str]]: first element is list of import-only regions,
                second element is list of generation-only regions
        """

        self.ba_ref = pd.read_csv(
            reference_table_folder("ba_reference.csv"), index_col="ba_code"
        )

        # Get generation-only regions
        generation_only = list(
            self.ba_ref[self.ba_ref["ba_category"] == "generation_only"].index
        )

        # Get regions that interchange with U.S. regions but whose generation and
        # emissions need to be derived from EIA-930
        import_only = [
            b
            for b in self.ba_ref[self.ba_ref["us_ba"] == "No"].index
            if b in self.eia930_regions
        ]
        return import_only, generation_only

    def _build_results(self) -> dict[str, pd.DataFrame]:
        """Builds result data frame per output file.

        Returns:
            dict[str, pd.DataFrame]: mapping of BA code to result dataframe.
        """
        results = {}
        cols = []
        for pol in POLLUTANTS:
            for adj in ADJUSTMENTS:
                cols.append(
                    get_pollutant_rate_column_name(pol, adjustment=adj, generated=False)
                )
                cols.append(get_pollutant_mass_column_name(pol, adjustment=adj))
        cols.append("net_consumed_mwh")
        for ba in self.regions:
            results[ba] = pd.DataFrame(
                index=self.generation.index, columns=cols, dtype=np.float64
            )
        return results

    def output_results(self):
        """Write results to files.

        In order to calculate the monthly and annual consumed emission rates correctly,
        we need to first calculate the hourly consumed electricity and the hourly
        consumed pollutant mass emissions.

        * Consumed electricity is calculated from EIA-930 demand
        * Consumed pollutant mass emissions is calculated as:
              consumed electricity * consumed pollutant emission rate
          where the consumed pollutant emission rate is solved in the linear system of
          equations.
        """
        for ba in self.regions:
            if (ba in self.import_regions) or (ba in self.generation_regions):
                continue

            # Add hourly consumed electricity taken as hourly EIA-930 demand to data
            # frame. This will be used to calculate monthly and annual consumed
            # emission rates from hourly ones.
            # Will be dropped from final output files
            self.results[ba]["net_consumed_mwh"] = self.demand[ba][
                self.generation.index
            ]

            # Add hourly mass emissions taken as hourly consumed
            # electricity multiplied by consumed emission rate to data frame.
            # Will be dropped from final output files
            for pol in POLLUTANTS:
                for adj in ADJUSTMENTS:
                    self.results[ba][
                        get_pollutant_mass_column_name(pol, adjustment=adj)
                    ] = (
                        self.results[ba][
                            get_pollutant_rate_column_name(
                                pol, adjustment=adj, generated=False
                            )
                        ]
                        * self.results[ba]["net_consumed_mwh"]
                    )

            # Although we directly calculate rates, to calculate annual average rates
            # we sum emissions and generation then divide.
            for time_resolution in TIME_RESOLUTIONS:
                time_dat = self.results[ba].copy(deep=True)

                # Get local timezone
                assert not pd.isnull(self.ba_ref.loc[ba, "timezone_local"])
                time_dat["datetime_local"] = time_dat.index.tz_convert(
                    self.ba_ref.loc[ba, "timezone_local"]
                )
                time_dat = time_dat.reset_index()  # move datetime_utc to column
                time_dat = time_dat[
                    time_dat["datetime_local"].dt.year == self.year
                ]  # keep year of local data

                if time_resolution == "hourly":
                    # No resampling needed; keep timestamp cols in output
                    time_cols = ["datetime_utc", "datetime_local"]
                    missing_hours = time_dat[time_dat.isna().any(axis=1)]
                    if len(missing_hours) > 0:
                        logger.warning(
                            f"{len(missing_hours)} hours are missing in {ba} consumed data"
                        )
                elif time_resolution == "monthly":
                    time_dat["month"] = time_dat["datetime_local"].dt.month
                    # Aggregate to appropriate resolution
                    time_dat = (
                        time_dat.groupby("month")[EMISSION_COLS + ["net_consumed_mwh"]]
                        .sum()
                        .reset_index()  # move "month" to column
                    )
                    time_cols = ["month"]
                elif time_resolution == "annual":
                    time_dat["year"] = time_dat["datetime_local"].dt.year
                    # Aggregate to appropriate resolution
                    time_dat = (
                        time_dat.groupby("year")[EMISSION_COLS + ["net_consumed_mwh"]]
                        .sum()
                        .reset_index()  # move "year" to column
                    )
                    time_cols = ["year"]

                # Calculate rates from summed mass emissions and consumed electricity
                # for aggregated time resolutions.
                if time_resolution != "hourly":
                    for pol in POLLUTANTS:
                        for adj in ADJUSTMENTS:
                            rate_col = get_pollutant_rate_column_name(
                                pol, adj, generated=False
                            )
                            emission_col = get_pollutant_mass_column_name(pol, adj)
                            time_dat[rate_col] = (
                                time_dat[emission_col] / time_dat["net_consumed_mwh"]
                            )

                # Output
                output_to_results(
                    time_dat[time_cols + CONSUMED_EMISSION_RATE_COLS],
                    self.year,
                    ba,
                    f"/carbon_accounting/{time_resolution}/",
                    self.prefix,
                    skip_outputs=self.skip_outputs,
                )

    def _impute_border_hours(self, temp: pd.Series) -> pd.Series:
        """Add border hours to time series.

        Add three hours to beginning and end of series. Impute hours by taking same
        hour from previous (or next) day of series. This matches EIA-930's imputation
        approach.

        Args:
            temp (pd.Series): time series (generation or emissions)

        Returns:
            pd.Series: time series with imputed border hours
        """
        temp = temp.dropna()
        last_day = temp.index.max()
        first_day = temp.index.min()
        for hour in range(1, 4):
            new_hour = first_day - pd.DateOffset(hours=hour)
            # Get closest hour
            best = new_hour + pd.DateOffset(days=1)
            # Fill in
            temp[new_hour] = temp[best]

            # Now at end of time series
            new_hour = last_day + pd.DateOffset(hours=hour)
            best = new_hour - pd.DateOffset(days=1)
            temp[new_hour] = temp[best]

        return temp.sort_index()

    def _load_pollutant_mass_emissions_and_generation_for_import_only_regions(
        self,
    ) -> tuple[dict[tuple[str, str], pd.DataFrame], pd.DataFrame]:
        """Get pollutant mass emissions and generation for import-only BAs.

        Returns:
            tuple[dict[tuple[str, str], pd.DataFrame], pd.DataFrame]: first element is a
                dictionary mapping (adjustment, pollutant) to data frame of mass
                emissions by BA; second element is a data frame of generation by BA.
        """
        mass_emissions = {}  # (adj, pol) -> {(BA, mass emission series)}
        ba_generation = {}
        for f in os.listdir(
            results_folder(f"{self.prefix}/power_sector_data/hourly/us_units/")
        ):
            # TODO: delete this message, for testing
            logger.info(f"Loading {f}")
            if ".DS_Store" in f:
                continue
            ba = pd.read_csv(
                results_folder(f"{self.prefix}/power_sector_data/hourly/us_units/") + f,
                index_col="datetime_utc",
                parse_dates=True,
            )
            ba = ba[ba["fuel_category"] == "total"]
            ba_name = f.replace(".csv", "")
            for adj in ADJUSTMENTS:
                for pol in POLLUTANTS:
                    ba_mass_emissions = mass_emissions.get((adj, pol), {})
                    ba_mass_emissions[ba_name] = self._impute_border_hours(
                        ba[get_pollutant_mass_column_name(pol, adjustment=adj)]
                    )
                    mass_emissions[(adj, pol)] = ba_mass_emissions
            ba_generation[ba_name] = self._impute_border_hours(ba["net_generation_mwh"])

        # Make each mass emission into a data frame and add mass emissions for
        # import-only regions
        for pol in POLLUTANTS:
            for adj in ADJUSTMENTS:
                all_ba_mass_emissions = pd.DataFrame(mass_emissions[(adj, pol)])

                # Add import regions to data frame
                for ba in self.import_regions:
                    # sum over energy sources, in the order of ENERGY_SOURCES
                    ba_mass_emissions = pd.Series(0, index=self.demand.index)
                    for src in ENERGY_SOURCES:
                        if (ba, src) in self.net_generation_by_energy_source.columns:
                            ba_mass_emissions = (
                                ba_mass_emissions
                                + self.default_factors[pol][adj][src]
                                * self.net_generation_by_energy_source[(ba, src)]
                            )
                    all_ba_mass_emissions.loc[:, ba] = ba_mass_emissions

                # Cut off emissions at 8 hours after UTC year
                all_ba_mass_emissions = all_ba_mass_emissions[
                    : f"{self.year + 1}-01-01 08:00:00+00:00"
                ]
                mass_emissions[(adj, pol)] = all_ba_mass_emissions

        # Make generation data frame
        generation = pd.DataFrame(data=ba_generation)
        generation = generation[: f"{self.year + 1}-01-01 08:00:00+00:00"]

        return mass_emissions, generation

    def build_matrices(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build emission, generation, and interchange arrays for every hour.

        The arrays cover every hour of `self.generation.index` and every region of
        `self.regions`, in that order, so each hour's linear system is a slice of them.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: emission array (hours, regions,
                series) with one series per entry of `POLLUTANT_ADJUSTMENTS`,
                generation array (hours, regions), and interchange array (hours,
                regions, regions).
        """
        hours = self.generation.index
        position = {r: i for i, r in enumerate(self.regions)}

        # Build transmission matrices from cleaned EIA-930 data. Only pairs of regions
        # that both take part in the linear system and have interchange data are
        # filled; every other pair is zero.
        pairs = [
            (ri, rj)
            for (ri, rj) in self.interchange.columns
            if (ri in position) and (rj in position)
        ]
        ID = np.zeros((len(hours), len(self.regions), len(self.regions)))
        if len(pairs) > 0:
            rows = [position[ri] for ri, _ in pairs]
            cols = [position[rj] for _, rj in pairs]
            ID[:, rows, cols] = self.interchange.loc[hours, pairs].to_numpy()

        # Build pollutant mass arrays, using default per-fuel factors for import-only
        # regions, with one series per (pollutant, adjustment) along the last axis.
        E = np.stack(
            [
                self.pollutant_mass_emissions[(adj, pol)]
                .loc[hours, self.regions]
                .to_numpy()
                for pol, adj in POLLUTANT_ADJUSTMENTS
            ],
            axis=-1,
        )

        # Build generation array, using 930 for import-only regions
        G = np.zeros((len(hours), len(self.regions)))
        for i, r in enumerate(self.regions):
            if r in self.import_regions:
                G[:, i] = self.net_generation.loc[hours, r].to_numpy()
            else:
                G[:, i] = self.generation.loc[hours, r].to_numpy()

        E = np.nan_to_num(E)
        G = np.nan_to_num(G)
        ID = np.nan_to_num(ID)

        # In some cases, we have zero generation but non-zero transmission. For those
        # region-hours, drop the region's row and column of the interchange matrix.
        to_fix = (ID.sum(axis=2) > 0) & (G == 0)
        ID[to_fix[:, :, np.newaxis] | to_fix[:, np.newaxis, :]] = 0.0

        return E, G, ID

    def run(self):
        """Orchestration method to calculate consumed emission rates"""
        logger.info("Solving consumed emissions for all pollutants and adjustments...")
        E, G, ID = self.build_matrices()

        # Solve each hour's linear system once for all emission series, since the
        # system matrix only depends on generation and interchange. Hours that can't
        # be solved stay NaN.
        consumed_emissions = np.full(E.shape, np.nan)
        total_failed = 0
        for t in range(len(G)):
            # If we have NaNs, don't bother
            if np.isnan(G[t]).any() or np.isnan(E[t]).any() or np.isnan(ID[t]).any():
                continue
            try:
                consumed_emissions[t], _ = consumption_emissions(E[t], G[t], ID[t])
            except np.linalg.LinAlgError:
                # These issues happen at boundary hours (beginning and end of year)
                # where we don't have full data for all BAs
                total_failed += 1
        if total_failed > 0:
            logger.warning(
                f"{total_failed} hours failed to solve for consumed emissions, for "
                "every pollutant and adjustment."
            )

        # Export each region's consumed emission rates as whole columns
        for s, (pol, adj) in enumerate(POLLUTANT_ADJUSTMENTS):
            col = get_pollutant_rate_column_name(pol, adjustment=adj, generated=False)
            for i, r in enumerate(self.regions):
                self.results[r][col] = consumed_emissions[:, i, s]
