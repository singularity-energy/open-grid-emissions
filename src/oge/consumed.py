import numpy as np
import pandas as pd
import os
import sys

from packaging.version import Version

from gridemissions.load import BaData
from gridemissions.eia_api import KEYS, SRC

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
    "GEO": "total",
    "BIO": "biomass",
}

POLLUTANTS = ["CO2", "CH4", "N2O", "CO2E", "NOX", "SO2"]

ADJUSTMENTS = ["for_electricity", "for_electricity_adjusted"]


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
            for fuel in SRC:
                column = get_pollutant_rate_column_name(pol, adjustment, generated=True)
                if FUEL_TYPE_MAP[fuel] not in genavg.index:
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

    Args:
        F (np.ndarray): emission vector.
        P (np.ndarray): generation vector.
        ID (np.ndarray): interchange matrix.

    Returns:
        tuple[np.ndarray, int]: intensity consumption emissions vector and number of
            perturbed nodes

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
            logger.warning("\n" + b[j].to_string())
            logger.warning("\n" + np.abs(A[j, :]).sum().to_string())
            logger.warning("\n" + np.abs(A[:, j]).sum().to_string())
            raise ValueError("X[%d] is %.2f instead of 0" % (j, X[j]))

    return X, len(perturbed)


class HourlyConsumed:
    """Class to load data, calculate consumed rates, and write results to files."""

    def __init__(
        self,
        eia930_file: str,
        prefix: str,
        year: int,
        small: bool = False,
        skip_outputs: bool = False,
    ):
        self.prefix = prefix
        self.year = year
        self.small = small
        self.skip_outputs = skip_outputs

        # 930 data
        self.eia930 = BaData(eia930_file)

        # Round 930 data to zero to clear out imputed values.
        # to mimic behavior of eia930.remove_imputed_ones,
        # round all values with abs(x) < 1.5 to 0
        self.eia930.df[self.eia930.df.abs() < 1.5] = 0

        # Emission factors for non-US bas
        self.default_factors = get_average_emission_factors(prefix)

        # Look up lists of BAs with specific requirements
        self.import_regions, self.generation_regions = self._get_special_regions()

        # Load generated pollutant mass emissions and generation for import-only BAs
        self.pollutant_mass_emissions, self.generation = (
            self._load_pollutant_mass_emissions_and_generation_for_import_only_regions()
        )

        # Identify shared BAs
        regions = set(self.eia930.regions)
        regions = regions.intersection(set(self.generation.columns))
        # Add back import-only regions
        regions = regions.union(set(self.import_regions))
        self.regions = list(regions)

        # Build result df
        self.results = self._build_results()

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
            if b in self.eia930.regions
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
            self.results[ba]["net_consumed_mwh"] = self.eia930.df[KEYS["E"]["D"] % ba][
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
                    gen_cols = [(src, KEYS["E"]["SRC_%s" % src] % ba) for src in SRC]
                    gen_cols = [
                        (src, col)
                        for src, col in gen_cols
                        if col in self.eia930.df.columns
                    ]
                    all_ba_mass_emissions.loc[:, ba] = self.eia930.df.apply(
                        lambda x: sum(
                            self.default_factors[pol][adj][src] * x[col]
                            for src, col in gen_cols
                        ),
                        axis=1,
                    )

                # Cut off emissions at 8 hours after UTC year
                all_ba_mass_emissions = all_ba_mass_emissions[
                    : f"{self.year + 1}-01-01 08:00:00+00:00"
                ]
                mass_emissions[(adj, pol)] = all_ba_mass_emissions

        # Make generation data frame
        generation = pd.DataFrame(data=ba_generation)
        generation = generation[: f"{self.year + 1}-01-01 08:00:00+00:00"]

        return mass_emissions, generation

    def build_matrices(
        self, pol: str, adj: str, date: pd.Timestamp
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build interchange matrix, emission vector, and generation vector for a given
        timestamp.

        Args:
            pol (str): pollutant name.
            adj (str): adjustment type.
            date (pd.Timestamp): timestamp for which to build matrices.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: interchange matrix, emission vector,
                and generation vector.
        """
        # Build transmission matrix from cleaned EIA-930 data
        ID = np.zeros((len(self.regions), len(self.regions)))
        for i, ri in enumerate(self.regions):
            for j, rj in enumerate(self.regions):
                if KEYS["E"]["ID"] % (ri, rj) in self.eia930.df.columns:
                    ID[i][j] = self.eia930.df.loc[date, KEYS["E"]["ID"] % (ri, rj)]

        # Build pollutant mass vector, using default per-fuel factors for import-only
        # regions.
        E = self.pollutant_mass_emissions[(adj, pol)].loc[date, self.regions].to_numpy()

        # Build generation vector, using 930 for import-only regions
        G = np.zeros(len(self.regions))
        for i, r in enumerate(self.regions):
            if r in self.import_regions:
                G[i] = self.eia930.df.loc[date, KEYS["E"]["NG"] % r]
            else:
                G[i] = self.generation.loc[date, r]

        E = np.nan_to_num(E)
        G = np.nan_to_num(G)
        ID = np.nan_to_num(ID)

        # In some cases, we have zero generation but non-zero transmission
        # usually due to imputed zeros during physics-based cleaning being set to 1.0
        # but sometimes due to ok values being set to 1.0
        to_fix = (ID.sum(axis=1) > 0) & (G == 0)
        ID[:, to_fix] = 0
        ID[to_fix, :] = 0

        return E, G, ID

    def run(self):
        """Orchestration method to calculate consumed emission rates"""
        for pol in POLLUTANTS:
            for adj in ADJUSTMENTS:
                total_failed = 0
                col = get_pollutant_rate_column_name(
                    pol, adjustment=adj, generated=False
                )
                logger.info(f"Solving consumed {pol} {adj} emissions...")
                # Calculate emissions
                for date in self.generation.index:
                    if self.small and (
                        date - self.generation.index[0] > pd.Timedelta(weeks=1)
                    ):
                        break  # only calculate one week if small run

                    E, G, ID = self.build_matrices(pol, adj, date)
                    # If we have NaNs, don't bother
                    if (
                        (np.isnan(G).sum() > 0)
                        or (np.isnan(E).sum() > 0)
                        or (np.isnan(ID).sum() > 0)
                    ):
                        consumed_emissions = np.full(len(self.regions), np.nan)
                    else:
                        # Run
                        try:
                            consumed_emissions, _ = consumption_emissions(E, G, ID)
                        except np.linalg.LinAlgError:
                            # These issues happen at boundary hours (beginning and end
                            # of year) where we don't have full data for all BAs
                            total_failed += 1
                            consumed_emissions = np.full(len(self.regions), np.nan)

                    # Export
                    for i, r in enumerate(self.regions):
                        self.results[r].loc[date, col] = consumed_emissions[i]
                if total_failed > 0:
                    logger.warning(
                        f"{total_failed} hours failed to solve for consumed "
                        f"{pol} {adj} emissions."
                    )
