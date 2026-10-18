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
    emissions: np.ndarray, generation: np.ndarray, interchange: np.ndarray
) -> tuple[np.ndarray, int]:
    """Solve for the consumed emission rate of every region in one hour.

    Adapted from gridemissions (https://github.com/jdechalendar/gridemissions), which
    implements the multi-regional input-output method of de Chalendar et al. (2019),
    https://www.pnas.org/doi/full/10.1073/pnas.1912950116 (Equations 1 to 4).

    Electricity available in region i comes from its own generation and from its
    imports, and the emissions embodied in it come from region i's own emissions and
    the emissions embodied in those imports. The consumed emission rate x_i of each
    region therefore solves:

        x_i * (p_i + U_i) - sum_j u_ij * x_j = f_i

    where (using the notation of the consumption-based emissions methodology docs):
        x_i: consumed emission rate of region i (mass per MWh)
        p_i: net generation in region i (MWh)
        u_ij: electricity imported into region i from region j (MWh)
        U_i: total imports into region i, sum_j u_ij (MWh)
        f_i: emissions from generation in region i (mass)

    In matrix form this is `system_matrix @ x = emissions`, with
    `system_matrix = diag(generation + total_imports) - imports`.

    The system matrix only depends on generation and interchange, so several emission
    series (e.g. different pollutants) for the same hour can be solved at once by
    passing them as the columns of `emissions`.

    Args:
        emissions (np.ndarray): emissions from generation in each region, as a vector
            (n_regions,), or as a matrix (n_regions, n_series) with one column per
            emission series.
        generation (np.ndarray): net generation in each region (n_regions,).
        interchange (np.ndarray): interchange between regions (n_regions, n_regions),
            where `interchange[i, j]` is the net flow from region i to region j as
            reported by region i: positive for exports, negative for imports.

    Returns:
        tuple[np.ndarray, int]: consumed emission rates, with the same shape as
            `emissions`, and the number of isolated regions whose rate was set to zero.

    Raises:
        np.linalg.LinAlgError: if the system matrix is singular, e.g. a region with no
            generation and no imports that other regions report importing from.
        ValueError: if an isolated region's solved rate isn't zero.
    """

    # numpy version must be high enough, otherwise np.linalg.cond fails on a matrix with
    # only zeros.
    assert Version(np.__version__) >= Version("1.15.1")

    # interchange reports flows from each region's point of view; keep only the
    # imports, as non-negative flows into region i from region j
    imports = (-interchange).clip(min=0)
    total_imports = imports.sum(axis=1)

    # electricity available in each region is its generation plus its imports, and
    # each import brings in the consumed emission rate of the region it comes from
    system_matrix = np.diag(generation + total_imports) - imports
    # copy so that zeroing isolated regions below doesn't modify the caller's array
    rhs = np.array(emissions, dtype=np.float64)

    # A region with no generation, no imports, and no other region importing from it
    # has an all-zero row and column, which makes the system singular. When the
    # system is ill-conditioned, give each such isolated region the equation
    # x_i = 0 so that the rest of the system can be solved.
    isolated_regions = []
    if np.linalg.cond(system_matrix) > (1.0 / sys.float_info.epsilon):
        for i in range(len(system_matrix)):
            if (np.abs(system_matrix[:, i]).sum() == 0.0) & (
                np.abs(system_matrix[i, :]).sum() == 0.0
            ):
                system_matrix[i, i] = 1.0
                isolated_regions += [i]
                # zero the region's emissions in every series to match x_i = 0
                rhs[i] = 0.0

    consumed_rates = np.linalg.solve(system_matrix, rhs)

    # isolated regions don't interact with the rest of the system, so their solved
    # rate must be exactly the zero they were given
    for i in isolated_regions:
        if np.any(consumed_rates[i] != 0.0):
            logger.warning(f"emissions of region {i}: {rhs[i]}")
            logger.warning(
                f"sum of abs(system_matrix[{i}, :]): "
                f"{np.abs(system_matrix[i, :]).sum()}"
            )
            logger.warning(
                f"sum of abs(system_matrix[:, {i}]): "
                f"{np.abs(system_matrix[:, i]).sum()}"
            )
            raise ValueError(
                f"consumed rate of isolated region {i} is {consumed_rates[i]} "
                "instead of 0"
            )

    return consumed_rates, len(isolated_regions)


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
        """Build the inputs of `consumption_emissions` for every hour at once.

        Each array covers every hour of `self.generation.index` (first axis) and every
        region of `self.regions` (in that order), so the inputs for hour t are
        `emissions_by_hour[t]`, `generation_by_hour[t]` and `interchange_by_hour[t]`.
        Missing values are treated as zero.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]:
                * emissions (hours, regions, series): emissions from generation in each
                  region, with one series per entry of `POLLUTANT_ADJUSTMENTS`. Uses
                  default per-fuel emission factors for import-only regions.
                * generation (hours, regions): net generation in each region. Uses
                  EIA-930 net generation for import-only regions.
                * interchange (hours, regions, regions): `[t, i, j]` is the net flow
                  from region i to region j reported by region i (positive for exports)
                  in EIA-930. Zero for pairs of regions that don't trade.
        """
        hours = self.generation.index
        n_regions = len(self.regions)
        region_index = {region: i for i, region in enumerate(self.regions)}

        # Build interchange arrays from cleaned EIA-930 data. Only pairs of regions
        # that both take part in the linear system and have interchange data are
        # filled; every other pair stays zero.
        region_pairs = [
            (from_region, to_region)
            for (from_region, to_region) in self.interchange.columns
            if (from_region in region_index) and (to_region in region_index)
        ]
        interchange_by_hour = np.zeros((len(hours), n_regions, n_regions))
        if len(region_pairs) > 0:
            from_index = [region_index[from_region] for from_region, _ in region_pairs]
            to_index = [region_index[to_region] for _, to_region in region_pairs]
            interchange_by_hour[:, from_index, to_index] = self.interchange.loc[
                hours, region_pairs
            ].to_numpy()

        # Build emission arrays, using default per-fuel factors for import-only
        # regions, with one series per (pollutant, adjustment) along the last axis.
        emissions_by_hour = np.stack(
            [
                self.pollutant_mass_emissions[(adj, pol)]
                .loc[hours, self.regions]
                .to_numpy()
                for pol, adj in POLLUTANT_ADJUSTMENTS
            ],
            axis=-1,
        )

        # Build generation arrays, using 930 for import-only regions
        generation_by_hour = np.zeros((len(hours), n_regions))
        for i, region in enumerate(self.regions):
            if region in self.import_regions:
                generation_by_hour[:, i] = self.net_generation.loc[hours, region]
            else:
                generation_by_hour[:, i] = self.generation.loc[hours, region]

        # treat missing values as zero so that every hour has a complete system
        emissions_by_hour = np.nan_to_num(emissions_by_hour)
        generation_by_hour = np.nan_to_num(generation_by_hour)
        interchange_by_hour = np.nan_to_num(interchange_by_hour)

        # In some cases, a region has zero generation but reports net exports, which
        # can't be traced to any emissions. For those region-hours, drop the region's
        # row and column of the interchange matrix so that it doesn't trade.
        zero_generation_exporters = (interchange_by_hour.sum(axis=2) > 0) & (
            generation_by_hour == 0
        )
        interchange_by_hour[
            zero_generation_exporters[:, :, np.newaxis]
            | zero_generation_exporters[:, np.newaxis, :]
        ] = 0.0

        return emissions_by_hour, generation_by_hour, interchange_by_hour

    def run(self):
        """Calculate hourly consumed emission rates for every region and series.

        Solves one linear system per hour with `consumption_emissions`, for all
        (pollutant, adjustment) series at once, and writes the rates to the
        consumed rate columns of `self.results`. Hours whose system can't be solved
        are left as NaN.
        """
        logger.info("Solving consumed emissions for all pollutants and adjustments...")
        emissions_by_hour, generation_by_hour, interchange_by_hour = (
            self.build_matrices()
        )

        # Solve each hour's linear system once for all emission series, since the
        # system matrix only depends on generation and interchange
        consumed_rates = np.full(emissions_by_hour.shape, np.nan)
        total_failed = 0
        for t in range(len(generation_by_hour)):
            # If we have NaNs, don't bother
            if (
                np.isnan(generation_by_hour[t]).any()
                or np.isnan(emissions_by_hour[t]).any()
                or np.isnan(interchange_by_hour[t]).any()
            ):
                continue
            try:
                consumed_rates[t], _ = consumption_emissions(
                    emissions_by_hour[t], generation_by_hour[t], interchange_by_hour[t]
                )
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
            for i, region in enumerate(self.regions):
                self.results[region][col] = consumed_rates[:, i, s]
