import numpy as np
import pandas as pd
import os
import sys

from gridemissions.load import BaData
from gridemissions.eia_api import KEYS, SRC
from filepaths import outputs_folder, manual_folder, results_folder

from output_data import (
    GENERATED_EMISSION_RATE_COLS,
    CONSUMED_EMISSION_RATE_COLS,
    output_to_results,
    TIME_RESOLUTIONS,
)

""" For these BAs, there are significant and systematic differences 
between our net_generation_mwh and EIA-930 net generation and interchange, 
so we cannot combine our net generation and 930 interchange to get net_consumed. 
Instead, we use 930 demand as net_consumed. Note: there may be issues with the 930 
demand! But it is better than combining inconsistent generation and interchange, 
which results in unreasonable profiles with many negative hours.
"""
# original values ["SPA", "CPLW", "GCPD", "AZPS", "EEI"]
BA_930_INCONSISTENCY = {
    2019: ["CPLW", "EEI"],
    2020: ["CPLW", "EEI", "SEC"],
    2021: ["CPLW", "GCPD"],
}

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


def get_column(poll: str, adjustment: str, ba: str = ""):
    """
    Return output file column name for a poll and adjustment type
    Returns mass columns, not rate columns
    """
    assert poll in POLLUTANTS
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


def get_average_emission_factors(prefix: str, year: int):
    """
    Locate per-fuel, per-adjustment, per-pollutant emission factors.
    Used to fill in emissions from BAs outside of US, where we have generation by
    fuel (from gridemissions) but no open-grid-emissions data

    We use `gridemissions` assumptions for fuel mix for non-US BAs, which are simple and not time-varying

    Structure: EMISSIONS_FACTORS[poll][adjustment][fuel]
    """
    genavg = pd.read_csv(
        outputs_folder(f"{prefix}annual_generation_averages_by_fuel_{year}.csv"),
        index_col="fuel_category",
    )
    efs = {}
    for pol in POLLUTANTS:
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


class HourlyConsumed:
    """
        `HourlyConsumed`
    Class to load data, calculate consumed rates, and output per-BA
    """

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
        self.default_factors = get_average_emission_factors(prefix, year)

        # Look up lists of BAs with specific requirements
        self.import_regions, self.generation_regions = self._get_special_regions()

        # Load generated rates, save to self.generated
        self.rates, self.generation = self._load_rates()

        # Identify shared BAs
        regions = set(self.eia930.regions)
        regions = regions.intersection(set(self.generation.columns))
        # Add back import-only regions
        regions = regions.union(set(self.import_regions))
        self.regions = list(regions)

        # Build result df
        self.results = self._build_results()

    def _get_special_regions(self):
        """
        Get import only regions:
            Get regions that interchange with US regions but whose generation and emissions
            we need to fill from EIA-930
        And generation only regions:
            We won't export files for these
        """
        self.ba_ref = pd.read_csv(
            manual_folder("ba_reference.csv"), index_col="ba_code"
        )
        generation_only = list(
            self.ba_ref[self.ba_ref.ba_category == "generation_only"].index
        )

        # Get import-only regions
        import_only = [
            b
            for b in self.ba_ref[self.ba_ref["us_ba"] == "No"].index
            if b in self.eia930.regions
        ]
        return import_only, generation_only

    def _build_results(self):
        """Builds result dataframe per output file."""
        results = {}
        cols = []
        for pol in POLLUTANTS:
            for adj in ADJUSTMENTS:
                cols.append(get_rate_column(pol, adjustment=adj, generated=False))
                cols.append(get_column(pol, adjustment=adj))
        cols.append("net_consumed_mwh")
        for ba in self.regions:
            results[ba] = pd.DataFrame(
                index=self.generation.index, columns=cols, dtype=np.float64
            )
        return results

    def output_results(self):
        """
            HourlyConsumed.output_results
        After running HourlyConsumed.run(), results will be saved in a map of BA -> result df
        Only rate is calculated with matrix calc, other cols calced here:
            * Consumed elec is calculated from 930 total interchange + our gen estimate
            * Consumed carbon is calculated as consumed elec * consumed CI
        Here we output each df to a file in `carbon_accounting`

        Note that we are calculating consumed carbon and MWh so can aggregate correctly,
        but we are dropping from final outputs for simplicity.
        """
        for ba in self.regions:
            if (ba in self.import_regions) or (ba in self.generation_regions):
                continue
            if ba in BA_930_INCONSISTENCY[self.year]:
                print(f"Using D instead of (G-TI) for consumed calc in {ba}")
                self.results[ba]["net_consumed_mwh"] = self.eia930.df[
                    KEYS["E"]["D"] % ba
                ][self.generation.index]
            else:
                self.results[ba]["net_consumed_mwh"] = (
                    self.generation[ba] - self.eia930.df[KEYS["E"]["TI"] % ba]
                )[self.generation.index]
            for pol in POLLUTANTS:
                for adj in ADJUSTMENTS:
                    self.results[ba][get_column(pol, adjustment=adj)] = (
                        self.results[ba][
                            get_rate_column(pol, adjustment=adj, generated=False)
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
                    time_dat.datetime_local.dt.year == self.year
                ]  # keep year of local data

                if time_resolution == "hourly":
                    # No resampling needed; keep timestamp cols in output
                    time_cols = ["datetime_utc", "datetime_local"]
                    missing_hours = time_dat[time_dat.isna().any(axis=1)]
                    if len(missing_hours) > 0:
                        print(
                            f"WARNING: {len(missing_hours)} hours are missing in {ba} consumed data"
                        )
                elif time_resolution == "monthly":
                    time_dat["month"] = time_dat.datetime_local.dt.month
                    # Aggregate to appropriate resolution
                    time_dat = (
                        time_dat.groupby("month")
                        .sum()[EMISSION_COLS + ["net_consumed_mwh"]]
                        .reset_index()  # move "month" to column
                    )
                    time_cols = ["month"]
                elif time_resolution == "annual":
                    time_dat["year"] = time_dat.datetime_local.dt.year
                    # Aggregate to appropriate resolution
                    time_dat = (
                        time_dat.groupby("year")
                        .sum()[EMISSION_COLS + ["net_consumed_mwh"]]
                        .reset_index()  # move "year" to column
                    )
                    time_cols = ["year"]

                # Calculate rates from summed emissions, consumption
                for pol in POLLUTANTS:
                    for adj in ADJUSTMENTS:
                        rate_col = get_rate_column(pol, adj, generated=False)
                        emission_col = get_column(pol, adj)
                        time_dat[rate_col] = (
                            time_dat[emission_col] / time_dat["net_consumed_mwh"]
                        )

                # Output
                output_to_results(
                    time_dat[time_cols + CONSUMED_EMISSION_RATE_COLS],
                    ba,
                    f"/carbon_accounting/{time_resolution}/",
                    self.prefix,
                    skip_outputs=self.skip_outputs,
                )
        return

    def _impute_border_hours(self, temp):
        """
        Add three hours to beginning and end of series.
        Impute hours by taking same hour from previous (or next) day of series
        This matches EIA-930's imputation approach
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

    def _load_rates(self):
        # Load all rates
        rates = {}  # (adj, pol) -> {(BA, rate series)}
        gens = {}
        for f in os.listdir(
            results_folder(f"{self.prefix}/power_sector_data/hourly/us_units/")
        ):
            if ".DS_Store" in f:
                continue
            this_ba = pd.read_csv(
                results_folder(f"{self.prefix}/power_sector_data/hourly/us_units/") + f,
                index_col="datetime_utc",
                parse_dates=True,
            )
            this_ba = this_ba[this_ba.fuel_category == "total"]
            ba_name = f.replace(".csv", "")
            for adj in ADJUSTMENTS:
                for pol in POLLUTANTS:
                    this_rate = rates.get((adj, pol), {})
                    this_rate[ba_name] = self._impute_border_hours(
                        this_ba[get_column(pol, adjustment=adj)]
                    )
                    rates[(adj, pol)] = this_rate
            gens[ba_name] = self._impute_border_hours(this_ba["net_generation_mwh"])

        # Make each rate into a DF and add emissions for import-only regions
        for pol in POLLUTANTS:
            for adj in ADJUSTMENTS:
                # Most rates we already loaded above
                emissions = pd.DataFrame(rates[(adj, pol)])

                # Add import regions to emissions DF
                for ba in self.import_regions:
                    gen_cols = [(src, KEYS["E"]["SRC_%s" % src] % ba) for src in SRC]
                    gen_cols = [
                        (src, col)
                        for src, col in gen_cols
                        if col in self.eia930.df.columns
                    ]
                    emissions.loc[:, ba] = self.eia930.df.apply(
                        lambda x: sum(
                            self.default_factors[pol][adj][src] * x[col]
                            for src, col in gen_cols
                        ),
                        axis=1,
                    )

                # Cut off emissions at 9 hours after UTC year
                emissions = emissions[:f"{self.year+1}-01-01 09:00:00+00:00"]
                rates[((adj, pol))] = emissions

        # Make generation data frame
        generation = pd.DataFrame(data=gens)
        generation = generation[:f"{self.year+1}-01-01 09:00:00+00:00"]

        return rates, generation

    def build_matrices(self, pol: str, adj: str, date):
        # return emissions, interchange, and generation
        # Build transmission matrix from cleaned
        ID = np.zeros((len(self.regions), len(self.regions)))
        for i, ri in enumerate(self.regions):
            for j, rj in enumerate(self.regions):
                if KEYS["E"]["ID"] % (ri, rj) in self.eia930.df.columns:
                    ID[i][j] = self.eia930.df.loc[date, KEYS["E"]["ID"] % (ri, rj)]

        # Build emission array, using default per-fuel factors for import-only regions (above)
        E = self.rates[(adj, pol)].loc[date, self.regions].to_numpy()

        # Build generation array, using 930 for import-only regions
        G = np.zeros(len(self.regions))
        for (i, r) in enumerate(self.regions):
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
        for pol in POLLUTANTS:
            for adj in ADJUSTMENTS:
                total_failed = 0
                col = get_rate_column(pol, adjustment=adj, generated=False)
                print(f"{pol}, {adj}", end="...")
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
                            # These issues happen at boundary hours (beginning and end of year)
                            # where we don't have full data for all BAs
                            total_failed += 1
                            consumed_emissions = np.full(len(self.regions), np.nan)

                    # Export
                    for (i, r) in enumerate(self.regions):
                        self.results[r].loc[date, col] = consumed_emissions[i]
                if total_failed > 0:
                    print(
                        f"Warning: {total_failed} hours failed to solve for consumed {pol} {adj} emissions."
                    )
