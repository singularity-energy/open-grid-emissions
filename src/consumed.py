import numpy as np
import pandas as pd
import os

from gridemissions.emissions import BaDataEmissionsCalc, EMISSIONS_FACTORS
from gridemissions.load import BaData
from gridemissions.eia_api import KEYS, SRC

from src.output_data import (
    GENERATED_EMISSION_RATE_COLS,
    output_to_results,
    TIME_RESOLUTIONS,
)

# # From workers/consumed_carbon.py: factors for inport-only regions
# DEFAULT_EMISSION_FACTORS = {
#     "CO2": {
#         "for_electricity": {
#             "IESO": 79.37,  # from a report published by ISONE
#             "HQT": 2.2,
#             "AESO": 1350.5,  # default_mix + emissions API for EGRID adjusted 2018
#             "NBSO": 546.17,  # default_mix + emissions API for EGRID adjusted 2018
#             "BCHA": 116.66,  # default_mix + emissions API for EGRID adjusted 2018
#             "MHEB": 16.7,  # default_mix + emissions API for EGRID adjusted 2018
#             "CFE": 949.37,  # default_mix + emissions API for EGRID adjusted 2018
#             "CEN": 949.37,
#         },  # default_mix + emissions API for EGRID adjusted 2018
#         "adjusted": {
#             "IESO": 79.37,  # from a report published by ISONE
#             "HQT": 2.2,
#             "AESO": 1350.5,  # default_mix + emissions API for EGRID adjusted 2018
#             "NBSO": 546.17,  # default_mix + emissions API for EGRID adjusted 2018
#             "BCHA": 116.66,  # default_mix + emissions API for EGRID adjusted 2018
#             "MHEB": 16.7,  # default_mix + emissions API for EGRID adjusted 2018
#             "CFE": 949.37,  # default_mix + emissions API for EGRID adjusted 2018
#             "CEN": 949.37,
#         },  # default_mix + emissions API for EGRID adjusted 2018
#     }
# }

# Regions outside the US that are import-only.
IMPORT_REGIONS = [
    "IESO",
    "HQT",
    "AESO",
    "NBSO",
    "BCHA",
    "MHEB",
    "CFE",
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

POLLUTANTS = ["CO2", "CH4", "N2O", "CO2E", "NOX", "SO2"]

ADJUSTMENTS = ["for_electricity", "for_electricity_adjusted"]

# Unused by us, but parent class wants to know
for pol in POLLUTANTS:
    EMISSIONS_FACTORS[pol] = np.nan


def get_column(pollutant: str, adjustment: str, ba: str = ""):
    """
    Return output file column name for a pollutant and adjustment type
    Returns mass columns, not rate columns
    """
    assert pollutant in POLLUTANTS
    assert adjustment in ADJUSTMENTS
    if ba == "":  # no BA, looking for output file column
        column = pollutant.lower() + "_mass_lb_" + adjustment
        assert column in EMISSION_COLS
    else:
        column = ba + "_" + pollutant.lower() + "_mass_lb_" + adjustment
    return column


def get_rate_column(pollutant: str, adjustment: str, generated: bool = True, ba: str = ""):
    """
    Return either generated or consumed output file rate column
    for pollutant `pollutant` and adjustment `adjustment`
    """
    assert pollutant in POLLUTANTS
    assert adjustment in ADJUSTMENTS
    if generated:
        column = "generated_" + pollutant.lower() + "_rate_lb_per_mwh_" + adjustment
        assert column in GENERATED_EMISSION_RATE_COLS
    else:
        column = "consumed_" + pollutant.lower() + "_rate_lb_per_mwh_" + adjustment
        assert column in CONSUMED_EMISSION_RATE_COLS
    if ba != "":  # For internal column use, add ba
        column = ba + "_" + column
    return column


# Lookups used by BaDataEmissionsCalc to build matrix calculation
# Add some keys for our unique pollutants
KEYS["CH4"] = {
    "D": "CH4_%s_D",
    "NG": "CH4_%s_NG",
    "TI": "CH4_%s_TI",
    "ID": "CH4_%s-%s_ID",
}
KEYS["N2O"] = {
    "D": "N2O_%s_D",
    "NG": "N2O_%s_NG",
    "TI": "N2O_%s_TI",
    "ID": "N2O_%s-%s_ID",
}
KEYS["CO2E"] = {
    "D": "CO2E_%s_D",
    "NG": "CO2E_%s_NG",
    "TI": "CO2E_%s_TI",
    "ID": "CO2E_%s-%s_ID",
}


def get_average_emission_factors(prefix: str = "", year: int = 2020):
    """
    Edit EMISSIONS dict with per-fuel, per-adjustment, per-pollutant emission factors.
    Used to fill in emissions from BAs outside of US, where we have generation by
    fuel (from gridemissions) but no hourly-egrid data

    Structure: EMISSIONS_FACTORS[pollutant][adjustment][fuel]
    """
    genavg = pd.read_csv(
        f"../data/outputs/{prefix}annual_generation_averages_by_fuel_{year}.csv",
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


class HourlyBaDataEmissionsCalc(BaDataEmissionsCalc):
    def __init__(
        self,
        path_930,
        pollutant="CO2",
        adjustment="for_electricity",
        year: int = 2020,
        small: bool = False,
        path_prefix: str = "",
    ):
        """
        TODO: take in out types, convert to gridemissions types
        """
        self.original = BaData(
            path_930
        )  # todo pass a dataframe with our generation columns, per-fuel generation dropped
        self.year = year
        self.small = small
        self.prefix = path_prefix
        self.pollutant = pollutant
        self.adjustment = adjustment  # "for_electricity" or "adjusted"

        data, rates, regions = self._replace_generation()
        self.rates = rates
        fixed = BaData(df=data)
        fixed.regions = regions  # make sure we only use our data regions. TODO: not needed if clean up cols in _replace_generation

        # Overwrite emission factors object
        self.emissions_factors = get_average_emission_factors(self.prefix, self.year)

        super().__init__(fixed, pollutant)

    def _drop_pol_cols(self, pollutant):
        """
        For repeated processing with different pols, need to drop pollutant columns

        """
        pol_cols = [c for c in self.df.columns if c[0 : len(pollutant)] == pollutant]
        self.df.drop(columns=pol_cols, inplace=True)

    def process(self):
        """
        Run process for all adjustments and pols
        """
        # Set up outputs: map from BA to timeseries of consumed emission rates
        outputs = {}
        for ba in self.output_regions:
            outputs[ba] = pd.DataFrame(
                index=self.df.index,
                columns=EMISSION_COLS + CONSUMED_EMISSION_RATE_COLS,
            )
        # Run process for all adjustments and pollutants

        for pol in POLLUTANTS:
            for adj in ADJUSTMENTS:
                self.adjustment = adj
                self.pollutant = pol
                self.KEY_pollutant = KEYS[pol]
                print(f"Running {adj}, {pol}")
                super().process()
                # After processing, add columns to output
                for ba in self.output_regions:
                    rate_col = get_rate_column(pol, adj, generated=False)
                    emission_col = get_column(pol, adj)
                    outputs[ba].loc[:, rate_col] = self.df["%si_%s_D" % (pol, ba)]
                    outputs[ba].loc[:, emission_col] = self.df["%s_%s_D" % (pol, ba)]
                # Drop consumed rate columns (will overwrite with same pol, next adjustment)
                print("dropping cols")
                self._drop_pol_cols(pol)
        # Add demand column to output
        for ba in self.output_regions:
            outputs[ba].loc[:, "net_consumed_mwh"] = self.df[KEYS["E"]["D"] % ba]
        self.output_dat = outputs
        return outputs

    def output_data(self, path_prefix: str):
        """
        Run after process.
        Follows output_data format for results files.
        Outputs per-ba consumed emissions and rate data
        """
        for ba in self.output_regions:
            dat = self.output_dat[ba]
            dat = dat.reset_index()
            dat = dat.rename(columns={"index": "datetime_utc"})
            for time_resolution in TIME_RESOLUTIONS:
                time_dat = dat.copy(deep=True)

                # Aggregate to appropriate resolution
                time_dat = time_dat.resample(
                    TIME_RESOLUTIONS[time_resolution],
                    closed="left",
                    label="left",
                    on="datetime_utc",
                ).sum()[EMISSION_COLS + ["net_consumed_mwh"]]

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
                    dat, ba, f"carbon_accounting/{time_resolution}/", path_prefix
                )
        return

    def _replace_generation(self):
        """
        Helper function to set up generation and rate data.
        1) Find list of regions to use
        2) Drop columns: all columns from regions we are not using, all generation columns
        3) Load default values from eGRID for zero rates (TODO: Fix once we know that zero rates are real)
        4/5) Load hourly generation data into BaData structure using KEYS as columns
             Load hourly rate data into rates data structure using `<ba>_<GENERATED_EMISSION_RATE_COLS>` as columns

        TODO: make demand sum of generation and net interchange
        """

        # 1: Find region list
        # We have some regions not reported in 930 data ('AMPL', 'CEA', 'CSTO', 'GRIS', 'HECO')
        # Drop all these, use intersection. TODO: Ensure `gridemissions` BA list is up to date
        regions930 = set(self.original.regions)
        our_regions = {
            f.replace(".csv", "")
            for f in os.listdir(
                f"../data/results/{self.prefix}power_sector_data/hourly/us_units/"
            )
        }
        self.output_regions = regions930.intersection(
            our_regions
        )  # will want to output just regions with carbon_accounting csvs
        our_regions.update(
            set(IMPORT_REGIONS)
        )  # Add input-only regions (outside of US)
        regions_to_use = regions930.intersection(our_regions)
        regions_to_drop = (
            regions930 - our_regions
        )  # Regions included in 930 output but not in our regions

        # 2: Drop columns
        data = self.original.df.copy(deep=True)  # start with original data
        cols_to_drop = []
        for col in data.columns:
            ba1 = col.split(".")[1].split("-")[0]  # primary column BA
            ba2 = col.split(".")[1].split("-")[1]  # "ALL" unless this is an ID col
            if ("NG" in col) and (
                ba1 not in IMPORT_REGIONS
            ):  # Drop net generation columns, will use our own data
                cols_to_drop.append(col)
                continue
            if (ba1 in regions_to_drop) or (
                ba2 in regions_to_drop
            ):  # Drop BAs we don't use
                cols_to_drop.append(col)
        data = data.drop(columns=cols_to_drop)
        data = data[data.index.year == self.year]
        if self.small:
            data = data[
                (data.index.month == 2) & (data.index.day == 3)
            ]  # just one day (24 hours)

        # 3: Load default values
        # For now, fill zero rates with eGRID annual average
        # TODO not a longterm solution!
        egrid = pd.read_excel(
            f"../data/downloads/egrid/egrid{2020}_data.xlsx",
            sheet_name=f"BA{str(2020)[-2:]}",
            header=1,
            index_col="BACODE",
            usecols=["BACO2RTA", "BACODE"],
        )
        egrid = egrid.replace(0, np.nan)
        avg = egrid[
            "BACO2RTA"
        ].mean()  # some eGRID rates are zero -- make these eGRID overall avg
        egrid = egrid.replace(np.nan, avg)

        # 4/5: For each BA, load hourly rate and generation from our data into BaDataEmissionsCalc data structure
        rates_cols = {}
        for ba in regions_to_use:
            # Import-only regions. Use 930 data and constant emission factors
            if (
                ba in IMPORT_REGIONS
            ):  # Skip. Did not drop these generation cols above, no rate cols
                continue

            # Normal case: load hourly generation and emission factor data
            else:
                new = pd.read_csv(
                    f"../data/results/{self.prefix}power_sector_data/hourly/us_units/{ba}.csv",
                    index_col="datetime_utc",
                    parse_dates=True,
                )
                new = new[
                    new["fuel_category"] == "total"
                ]  # Don't want per-fuel gen and rate
                data.loc[:, f"EBA.{ba}-ALL.NG.H"] = new["net_generation_mwh"]
                # Some BAs have large regions of zeros in our data that are not zero in 930 data. Fill. TODO: fix our data
                zeros = data.index[data[f"EBA.{ba}-ALL.NG.H"] == 0]
                data.loc[zeros, f"EBA.{ba}-ALL.NG.H"] = self.original.df.loc[
                    zeros, f"EBA.{ba}-ALL.NG.H"
                ]

                for pollutant in POLLUTANTS:
                    for adjustment in ADJUSTMENTS:
                        col_name = get_rate_column(pollutant, adjustment, ba=ba)
                        rates_cols[col_name] = new[
                            get_rate_column(self.pollutant, self.adjustment, generated=True)
                        ]
                        # TODO Some of our rates are eroneously zero. Replace with eGRID rates for now
                        rates_cols[col_name].replace(
                            0, egrid.get(ba, avg), inplace=True
                        )

        return data, pd.concat(rates_cols, axis="columns"), regions_to_use

    # TODO allow for use of parent data when missing our data
    def _add_production_emissions(self):
        """

        Overwrite production emission calculation
        to use hourly emission rates

        Creates columns with name `self.KEY_pollutant["NG"] % ba` with hourly rates
        """
        for ba in self.regions:
            # For import regions, we don't use our own data, so use same EF logic as parent class.
            if ba in IMPORT_REGIONS:
                gen_cols = [(src, self.KEY_E["SRC_%s" % src] % ba) for src in SRC]
                gen_cols = [
                    (src, col) for src, col in gen_cols if col in self.df.columns
                ]
                self.df.loc[:, self.KEY_pollutant["NG"] % ba] = self.df.apply(
                    lambda x: sum(
                        self.emissions_factors[self.pollutant][self.adjustment][src] * x[col]
                        for src, col in gen_cols
                    ),
                    axis=1,
                )
            # Otherwise, use hourly rates and net generation numbers
            else:
                rate_col = get_rate_column(self.pollutant, self.adjustment, ba=ba)
                self.df.loc[:, self.KEY_pollutant["NG"] % ba] = (
                    self.rates[rate_col] * self.df[f"EBA.{ba}-ALL.NG.H"]
                )
