import numpy as np
import pandas as pd
import os

from gridemissions.emissions import BaDataEmissionsCalc
from gridemissions.load import BaData

# From workers/consumed_carbon.py: factors for inport-only regions
DEFAULT_EMISSION_FACTORS = {
    # from a report published by ISONE
    "IESO": 79.37,
    "HQT": 2.2,
    "AESO": 1350.5,  # default_mix + emissions API for EGRID adjusted 2018
    "NBSO": 546.17,  # default_mix + emissions API for EGRID adjusted 2018
    "BCHA": 116.66,  # default_mix + emissions API for EGRID adjusted 2018
    "MHEB": 16.7,  # default_mix + emissions API for EGRID adjusted 2018
    "CFE": 949.37,  # default_mix + emissions API for EGRID adjusted 2018
    "CEN": 949.37,  # default_mix + emissions API for EGRID adjusted 2018
}


class HourlyBaDataEmissionsCalc(BaDataEmissionsCalc):
    def __init__(self, path_930, poll="CO2", year: int = 2020):
        """
        TODO: take in out types, convert to gridemissions types
        """
        self.original = BaData(
            path_930
        )  # todo pass a dataframe with our generation columns, per-fuel generation dropped
        self.year = year
        self.poll = poll

        data, rates, regions = self._replace_generation()
        self.rates = rates
        fixed = BaData(df=data)
        fixed.regions = regions  # make sure we only use our data regions. TODO: not needed if clean up cols in _replace_generation

        super().__init__(fixed, poll)

    def _replace_generation(self):
        # We have some regions not reported in 930 data ('AMPL', 'CEA', 'CSTO', 'GRIS', 'HECO')
        # gridemissions reports 930 data with some old regions ('AESO', 'BCHA', 'CEN', 'CFE', 'HQT', 'IESO', 'MHEB', 'SPC'),
        # since it also works for older years that have those regions
        # Drop all these, use intersection. TODO: Ensure `gridemissions` is up to date
        regions930 = set(self.original.regions)
        our_regions = {
            f.replace(".csv", "")
            for f in os.listdir("../data/results/carbon_accounting/")
        }
        self.output_regions = (
            our_regions  # will want to output just regions with carbon_accounting csvs
        )
        our_regions.update(
            set(DEFAULT_EMISSION_FACTORS.keys())
        )  # Add input-only regions
        regions_to_use = regions930.intersection(our_regions)
        regions_to_drop = regions930 - our_regions

        data = self.original.df.copy(deep=True)  # data that we'll edit
        # data = data.drop(columns=[c for c in data.columns if ("NG" in c)]) # drop all generation columns (per-fuel and total)
        cols_to_drop = []
        for col in data.columns:
            ba1 = col.split(".")[1].split("-")[0]  # primary column BA
            ba2 = col.split(".")[1].split("-")[1]  # "ALL" unless this is an ID col
            # if (ba1 in DEFAULT_EMISSION_FACTORS.keys()) or (
            #     ba2 in DEFAULT_EMISSION_FACTORS.keys()
            # ):
            #     continue  # Don't drop import-only columns, we will use default rates here
            if "NG" in col:
                cols_to_drop.append(col)
                continue
            if (ba1 in regions_to_drop) or (ba2 in regions_to_drop):
                cols_to_drop.append(col)
        data = data.drop(columns=cols_to_drop)
        data = data[data.index.year == self.year]

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

        rates = pd.DataFrame()
        for ba in regions_to_use:
            if ba in {
                "GRID",
                "GRMA",
                "HGMA",
                "NSB",
            }:  # BAs where our generation data is zero. TODO fix
                data[f"EBA.{ba}-ALL.NG.H"] = self.original.df[f"EBA.{ba}-ALL.NG.H"]
            elif ba in DEFAULT_EMISSION_FACTORS.keys():
                data[f"EBA.{ba}-ALL.NG.H"] = DEFAULT_EMISSION_FACTORS[ba]
            else:  # Load hourly emission factor data
                new = pd.read_csv(
                    f"../data/results/carbon_accounting/{ba}.csv",
                    index_col=0,
                    parse_dates=True,
                )
                data[f"EBA.{ba}-ALL.NG.H"] = new["net_generation_mwh_total"]
                zeros = data.index[data[f"EBA.{ba}-ALL.NG.H"] == 0]
                data.loc[zeros, f"EBA.{ba}-ALL.NG.H"] = self.original.df.loc[
                    zeros, f"EBA.{ba}-ALL.NG.H"
                ]
            # rates[f"{ba}_{self.poll}_RATE"] = new["co2_mass_lb_total"] / new["net_generation_mwh_total"]
            rates[f"{ba}_{self.poll}_RATE"] = (
                new["generated_co2_rate_lb_per_mwh_total"] / 2000
            )  # TODO also allow adjusted column 'adjusted_generated_co2_rate_lb_per_mwh_total'
            rates[f"{ba}_{self.poll}_RATE"].replace(0, egrid.get(ba, avg), inplace=True)

        return data, rates, regions_to_use

    # TODO allow for use of parent data when missing our data
    def _add_production_emissions(self):
        """

        Overwrite production emission calculation
        to use hourly emission rates

        Creates columns with name `self.KEY_poll["NG"] % ba` with hourly rates
        """
        for ba in self.regions:
            self.df.loc[:, self.KEY_poll["NG"] % ba] = (
                self.rates[f"{ba}_{self.poll}_RATE"] * self.df[f"EBA.{ba}-ALL.NG.H"]
            )
        # # Constant rates for import-only regions
        # for ba in DEFAULT_EMISSION_FACTORS.keys():
        #     self.df.loc[:, self.KEY_poll["NG"] % ba] = (
        #         DEFAULT_EMISSION_FACTORS[ba] * self.df[f"EBA.{ba}-ALL.NG.H"]
        #     )
