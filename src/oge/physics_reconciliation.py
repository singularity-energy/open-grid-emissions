"""Physics-based cleaning of EIA-930 data, ported from the `gridemissions` package.

This module is a copy of the parts of `gridemissions` that OGE uses to clean EIA-930
data: the basic cleaner, the rolling window filter, the optimization-based (cvxpy)
reconciliation, and the `BaData` container they operate on. It was copied from the
Singularity fork of `gridemissions` at commit
0cbe7620620f02ab91227940ef107b53bd1f9bb4
(https://github.com/singularity-energy/gridemissions), which is itself a fork of
https://github.com/jdechalendar/gridemissions.

The code is intentionally kept as close as possible to the original so that its
outputs can be compared byte for byte with the `gridemissions` package. Parts of the
package that OGE does not use (the EIA API scraper, the pyomo cleaner, the emissions
calculation, the AMPD and visualization tools, and the config file handling) were not
copied.

The original code is distributed under the following license:

    MIT License

    Copyright (c) 2020-2021 Jacques de Chalendar

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
"""

import os
from os.path import join
import time
import logging
import re
import pandas as pd
import numpy as np
from collections import defaultdict
import cvxpy as cp
import dask


################################################################################
# From gridemissions/eia_api.py
################################################################################

KEYS = {
    "E": {
        "D": "EBA.%s-ALL.D.H",
        "NG": "EBA.%s-ALL.NG.H",
        "TI": "EBA.%s-ALL.TI.H",
        "ID": "EBA.%s-%s.ID.H",
    },
}

SRC = ["COL", "NG", "NUC", "OIL", "OTH", "SUN", "UNK", "WAT", "WND", "GEO", "BIO"]

for src in SRC:
    KEYS["E"]["SRC_%s" % src] = "EBA.%s-ALL.NG." + src + ".H"


################################################################################
# From gridemissions/load.py
################################################################################


class BaData(object):
    """Class to handle BA-level data. The EBA class provides generation,
    consumption, the trade matrix and total interchange either at the BA or at
    the regional (IEA-defined) level. User guide:
        https://www.eia.gov/realtime_grid/docs/userguide-knownissues.pdf

    The class is a light wrapped around the pd.DataFrame object which is the
    format in which the underlying data are stored.

    The main purpose of the class is to provide convenience functions that make
    handling the data easier.

    Timestamps are in UTC.

    EBA data columns
    ----------------
    D: Demand
    NG: Net Generation
    TI: Total Interchange - (positive if exports)
    ID: Interchange with directly connected balancing authorities - (positive
        if exports)

    Consistency requirements
    ------------------------
    - Interchange data is antisymmetric: ID[i,j] == -ID[j,i]
    - Total trade with interchange data: TI == sum(ID[i,:])
    - Balance equation for total trade, demand, generation: TI + D == NG

    Methods
    -------
    get_cols(self, r) : generate column names for regions r for a given field.

    Attributes
    ----------
    regions : are in alphabetical order
    df : raw dataframe
    """

    def __init__(self, fileNm=None, df=None, variable="E"):
        """
        Initialize the BaData object

        There are two ways to initialize this object: by passing a `pd.DataFrame`,
        or by passing a file name from which to read the data using `pd.read_csv`.


        Parameters
        ----------
        fileNm: str, default None
            fileNm from which to read the data
        df: pd.DataFrame, default None
            data

        """
        self.logger = logging.getLogger("load")

        if df is not None:
            self.df = df
        else:
            self.df = pd.read_csv(fileNm, index_col=0, parse_dates=True)

        self.variable = variable
        self.regions = self._parse_data_cols()
        self.fileNm = fileNm
        self.KEY = KEYS[variable]

    def get_cols(self, r=None, field="D"):
        """
        Retrieve column name(s) corresponding to region(s) and a field

        Parameters
        ----------
        r: str or list of str, default None
            regions. If None, data is returned for all regions
        field: str
            field for which to load columns. Used to index in self.KEY

        Returns
        -------
        cols: list of str
        """
        if field not in self.KEY:
            raise ValueError(f"{field} not in str(list(self.KEY.keys()))")
        if field != "ID":
            if r is None:
                r = self.regions
            if isinstance(r, str):
                r = [r]
            return [self.KEY[field] % ir for ir in r]
        else:
            if r is None:
                r = self.regions
            if isinstance(r, str):
                r = [r]
            return [
                self.KEY[field] % (ir, ir2)
                for ir in r
                for ir2 in self.regions
                if self.KEY[field] % (ir, ir2) in self.df.columns
            ]

    def get_trade_partners(self, r):
        """
        Return list of regions that trade with a given region

        Parameter
        ---------
        r: str
            region for which to search
        """
        partners = []
        for r2 in self.regions:
            if (self.KEY["ID"] % (r, r2) in self.df.columns) and (
                self.KEY["ID"] % (r2, r) in self.df.columns
            ):
                partners += [r2]
        return partners

    def _parse_data_cols(self):
        """
        Checks:
        - Consistent number of regions for demand / generation / total
            interchange / trade matrix
        Returns the list of regions
        """
        regions = set([re.split(r"\.|-|_", el)[1] for el in self.df.columns])
        D_cols = [
            re.split(r"\.|-|_", el)[1]
            for el in self.df.columns
            if "D" in re.split(r"\.|-|_", el)
        ]
        # The check in [3, 5] was added to filter out the generation columns
        # by source in electricity
        NG_cols = [
            re.split(r"\.|-|_", el)[1]
            for el in self.df.columns
            if (
                ("NG" in re.split(r"\.|-|_", el))
                and (len(re.split(r"\.|-|_", el)) in [3, 5])
            )
        ]
        TI_cols = [
            re.split(r"\.|-|_", el)[1]
            for el in self.df.columns
            if "TI" in re.split(r"\.|-|_", el)
        ]
        ID_cols = [
            re.split(r"\.|-|_", el)[1]
            for el in self.df.columns
            if "ID" in re.split(r"\.|-|_", el)
        ]
        ID_cols2 = [
            re.split(r"\.|-|_", el)[2]
            for el in self.df.columns
            if "ID" in re.split(r"\.|-|_", el)
        ]

        self.D_cols = D_cols
        self.NG_cols = NG_cols
        self.TI_cols = TI_cols
        self.ID_cols = ID_cols
        self.ID_cols2 = ID_cols2

        if len(NG_cols) != len(D_cols):
            self.logger.warn("Inconsistent columns: len(NG_cols) != len(D_cols)")
        if set(NG_cols) != regions:
            self.logger.warn("Inconsistent columns: set(NG_cols) != regions")

        if "i" not in self.variable:
            if len(NG_cols) != len(TI_cols):
                self.logger.warn("Inconsistent columns: len(NG_cols) != len(TI_cols)")
            if set(NG_cols) != set(ID_cols):
                self.logger.warn("Inconsistent columns: set(NG_cols) != set(ID_cols)")
            if set(NG_cols) != set(ID_cols2):
                self.logger.warn("Inconsistent columns: set(NG_cols) != set(ID_cols2)")

        return sorted(list(regions))

    def get_trade_out(self, r=None):
        if r is None:
            r = self.regions
        if isinstance(r, str):
            r = [r]
        cols = []
        for ir2 in self.regions:
            cols += [self.KEY["ID"] % (ir, ir2) for ir in r]
        return [c for c in cols if c in self.df.columns]

    def checkBA(self, ba, tol=1e-2, log_level=logging.INFO):
        """
        Sanity check function
        TODO: add check for different generation sources, if data is present
        """
        logger = self.logger
        log_level_old = logger.level
        logger.setLevel(log_level)
        logger.debug("Checking %s" % ba)
        partners = self.get_trade_partners(ba)

        # NaNs
        for field in ["D", "NG", "TI"]:
            ind_na = self.df.loc[:, self.get_cols(r=ba, field=field)[0]].isna()
            cnt_na = ind_na.sum()
            if cnt_na != 0:
                logger.error(
                    "There are still %d nans for %s field %s" % (cnt_na, ba, field)
                )

        for ba2 in partners:
            cnt_na = self.df.loc[:, self.KEY["ID"] % (ba, ba2)].isna().sum()
            if cnt_na != 0:
                logger.error("There are still %d nans for %s-%s" % (cnt_na, ba, ba2))

        # TI+D == NG
        res1 = self.df.loc[:, self.get_cols(r=ba, field="NG")[0]] - (
            self.df.loc[:, self.get_cols(r=ba, field="D")[0]]
            + self.df.loc[:, self.get_cols(r=ba, field="TI")[0]]
        )
        if (res1.abs() > tol).sum() != 0:
            logger.error("%s: TI+D == NG violated" % ba)

        # TI == ID.sum()
        res2 = self.df.loc[:, self.get_cols(r=ba, field="TI")[0]] - self.df.loc[
            :, [self.KEY["ID"] % (ba, ba2) for ba2 in partners]
        ].sum(axis=1)
        if (res2.abs() > tol).sum() != 0:
            logger.error("%s: TI == ID.sum()violated" % ba)

        # ID[i,j] == -ID[j,i]
        for ba2 in partners:
            res3 = (
                self.df.loc[:, self.KEY["ID"] % (ba, ba2)]
                + self.df.loc[:, self.KEY["ID"] % (ba2, ba)]
            )
            if (res3.abs() > tol).sum() != 0:
                logger.error("%s-%s: ID[i,j] == -ID[j,i] violated" % (ba, ba2))

        # D and NG negative
        for field in ["D", "NG"]:
            ind_neg = self.df.loc[:, self.get_cols(r=ba, field=field)[0]] < 0
            cnt_neg = ind_neg.sum()
            if cnt_neg != 0:
                logger.error(
                    "%s: there are %d <0 values for field %s" % (ba, cnt_neg, field)
                )
        logger.setLevel(log_level_old)


################################################################################
# From gridemissions/clean.py
################################################################################

A = 1e4  # MWh
GAMMA = 10  # MWh
EPSILON = 1  # MWh


class BaDataCleaner(object):
    """
    Template class for data cleaning.

    This is mostly just a shell to show how cleaning classes should operate.
    """

    def __init__(self, ba_data):
        """
        Parameters
        ----------
        ba_data : BaData object
        """
        self.d = ba_data
        self.logger = logging.getLogger("clean")

    def process(self):
        pass


class BaDataBasicCleaner(BaDataCleaner):
    """
    Basic data cleaning class.

    We run this as the first step of the cleaning process.
    """

    def __init__(self, ba_data, add_ca_fuels: bool = True):
        """
        Parameters
        ----------
        ba_data : BaData object
        add_ca_fuels : bool, Add GEO and BIO columns for CISO generation?
        """
        self.add_ca_fuels = add_ca_fuels
        super().__init__(ba_data)

    def process(self):
        self.logger.info("Running BaDataBasicCleaner")
        start = time.time()
        data = self.d
        missing_D_cols = [col for col in data.NG_cols if col not in data.D_cols]
        self.logger.info("Adding demand columns for %d bas" % len(missing_D_cols))
        for ba in missing_D_cols:
            data.df.loc[:, data.KEY["D"] % ba] = 1.0
            data.df.loc[:, data.KEY["NG"] % ba] -= 1.0
            data.df.loc[:, data.KEY["TI"] % ba] -= 1.0

        # AVRN only exports to BPAT - this is missing for now
        if "AVRN" not in data.ID_cols:
            self.logger.info("Adding trade columns for AVRN")
            ba = "AVRN"
            ba2 = "BPAT"
            data.df.loc[:, data.KEY["ID"] % (ba, ba2)] = (
                data.df.loc[:, data.KEY["NG"] % ba] - 1.0
            )
            data.df.loc[:, data.KEY["ID"] % (ba2, ba)] = (
                -data.df.loc[:, data.KEY["NG"] % ba] + 1.0
            )

        if self.add_ca_fuels:
            # Add columns for biomass and geothermal for CISO
            # We are assuming constant generation for each of these sources
            # based on historical data. Before updating this, need to
            # contact the EIA API maintainers to understand why this isn't
            # reported and where to find it
            self.logger.info("Adding GEO and BIO columns for CISO")
            data.df.loc[:, "EBA.CISO-ALL.NG.GEO.H"] = 900.0
            data.df.loc[:, "EBA.CISO-ALL.NG.BIO.H"] = 600.0
        #           data.df.loc[:, "EBA.CISO-ALL.NG.H"] += 600.0 + 900.0

        # Add columns for the BAs that are outside of the US
        foreign_bas = list(
            set([col for col in data.ID_cols2 if col not in data.NG_cols])
        )
        self.logger.info(
            "Adding demand, generation and TI columns for %d foreign bas"
            % len(foreign_bas)
        )
        for ba in foreign_bas:
            trade_cols = [col for col in data.df.columns if "%s.ID.H" % ba in col]
            TI = -data.df.loc[:, trade_cols].sum(axis=1)
            data.df.loc[:, data.KEY["TI"] % ba] = TI
            exports = TI.apply(lambda x: max(x, 0))
            imports = TI.apply(lambda x: min(x, 0))
            data.df.loc[:, data.KEY["D"] % ba] = -imports
            data.df.loc[:, data.KEY["NG"] % ba] = exports
            if ba in ["BCHA", "HQT", "MHEB"]:
                # Assume for these Canadian BAs generation is hydro
                data.df.loc[:, data.KEY["SRC_WAT"] % ba] = exports
            else:
                # And all others are OTH (other)
                data.df.loc[:, data.KEY["SRC_OTH"] % ba] = exports
            for col in trade_cols:
                ba2 = re.split(r"\.|-|_", col)[1]
                data.df.loc[:, data.KEY["ID"] % (ba, ba2)] = -data.df.loc[:, col]

        # Make sure that trade columns exist both ways
        for col in data.get_cols(field="ID"):
            ba = re.split(r"\.|-|_", col)[1]
            ba2 = re.split(r"\.|-|_", col)[2]
            othercol = data.KEY["ID"] % (ba2, ba)
            if othercol not in data.df.columns:
                self.logger.info("Adding %s" % othercol)
                data.df.loc[:, othercol] = -data.df.loc[:, col]

        # Filter unrealistic values using self.reject_dict
        self._create_reject_dict()
        cols = (
            data.get_cols(field="D")
            + data.get_cols(field="NG")
            + data.get_cols(field="TI")
            + data.get_cols(field="ID")
        )
        for col in cols:
            s = data.df.loc[:, col]
            data.df.loc[:, col] = s.where(
                (s >= self.reject_dict[col][0]) & (s <= self.reject_dict[col][1])
            )

        # Do the same for the generation by source columns
        # If there is no generation by source, add one that is OTH
        # Edge case for solar:
        # There are a lot of values at -50 MWh or so during the night. We want
        # to set those to 0, but consider that very negative values (below
        # -1GW) are rejected
        for ba in data.regions:
            missing = True
            for src in SRC:
                col = data.KEY["SRC_%s" % src] % ba
                if col in data.df.columns:
                    missing = False
                    s = data.df.loc[:, col]
                    if src == "SUN":
                        self.reject_dict[col] = (-1e3, 200e3)
                    data.df.loc[:, col] = s.where(
                        (s >= self.reject_dict[col][0])
                        & (s <= self.reject_dict[col][1])
                    )
                    if src == "SUN":
                        data.df.loc[:, col] = data.df.loc[:, col].apply(
                            lambda x: max(x, 0)
                        )
            if missing:
                data.df.loc[:, data.KEY["SRC_OTH"] % ba] = data.df.loc[
                    :, data.KEY["NG"] % ba
                ]

        # Reinitialize fields
        self.logger.info("Reinitializing fields")
        data = BaData(df=data.df)

        self.r = data

        self.logger.info("Basic cleaning took %.2f seconds" % (time.time() - start))

    def _create_reject_dict(self):
        """
        Create a defaultdict to store ranges outside of which values are
        considered unrealistic.

        The default range is (-1., 200e3) MW. Manual ranges can be set for
        specific columns here if that range is not strict enough.
        """
        reject_dict = defaultdict(lambda: (-1.0, 200e3))
        for col in self.d.get_cols(field="TI"):
            reject_dict[col] = (-100e3, 100e3)
        for col in self.d.get_cols(field="ID"):
            reject_dict[col] = (-100e3, 100e3)
        reject_dict["EBA.AZPS-ALL.D.H"] = (1.0, 30e3)
        reject_dict["EBA.BANC-ALL.D.H"] = (1.0, 6.5e3)
        reject_dict["EBA.BANC-ALL.TI.H"] = (-5e3, 5e3)
        reject_dict["EBA.CISO-ALL.NG.H"] = (5e3, 60e3)
        self.reject_dict = reject_dict


def rolling_window_filter(
    df,
    offset=20 * 24,
    min_periods=100,
    center=True,
    replace_nan_with_mean=True,
    return_mean=False,
):
    """
    Apply a rolling window filter to a dataframe.

    Filter using dynamic bounds: reject points that are farther than 4 standard
    deviations from the mean, using a rolling window to compute the mean and
    standard deviation.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to filter
    offset : int
        Passed on to pandas' rolling function
    min_periods : int
        Passed on to pandas' rolling function
    center : bool
        Passed on to pandas' rolling function
    replace_nan_with_mean : bool
        Whether to replace NaNs with the mean of the timeseries at the end of
        the procedure

    Notes
    -----
    Keeps at least 200 MWh around the mean as an acceptance range.
    """
    for col in df.columns:
        rolling_ = df[col].rolling(offset, min_periods=min_periods, center=center)
        median_ = rolling_.median()
        # 84.1 % quantile = mean + 1 std in normal data
        iqr_ = (rolling_.quantile(0.95) - rolling_.quantile(0.05)).apply(
            lambda x: max(150, x)
        )
        ub = median_ + 3 * iqr_
        lb = median_ - 3 * iqr_
        idx_reject = (df[col] >= ub) | (df[col] <= lb)
        df.loc[idx_reject, col] = np.nan
        if replace_nan_with_mean:
            # First try interpolating linearly, but only for up to 3 hours
            df.loc[:, col] = df.loc[:, col].interpolate(limit=3)
            # If there is more than 3 hours of missing data, use rolling mean
            df.loc[df[col].isnull(), col] = median_.loc[df[col].isnull()]

    if return_mean:
        mean_ = df.rolling(offset, min_periods=min_periods, center=center).mean()
        return (df, mean_)
    return df


class BaDataRollingCleaner(BaDataCleaner):
    """
    Rolling window cleaning.

    This applies the `rolling_window_filter` function to the dataset. In order
    to apply this properly to the beginning of the dataset, we load past data
    that will be used for the cleaning - that is then dropped.
    """

    def process(self, file_name="", folder_hist="", nruns=2):
        """
        Processor function for the cleaner object.

        Parameters
        ----------
        file_name : str
            Base name of the file from which to read historical data.
            Data is read from "%s_basic.csv" % file_name
        folder_hist : str
            Folder from which to read historical data
        nruns : int
            Number of times to apply the rolling window procedure

        Notes
        -----
        If we are not processing a large amount of data at a time, we may not
        have enough data to appropriately estimate the rolling mean and
        standard deviation for the rolling window procedure. If values are
        given for `file_name` and `folder_hist`, data will be read from a
        historical dataset to estimate the rolling mean and standard deviation.
        If there are very large outliers, they can 'mask' smaller outliers.
        Running the rolling window procedure a couple of times helps with this
        issue.
        """
        self.logger.info("Running BaDataRollingCleaner (%d runs)" % nruns)
        start = time.time()
        data = self.d

        # Remember what part we are cleaning
        idx_cleaning = data.df.index

        try:
            # Load the data we already have in memory
            df_hist = pd.read_csv(
                os.path.join(folder_hist, "%s_basic.csv" % file_name),
                index_col=0,
                parse_dates=True,
            )

            # Only take the last 1,000 rows
            # Note that if df_hist has less than 1,000 rows,
            # pandas knows to select df_hist without raising an error.
            df_hist = df_hist.iloc[-1000:]

            # Overwrite with the new data
            old_rows = df_hist.index.difference(data.df.index)
            df_hist = pd.concat([data.df, df_hist.loc[old_rows, :]])
            df_hist.sort_index(inplace=True)

        except FileNotFoundError:
            self.logger.info("No history file")
            df_hist = data.df

        # Apply rolling horizon threshold procedure
        # 20200206 update: don't try replacing NaNs anymore, leave that to the
        # next step
        for _ in range(nruns):
            df_hist = rolling_window_filter(df_hist, replace_nan_with_mean=False)

        # Deal with NaNs
        # First deal with NaNs by taking the average of the previous day and
        # next day. In general we observe strong daily patterns so this seems
        # to work well. Limit the filling procedure to one day at a time. If
        # there are multiple missing days, this makes for a smoother transition
        # between the two valid days. If we had to do this more than 4 times,
        # give up and use forward and backward fills without limits
        for col in df_hist.columns:
            npasses = 0
            while (df_hist.loc[:, col].isna().sum() > 0) and (npasses < 4):
                npasses += 1
                df_hist.loc[:, col] = pd.concat(
                    [
                        df_hist.loc[:, col].groupby(df_hist.index.hour).ffill(limit=1),
                        df_hist.loc[:, col].groupby(df_hist.index.hour).bfill(limit=1),
                    ],
                    axis=1,
                ).mean(axis=1)
            if npasses == 4:
                self.logger.debug("A lot of bad data for %s" % col)
                df_hist.loc[:, col] = pd.concat(
                    [
                        df_hist.loc[:, col].groupby(df_hist.index.hour).ffill(),
                        df_hist.loc[:, col].groupby(df_hist.index.hour).bfill(),
                    ],
                    axis=1,
                ).mean(axis=1)

            # All bad data columns
            if df_hist.loc[:, col].isna().sum() == len(df_hist):
                df_hist.loc[:, col] = 0.0

        # Some NaNs will still remain - try using the rolling mean average
        df_hist, mean_ = rolling_window_filter(
            df_hist, replace_nan_with_mean=True, return_mean=True
        )

        if df_hist.isna().sum().sum() > 0:
            self.logger.warning("There are still some NaNs. Unexpected")

        # Just keep the indices we are working on currently
        data = BaData(df=df_hist.loc[idx_cleaning, :])

        self.r = data
        self.weights = mean_.loc[idx_cleaning, :].applymap(
            lambda x: A / max(GAMMA, abs(x))
        )

        self.logger.info(
            "Rolling window cleaning took %.2f seconds" % (time.time() - start)
        )


class BaDataCvxCleaner(BaDataCleaner):
    """
    Optimization-based cleaning class.

    Uses cvxpy.
    """

    def __init__(self, ba_data, weights=None):
        super().__init__(ba_data)
        self.weights = weights
        if weights is not None:
            self.d.df = pd.concat(
                [self.d.df, weights.rename(lambda x: x + "_W", axis=1)], axis=1
            )

    def process(self, debug=False, with_ng_src=True):
        start = time.time()
        self.logger.info("Running BaDataCvxCleaner for %d rows" % len(self.d.df))
        self.d.df = self.d.df.fillna(0)

        results = []

        def cvx_solve(row, regions, debug=False):
            if row.isna().sum() > 0:
                raise ValueError("Cannot call this method on data with NaNs")

            n_regions = len(regions)

            D = row[[KEYS["E"]["D"] % r for r in regions]].values
            D_W = [
                el**0.5
                for el in row[[KEYS["E"]["D"] % r + "_W" for r in regions]].values
            ]
            NG = row[[KEYS["E"]["NG"] % r for r in regions]].values
            NG_W = [
                el**0.5
                for el in row[[KEYS["E"]["NG"] % r + "_W" for r in regions]].values
            ]
            TI = row[[KEYS["E"]["TI"] % r for r in regions]].values
            TI_W = [
                el**0.5
                for el in row[[KEYS["E"]["TI"] % r + "_W" for r in regions]].values
            ]

            delta_D = cp.Variable(n_regions, name="delta_D")
            delta_NG = cp.Variable(n_regions, name="delta_NG")
            delta_TI = cp.Variable(n_regions, name="delta_TI")

            obj = (
                cp.sum_squares(cp.multiply(D_W, delta_D))
                + cp.sum_squares(cp.multiply(NG_W, delta_NG))
                + cp.sum_squares(cp.multiply(TI_W, delta_TI))
            )

            ID = {}
            ID_W = {}
            for i, ri in enumerate(regions):
                for j, rj in enumerate(regions):
                    if KEYS["E"]["ID"] % (ri, rj) in row.index:
                        ID[(ri, rj)] = row[KEYS["E"]["ID"] % (ri, rj)]
                        ID_W[(ri, rj)] = row[KEYS["E"]["ID"] % (ri, rj) + "_W"]
            delta_ID = {k: cp.Variable(name=f"{k}") for k in ID}
            constraints = [
                D + delta_D >= 1.0,
                NG + delta_NG >= 1.0,
                D + delta_D + TI + delta_TI - NG - delta_NG == 0.0,
            ]

            if with_ng_src:
                NG_SRC = {}
                NG_SRC_W = {}

                for i, src in enumerate(SRC):
                    for j, r in enumerate(regions):
                        if KEYS["E"][f"SRC_{src}"] % r in row.index:
                            NG_SRC[(src, r)] = row[KEYS["E"][f"SRC_{src}"] % r]
                            NG_SRC_W[(src, r)] = row[KEYS["E"][f"SRC_{src}"] % r + "_W"]
                delta_NG_SRC = {k: cp.Variable(name=f"{k}") for k in NG_SRC}

                for k in NG_SRC:
                    constraints += [NG_SRC[k] + delta_NG_SRC[k] >= 1.0]
                    obj += NG_SRC_W[k] * delta_NG_SRC[k] ** 2

            # Add the antisymmetry constraints twice is less efficient but not a huge deal.
            for ri, rj in ID:  # then (rj, ri) must also be in ID
                constraints += [
                    ID[(ri, rj)]
                    + delta_ID[(ri, rj)]
                    + ID[(rj, ri)]
                    + delta_ID[(rj, ri)]
                    == 0.0
                ]
                obj += ID_W[(ri, rj)] * delta_ID[(ri, rj)] ** 2

            for i, ri in enumerate(regions):
                if with_ng_src:
                    constraints += [
                        NG[i]
                        + delta_NG[i]
                        - cp.sum(
                            [
                                NG_SRC[(src, ri)] + delta_NG_SRC[(src, ri)]
                                for src in SRC
                                if (src, ri) in NG_SRC
                            ]
                        )
                        == 0.0
                    ]
                constraints += [
                    TI[i]
                    + delta_TI[i]
                    - cp.sum(
                        [
                            ID[(ri, rj)] + delta_ID[(ri, rj)]
                            for rj in regions
                            if (ri, rj) in ID
                        ]
                    )
                    == 0.0
                ]
            objective = cp.Minimize(obj)

            prob = cp.Problem(objective, constraints)
            prob.solve()

            if with_ng_src:
                r = pd.concat(
                    [
                        pd.Series(
                            NG + delta_NG.value,
                            index=[KEYS["E"]["NG"] % r for r in regions],
                        ),
                        pd.Series(
                            D + delta_D.value,
                            index=[KEYS["E"]["D"] % r for r in regions],
                        ),
                        pd.Series(
                            TI + delta_TI.value,
                            index=[KEYS["E"]["TI"] % r for r in regions],
                        ),
                        pd.Series(
                            {KEYS["E"]["ID"] % k: ID[k] + delta_ID[k].value for k in ID}
                        ),
                        pd.Series(
                            {
                                KEYS["E"][f"SRC_{s}"] % r: NG_SRC[(s, r)]
                                + delta_NG_SRC[(s, r)].value
                                for (s, r) in NG_SRC
                            }
                        ),
                        pd.Series({"CleaningObjective": prob.value}),
                    ]
                )
            else:
                r = pd.concat(
                    [
                        pd.Series(
                            NG + delta_NG.value,
                            index=[KEYS["E"]["NG"] % r for r in regions],
                        ),
                        pd.Series(
                            D + delta_D.value,
                            index=[KEYS["E"]["D"] % r for r in regions],
                        ),
                        pd.Series(
                            TI + delta_TI.value,
                            index=[KEYS["E"]["TI"] % r for r in regions],
                        ),
                        pd.Series(
                            {KEYS["E"]["ID"] % k: ID[k] + delta_ID[k].value for k in ID}
                        ),
                        pd.Series({"CleaningObjective": prob.value}),
                    ]
                )

            if not debug:
                return r

            if with_ng_src:
                deltas = pd.concat(
                    [
                        pd.Series(
                            delta_NG.value, index=[KEYS["E"]["NG"] % r for r in regions]
                        ),
                        pd.Series(
                            delta_D.value, index=[KEYS["E"]["D"] % r for r in regions]
                        ),
                        pd.Series(
                            delta_TI.value, index=[KEYS["E"]["TI"] % r for r in regions]
                        ),
                        pd.Series({KEYS["E"]["ID"] % k: delta_ID[k].value for k in ID}),
                        pd.Series(
                            {
                                KEYS["E"][f"SRC_{s}"] % r: delta_NG_SRC[(s, r)].value
                                for (s, r) in NG_SRC
                            }
                        ),
                    ]
                )
            else:
                deltas = pd.concat(
                    [
                        pd.Series(
                            delta_NG.value, index=[KEYS["E"]["NG"] % r for r in regions]
                        ),
                        pd.Series(
                            delta_D.value, index=[KEYS["E"]["D"] % r for r in regions]
                        ),
                        pd.Series(
                            delta_TI.value, index=[KEYS["E"]["TI"] % r for r in regions]
                        ),
                        pd.Series({KEYS["E"]["ID"] % k: delta_ID[k].value for k in ID}),
                    ]
                )
            return pd.concat([r, deltas.rename(lambda x: x + "_Delta")])

        cvx_solve = dask.delayed(cvx_solve)
        for idx, row in self.d.df.iterrows():
            results.append(cvx_solve(row, self.d.regions, debug=debug))
        results = dask.compute(*results, scheduler="processes")
        df = pd.DataFrame(results, index=self.d.df.index)

        self.r = df.loc[
            :,
            [
                c
                for c in df.columns
                if "Delta" not in c and "CleaningObjective" not in c
            ],
        ]
        self.CleaningObjective = df.CleaningObjective
        self.deltas = df.loc[:, [c for c in df.columns if "Delta" in c]]

        # Make sure the cleaning step performed as expected
        self.r = BaData(df=self.r)
        self.logger.info("Checking BAs...")
        for ba in self.r.regions:
            self.r.checkBA(ba)
        self.logger.info("Execution took %.2f seconds" % (time.time() - start))


################################################################################
# From gridemissions/workflows.py
################################################################################

# Optimization-based cleaning is different pre and post July 2018
THRESH_DATE = pd.to_datetime("20180701", utc=True)


def make_dataset(
    file_name="EBA",
    tmp_folder="tmp",
    folder_hist=None,
    add_ca_fuels: bool = True,
):
    """
    Run data through the cleaning workflow.

    Reads "%s_raw.csv" % file_name from `tmp_folder`, which must already exist.

    Uses historical data if available.

    Option to add CISO fuels not reported in 930 (GEO, BIO)
    """
    start_time = time.time()

    logger = logging.getLogger("scraper")

    # Basic data cleaning
    logger.info("Basic data cleaning")
    data = BaData(fileNm=join(tmp_folder, "%s_raw.csv" % file_name))
    cleaner = BaDataBasicCleaner(data, add_ca_fuels=add_ca_fuels)
    cleaner.process()
    cleaner.r.df.to_csv(join(tmp_folder, "%s_basic.csv" % file_name))
    data = cleaner.r

    weights = None
    if folder_hist is not None:  # Rolling-window-based data cleaning
        logger.info("Rolling window data cleaning")
        data = BaData(fileNm=join(tmp_folder, "%s_basic.csv" % file_name))
        cleaner = BaDataRollingCleaner(data)
        cleaner.process(file_name, folder_hist)
        cleaner.r.df.to_csv(join(tmp_folder, "%s_rolling.csv" % file_name))
        data = cleaner.r
        weights = cleaner.weights
        weights.to_csv(join(tmp_folder, "%s_weights.csv" % file_name))
    else:
        logger.warning("No rolling window data cleaning!")

    if len(data.df.loc[:THRESH_DATE, :]) > 0:
        logger.info(f"Optimization-based cleaning without src data: pre {THRESH_DATE}")
        ba_data = BaData(df=data.df.loc[:THRESH_DATE, :])
        if weights is not None:
            cleaner = BaDataCvxCleaner(ba_data, weights=weights.loc[:THRESH_DATE, :])
        else:
            cleaner = BaDataCvxCleaner(ba_data)
        cleaner.process(debug=False, with_ng_src=False)
        cleaner.r.df.to_csv(join(tmp_folder, "%s_opt_no_src.csv" % file_name))
        cleaner.CleaningObjective.to_csv(
            join(tmp_folder, "%s_objective_no_src.csv" % file_name)
        )

    # Only keep going if we have data post THRESH_DATE
    if len(data.df.loc[THRESH_DATE:, :]) == 0:
        return

    logger.info(f"Optimization-based cleaning with src data: post {THRESH_DATE}")
    data.df = data.df.loc[THRESH_DATE:, :]
    if weights is not None:
        cleaner = BaDataCvxCleaner(data, weights=weights.loc[THRESH_DATE:, :])
    else:
        cleaner = BaDataCvxCleaner(data)
    cleaner.process(debug=False)
    cleaner.r.df.to_csv(join(tmp_folder, "%s_opt.csv" % file_name))
    cleaner.CleaningObjective.to_csv(join(tmp_folder, "%s_objective.csv" % file_name))

    # Post-processing (none for now)
    cleaner.r.df.to_csv(join(tmp_folder, "%s_elec.csv" % file_name))
    data = cleaner.r

    logger.info(
        "oge.physics_reconciliation.make_dataset took %.2f seconds"
        % (time.time() - start_time)
    )
