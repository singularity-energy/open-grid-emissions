from typing import Any, Optional

import numpy as np
import pandas as pd

from oge.logging_util import get_logger

logger = get_logger(__name__)


class AnomalyScreeningFirstStep:
    """Screen timeseries for anomalous value following the algorithm steps described
    in Tyler H. Ruggles et al. Developing reliable hourly electricity demand data
    through screening and imputation (2020)

    Note that the screening algorithms have been developed for demand timeseries and
    some of these algorithms might not be tailored to other time series.

    The screening is conducted in 2 steps. Step 1 removes the most egregious anomalies
    where few or no calculations are needed. Afterward, in Step 2, the most extreme
    values have been removed making calculations of local characteristics of the data
    more reasonable. Through this screening process hourly values can be recategorized
    from okay to other classifications based on the algorithms.
    """

    def __init__(
        self, df: pd.DataFrame, field: str, **kwargs: Optional[dict[str, Any]]
    ):
        """Constructor

        Arguments
        ---------
        * `df` : timeseries as a data frame.
        * `field` : Column to screen.
        * `kwargs` : Optional keyword arguments.
        """
        self.df = df.copy().reset_index()
        self.field = field
        self.global_cut_multiplier = (
            kwargs["global_cut_multiplier"] if "global_cut_multiplier" in kwargs else 10
        )

        self.add_category()

    def add_category(self):
        """Add category label to track which algorithm screens an hourly value"""
        self.df["category"] = np.where(self.df[self.field].isna(), "MISSING", "OKAY")
        self.df = self.df.assign(missing=self.df[self.field].isna())

    def flag_negative_values(self):
        """Flag negative values"""
        flag = self.df[self.field] < 0
        self.df["category"] = np.where(flag, "NEGATIVE", self.df["category"])
        filtered = np.where(flag, self.df[self.field], np.nan)

        self.df["negative_filtered"] = filtered
        self.df[self.field] = self.df[self.field].mask(flag)

    def flag_zero_values(self):
        """Flag zero values"""
        flag = self.df[self.field] == 0
        self.df["category"] = np.where(flag, "ZERO", self.df["category"])
        filtered = np.where(flag, self.df[self.field], np.nan)

        self.df["zero_filtered"] = filtered
        self.df[self.field] = self.df[self.field].mask(flag)

    def flag_runs(self):
        """flag the third hour onwards as anomalous in a series of identical values"""

        d1 = self.df[self.field].diff(periods=1)
        d2 = self.df[self.field].diff(periods=2)

        filtered = self.df[self.field].mask((d1 == 0) & (d2 == 0))
        self.df["run_filtered"] = np.where(
            self.df[self.field] != filtered, self.df[self.field], np.nan
        )
        self.df[self.field] = filtered
        self.df["category"] = np.where(
            self.df["run_filtered"].notna(), "IDENTICAL_RUN", self.df["category"]
        )

    def flag_global_extreme_values(self):
        """Flag entries value is above a multiplier of the median."""
        med = np.nanmedian(self.df[self.field])
        filtered = self.df[self.field].where(
            self.df[self.field] < med * self.global_cut_multiplier
        )
        self.df["global_extreme_filtered"] = np.where(
            self.df[self.field] != filtered, self.df[self.field], np.nan
        )
        self.df["category"] = self.df["category"].mask(
            ((self.df[self.field] != filtered) & (self.df[self.field].notna())),
            other="GLOBAL_EXTREME",
        )
        self.df[self.field] = filtered

    def flag_global_extreme_plus_minus_one_hour(self):
        """Flag entries +/- 1 hour from any extreme flagged hours"""
        filtered = [np.nan for _ in self.df.index]
        for i in self.df.index:
            if self.df.loc[i, "category"] == "GLOBAL_EXTREME":
                if i != 0 and self.df.loc[i - 1, "category"] == "OKAY":
                    self.df.loc[i - 1, "category"] = "GLOBAL_EXTREME_+/-_1H"
                    filtered[i - 1] = self.df.loc[i - 1, self.field]
                    self.df.loc[i - 1, self.field] = np.nan
                if i != len(self.df) - 1 and self.df.loc[i + 1, "category"] == "OKAY":
                    self.df.loc[i + 1, "category"] = "GLOBAL_EXTREME_+/-_1H"
                    filtered[i + 1] = self.df.loc[i + 1, self.field]
                    self.df.loc[i + 1, self.field] = np.nan
        self.df["global_extreme_filtered_+/-_1h"] = filtered

    def run_all_first_step(self):
        """Run all four algorithms"""
        if len(self.df.query("category != 'MISSING'")) > 0:
            self.flag_negative_values()
            self.flag_zero_values()
            self.flag_runs()
            self.flag_global_extreme_values()
            self.flag_global_extreme_plus_minus_one_hour()

    def summary(self):
        """Print for each screening steps the number of entries flagged"""
        summary = self.df.groupby("category").size().to_frame(name="count")
        logger.info(
            summary.assign(fraction=summary["count"] / len(self.df)).to_string()
        )

    def get_filtered_df(self):
        """Get filtered data frame"""
        return self.df


class AnomalyScreeningSecondStep(AnomalyScreeningFirstStep):
    """This class encloses the methods for performimg the second set of screening
    described in Tyler H. Ruggles et al. Developing reliable hourly electricity demand
    data through screening and imputation (2020)
    """

    def __init__(
        self, df: pd.DataFrame, field: str, **kwargs: Optional[dict[str, Any]]
    ):
        super().__init__(df, field)

        self.run_all_first_step(
            kwargs["global_cut_multiplier"]
        ) if "global_cut_multiplier" in kwargs else self.run_all_first_step()
        self.short_hour_window = (
            kwargs["short_hour_window"] if "short_hour_window" in kwargs else 48
        )
        self.iqr_hour_window = (
            kwargs["iqr_hour_window"] if "iqr_hour_window" in kwargs else 240
        )
        self.n_day = kwargs["n_day"] if "n_day" in kwargs else 10
        self.delta_multiplier = (
            kwargs["delta_multiplier"] if "delta_multiplier" in kwargs else 2
        )

        self.prepare_second_step_of_screening()

    def prepare_second_step_of_screening(self):
        """Calculate characteristics needed for step 2 screening"""
        self.calculate_rolling_median_short()
        self.calculate_rolling_median_long()
        self.subtract_rolling_median_short()
        self.subtract_rolling_iqr()
        self.calculate_delta()
        self.calculate_rolling_delta_iqr()
        self.calculate_hourly_median_deviation()
        self.calculate_hourly_relative_difference()
        self.calculate_hourly_relative_delta()

    def flag_local_deviation_from_expected_value(
        self, multiplier_up: float = 3.5, multiplier_down: float = 2.5
    ):
        """Construct an estimated value based on 48 hour moving median and the typical
        local daily cycle. Flag the hour for which the value deviates significantly
        from the estimate. Here, significance is defined with respect to the local
        interquartile range (IQR) of the difference between the value and the 48 hour
        moving median (IQR of the spread of the daily cycle).

        Arguments
        ---------
            * `multiplier_up`: Upwards threshold for cut. Defaults to 3.5.
            * `multiplier_down`: Downwards threshold for cut. Defaults to 2.5.
        """
        filtered = self.df[self.field].where(
            (
                self.df[self.field]
                < self.df["rolling_median_short"] * self.df["hourly_median_dev"]
                + multiplier_up * self.df[f"{self.field}_minus_rolling_IQR"]
            )
        )
        self.df[f"local_{self.field}_filtered_up"] = np.where(
            self.df[self.field] != filtered, self.df[self.field], np.nan
        )
        self.df["category"] = self.df["category"].mask(
            ((self.df[self.field] != filtered) & (self.df[self.field].notna())),
            other="LOCAL_UP",
        )
        self.df[self.field] = filtered

        filtered = self.df[self.field].where(
            (
                self.df[self.field]
                > self.df["rolling_median_short"] * self.df["hourly_median_dev"]
                - multiplier_down * self.df[f"{self.field}_minus_rolling_IQR"]
            )
        )
        self.df[f"local_{self.field}_filtered_down"] = np.where(
            self.df[self.field] != filtered, self.df[self.field], np.nan
        )
        self.df["category"] = self.df["category"].mask(
            ((self.df[self.field] != filtered) & (self.df[self.field].notna())),
            other="LOCAL_DOWN",
        )
        self.df[self.field] = filtered

    def flag_deltas(self):
        """Flag on a multiplier of the IQR. Hours with large deltas on both sides
        are considered"""
        filtered = self.df[self.field].mask(
            (
                (
                    self.df["delta_pre"]
                    > self.df["delta_rolling_IQR"] * self.delta_multiplier
                )
                & (
                    self.df["delta_post"]
                    > self.df["delta_rolling_IQR"] * self.delta_multiplier
                )
            )
            | (
                (
                    self.df["delta_pre"]
                    < -1.0 * self.df["delta_rolling_IQR"] * self.delta_multiplier
                )
                & (
                    self.df["delta_post"]
                    < -1.0 * self.df["delta_rolling_IQR"] * self.delta_multiplier
                )
            )
        )

        self.df["double_delta_filtered"] = np.where(
            self.df[self.field] != filtered, self.df[self.field], np.nan
        )
        self.df["category"] = self.df["category"].mask(
            ((self.df[self.field] != filtered) & (self.df[self.field].notna())),
            other="DOUBLE_DELTA",
        )
        self.df[self.field] = filtered

    def calculate_rolling_median_short(self):
        """Calculate median in short rolling windows."""
        self.df["rolling_median_short"] = (
            self.df[self.field]
            .rolling(self.short_hour_window, min_periods=1, center=True)
            .median()
        )

    def calculate_rolling_median_long(self):
        """Calculate median long rolling windows."""
        self.df["rolling_median_long"] = (
            self.df[self.field]
            .rolling(48 * self.n_day, min_periods=1, center=True)
            .median()
        )

    def subtract_rolling_median_short(self):
        """Subtract median calculated in the short rolling windows."""
        self.df.insert(
            len(self.df.columns),
            f"{self.field}_minus_rolling_median_short",
            self.df[self.field] - self.df["rolling_median_short"],
        )

    def subtract_rolling_iqr(self):
        """Subtract IQR in rolling windows."""
        rolling_window = self.df[f"{self.field}_minus_rolling_median_short"].rolling(
            self.iqr_hour_window, min_periods=1, center=True
        )
        self.df[f"{self.field}_minus_rolling_IQR"] = rolling_window.quantile(
            0.75
        ) - rolling_window.quantile(0.25)

    def calculate_delta(self):
        """Calculate deltas."""
        self.df = self.df.assign(delta_pre=self.df[self.field].diff())
        self.df = self.df.assign(delta_post=self.df[self.field].diff(periods=-1))

    def calculate_rolling_delta_iqr(self):
        """Calculate IQR values in rolling delta windows."""
        rolling_window = self.df["delta_pre"].rolling(
            self.iqr_hour_window, min_periods=1, center=True
        )
        self.df["delta_rolling_IQR"] = rolling_window.quantile(
            0.75
        ) - rolling_window.quantile(0.25)

    def calculate_hourly_median_deviation(self):
        """Calculate deviation to median value"""
        value_minus_rolling_median_short = self.df[
            f"{self.field}_minus_rolling_median_short"
        ]
        for i in range(-self.n_day, self.n_day + 1):
            # Already initialized with zero value
            if i == 0:
                continue
            value_minus_rolling_median_short = pd.concat(
                [
                    value_minus_rolling_median_short,
                    self.df.shift(periods=i * 24)[
                        f"{self.field}_minus_rolling_median_short"
                    ],
                ],
                axis=1,
            )

        self.df[f"{self.field}_minus_rolling_median_short"] = (
            value_minus_rolling_median_short.median(axis=1, skipna=True)
        )
        self.df = self.df.assign(
            hourly_median_dev=1.0
            + self.df[f"{self.field}_minus_rolling_median_short"]
            / self.df["rolling_median_long"]
        )

    def calculate_hourly_relative_difference(self):
        """Calculate hourly relative difference wrt expected value."""
        self.df = self.df.assign(
            hourly_relative_difference_short=self.df[self.field]
            / (self.df["rolling_median_short"] * self.df["hourly_median_dev"])
        )
        self.df = self.df.assign(
            hourly_relative_difference_long=self.df[self.field]
            / (self.df["rolling_median_long"] * self.df["hourly_median_dev"])
        )

    def calculate_hourly_relative_delta(self):
        """Calculate hourly relative delats wrt expected value."""
        self.df = self.df.assign(
            hourly_relative_delta_pre=self.df["hourly_relative_difference_short"].diff()
        )
        self.df = self.df.assign(
            hourly_relative_delta_post=self.df["hourly_relative_difference_short"].diff(
                periods=-1
            )
        )

    def calculate_relative_difference_iqr(self):
        return np.nanpercentile(
            self.df["hourly_relative_delta_pre"], 75
        ) - np.nanpercentile(self.df["hourly_relative_delta_pre"], 25)
