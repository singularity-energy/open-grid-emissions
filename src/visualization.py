"""
Helper functions for visualization

"""
import numpy as np
import pandas as pd

"""
From a timeseries, make a day-hour array with hours in EST. Also return day labels for plotting.

Assumes timeseries is hourly, because we don't know how to resample to hourly if it's not.
"""


def day_hour_heatmap(timeseries: pd.Series, year: int = 2021):
    timeseries = timeseries.rename("a")  # needs to be named so we can merge

    # Make a new dataframe that has all full days in year in EST
    timeseries.index = timeseries.index.tz_convert("EST")
    hours_index = pd.DataFrame(
        index=pd.date_range(
            f"{year}-01-01 T00:00", f"{year}-12-31 T23:00", freq="H"
        ).tz_localize("EST")
    )
    hours_index = hours_index.merge(
        timeseries, how="left", left_index=True, right_index=True
    )

    # Get times for x labels
    times = hours_index.resample("D").first().index.to_numpy()
    hours = hours_index.to_numpy()
    hours = hours.reshape((24, len(hours) // 24), order="F")

    return (hours, times)
