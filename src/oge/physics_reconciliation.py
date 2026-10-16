"""Physics-based reconciliation of EIA-930 data.

EIA-930 data are not always physically consistent: a BA's demand plus its total
interchange does not always equal its net generation, total interchange does not
always equal the sum of interchange with each neighbor, and two neighboring BAs do not
always agree on the interchange between them. This module cleans the data in three
steps:

1. Basic cleaning: drop values outside of realistic bounds and add series for BAs
   outside of the US.
2. Rolling window cleaning: drop anomalous values using a rolling median and
   percentile range, then fill gaps.
3. Reconciliation: for each hour, find the smallest weighted adjustments to the data
   that make it physically consistent, using a quadratic program.

The data enter and leave this module as three tables shaped like PUDL's EIA-930
tables (see `EIA930Tables`). Internally, each step operates on a "series" data frame
with one row per hour and one column per time series, where each column is identified
by a (field, balancing_authority_code_eia, other) tuple. `other` is the adjacent BA
for interchange series, the energy source for generation by source series, and "" for
the other series.

This code was ported from the `gridemissions` package
(https://github.com/jdechalendar/gridemissions), by way of the Singularity fork
(https://github.com/singularity-energy/gridemissions, commit
0cbe7620620f02ab91227940ef107b53bd1f9bb4). The cleaning logic is unchanged from that
version, but the reconciliation is solved with a more precise solver (see
`HourlyReconciliationProblem`), so reconciled values can differ from gridemissions' by
up to about 1 MWh. The original code is distributed under the following license:

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
import time
from collections import defaultdict
from typing import NamedTuple

import cvxpy as cp
import numpy as np
import pandas as pd
import scipy.sparse

from oge.logging_util import get_logger

logger = get_logger(__name__)

# fields that identify the type of each time series
DEMAND = "demand"
NET_GENERATION = "net_generation"
TOTAL_INTERCHANGE = "total_interchange"
INTERCHANGE = "interchange"
NET_GENERATION_BY_ENERGY_SOURCE = "net_generation_by_energy_source"

# names of the column levels of a series data frame
SERIES_LEVELS = ["field", "balancing_authority_code_eia", "other"]

# energy sources reported in EIA-930, in the order of their series in the
# reconciliation
ENERGY_SOURCES = [
    "coal",
    "gas",
    "nuclear",
    "oil",
    "other",
    "solar",
    "unknown",
    "hydro",
    "wind",
]

# Canadian BAs whose generation is assumed to be hydro when it is estimated from
# their interchange with US BAs
HYDRO_FOREIGN_BAS = ["BCHA", "HQT", "MHEB"]

# the optimization only reconciles generation by energy source after this date,
# because EIA-930 did not collect generation by energy source before it
GENERATION_BY_ENERGY_SOURCE_START_DATE = pd.to_datetime("20180701", utc=True)

# parameters used to compute the weight of each series in the optimization
WEIGHT_NUMERATOR = 1e4  # MWh
WEIGHT_MINIMUM_DENOMINATOR = 10  # MWh

# lowest value allowed for demand, net generation and generation by energy source in
# the optimization
OPTIMIZATION_LOWER_BOUND = 1.0  # MWh

# settings of the Clarabel solver. With these tolerances, the reconciled values are
# within about 0.01 MWh of the exact solution.
CLARABEL_SETTINGS = {"tol_gap_abs": 1e-10, "tol_gap_rel": 1e-10, "tol_feas": 1e-10}

# names of the parquet files written for each table of `EIA930Tables`
TABLE_FILE_NAMES = {
    "operations": "operations",
    "net_generation_by_energy_source": "net_generation_by_energy_source",
    "interchange": "interchange",
}


class EIA930Tables(NamedTuple):
    """EIA-930 data, in three tables shaped like PUDL's EIA-930 tables.

    Attributes:
        operations (pd.DataFrame): one row per BA and hour, with columns
            "datetime_utc", "balancing_authority_code_eia", "demand_mwh",
            "net_generation_mwh" and "interchange_mwh" (total interchange, positive
            for exports).
        net_generation_by_energy_source (pd.DataFrame): one row per BA, energy source
            and hour, with columns "datetime_utc", "balancing_authority_code_eia",
            "generation_energy_source" and "net_generation_mwh".
        interchange (pd.DataFrame): one row per pair of BAs and hour, with columns
            "datetime_utc", "balancing_authority_code_eia",
            "balancing_authority_code_adjacent_eia" and "interchange_mwh" (positive for
            exports from "balancing_authority_code_eia").
    """

    operations: pd.DataFrame
    net_generation_by_energy_source: pd.DataFrame
    interchange: pd.DataFrame


def reconcile_eia930(tables: EIA930Tables, output_folder: str) -> EIA930Tables | None:
    """This is the coordinating function for the physics-based reconciliation of
    EIA-930 data.

    Each step writes its results to `output_folder` as parquet files, so they can be
    inspected later:
        - "basic": after basic cleaning
        - "rolling": after rolling window cleaning
        - "weights": the weight of each series in the optimization (values are weights,
          not MWh)
        - "reconciled_without_energy_sources": reconciled data before 2018-07-01,
          when generation by energy source is not reconciled. Only written if the data
          includes hours before that date.
        - "reconciled": reconciled data from 2018-07-01 onwards
    The value of the optimization objective for each hour is written to
    "eia930_objective.parquet" (and "eia930_objective_without_energy_sources.parquet").

    Args:
        tables (EIA930Tables): EIA-930 data to reconcile.
        output_folder (str): folder to write the results of each step to.

    Returns:
        EIA930Tables | None: reconciled EIA-930 data from 2018-07-01 onwards. None if
            the data has no hours from 2018-07-01 onwards.
    """
    start_time = time.time()

    check_eia930_tables(tables)
    series = pivot_eia930_tables_to_series(tables)

    logger.info("Basic data cleaning")
    series = apply_basic_cleaning(series)
    write_eia930_tables(melt_series_to_eia930_tables(series), output_folder, "basic")

    logger.info("Rolling window data cleaning")
    series, weights = apply_rolling_window_cleaning(series)
    write_eia930_tables(melt_series_to_eia930_tables(series), output_folder, "rolling")
    write_eia930_tables(
        melt_series_to_eia930_tables(weights, value_suffix="_weight"),
        output_folder,
        "weights",
    )

    # before generation by energy source was collected, reconcile without it. Note
    # that the hour at GENERATION_BY_ENERGY_SOURCE_START_DATE is included in both
    # periods.
    before_start = series.loc[:GENERATION_BY_ENERGY_SOURCE_START_DATE, :]
    if len(before_start) > 0:
        logger.info(
            "Optimization-based cleaning without generation by energy source, before "
            f"{GENERATION_BY_ENERGY_SOURCE_START_DATE}"
        )
        reconciled, objective = reconcile_energy_balance(
            before_start,
            weights.loc[:GENERATION_BY_ENERGY_SOURCE_START_DATE, :],
            include_energy_sources=False,
        )
        write_eia930_tables(
            melt_series_to_eia930_tables(reconciled),
            output_folder,
            "reconciled_without_energy_sources",
        )
        objective.to_frame().to_parquet(
            os.path.join(
                output_folder, "eia930_objective_without_energy_sources.parquet"
            )
        )

    from_start = series.loc[GENERATION_BY_ENERGY_SOURCE_START_DATE:, :]
    if len(from_start) == 0:
        return None
    logger.info(
        "Optimization-based cleaning with generation by energy source, from "
        f"{GENERATION_BY_ENERGY_SOURCE_START_DATE}"
    )
    reconciled, objective = reconcile_energy_balance(
        from_start,
        weights.loc[GENERATION_BY_ENERGY_SOURCE_START_DATE:, :],
        include_energy_sources=True,
    )
    reconciled_tables = melt_series_to_eia930_tables(reconciled)
    write_eia930_tables(reconciled_tables, output_folder, "reconciled")
    objective.to_frame().to_parquet(
        os.path.join(output_folder, "eia930_objective.parquet")
    )

    logger.info(f"EIA-930 reconciliation took {time.time() - start_time:.0f} seconds")
    return reconciled_tables


def check_eia930_tables(tables: EIA930Tables):
    """Logs which BAs will be treated as outside of the US, and warns about BAs with
    incomplete data.

    BAs that only appear as the adjacent BA of interchange data are treated as being
    outside of the US: their demand, net generation and total interchange are estimated
    from their neighbors' interchange (see `apply_basic_cleaning`). A BA that reports
    generation by energy source or interchange but has no operations data is also
    treated this way, which discards the data it reports, so it gets a warning.

    Args:
        tables (EIA930Tables): EIA-930 data.
    """
    operations_bas = set(tables.operations["balancing_authority_code_eia"])
    generation_bas = set(
        tables.net_generation_by_energy_source["balancing_authority_code_eia"]
    )
    interchange_bas = set(tables.interchange["balancing_authority_code_eia"])
    adjacent_bas = set(tables.interchange["balancing_authority_code_adjacent_eia"])

    incomplete_bas = sorted((generation_bas | interchange_bas) - operations_bas)
    if len(incomplete_bas) > 0:
        logger.warning(
            f"{incomplete_bas} report generation by energy source or interchange but "
            "have no operations data, so they will be treated as outside of the US"
        )
    outside_us_bas = sorted(adjacent_bas - operations_bas)
    logger.info(f"Treating {outside_us_bas} as outside of the US")


def pivot_eia930_tables_to_series(tables: EIA930Tables) -> pd.DataFrame:
    """Pivots EIA-930 tables to a series data frame, with one column per time series.

    Each (BA, hour) row in `operations` becomes a value in three series (demand, net
    generation and total interchange), so every BA in `operations` has all three.

    Args:
        tables (EIA930Tables): EIA-930 data.

    Returns:
        pd.DataFrame: series data frame, indexed by "datetime_utc", with columns
            identified by (field, balancing_authority_code_eia, other) tuples.
    """
    operations = tables.operations.pivot(
        index="datetime_utc",
        columns="balancing_authority_code_eia",
        values=["net_generation_mwh", "interchange_mwh", "demand_mwh"],
    )
    operations_fields = {
        "net_generation_mwh": NET_GENERATION,
        "interchange_mwh": TOTAL_INTERCHANGE,
        "demand_mwh": DEMAND,
    }
    operations.columns = pd.MultiIndex.from_tuples(
        [(operations_fields[value], ba, "") for value, ba in operations.columns],
        names=SERIES_LEVELS,
    )

    net_generation_by_energy_source = tables.net_generation_by_energy_source.pivot(
        index="datetime_utc",
        columns=["balancing_authority_code_eia", "generation_energy_source"],
        values="net_generation_mwh",
    )
    net_generation_by_energy_source.columns = pd.MultiIndex.from_tuples(
        [
            (NET_GENERATION_BY_ENERGY_SOURCE, ba, energy_source)
            for ba, energy_source in net_generation_by_energy_source.columns
        ],
        names=SERIES_LEVELS,
    )

    interchange = tables.interchange.pivot(
        index="datetime_utc",
        columns=[
            "balancing_authority_code_eia",
            "balancing_authority_code_adjacent_eia",
        ],
        values="interchange_mwh",
    )
    interchange.columns = pd.MultiIndex.from_tuples(
        [(INTERCHANGE, ba, adjacent_ba) for ba, adjacent_ba in interchange.columns],
        names=SERIES_LEVELS,
    )

    # series that are missing for some hours (e.g. after a timestamp shift) are blank
    # for those hours. The column names keep the levels when a table is empty, e.g.
    # generation by energy source before 2018-07-01.
    series = pd.concat(
        [net_generation_by_energy_source, operations, interchange], axis="columns"
    ).sort_index()
    series.index.name = "datetime_utc"
    series.columns.names = SERIES_LEVELS
    return series


def melt_series_to_eia930_tables(
    series: pd.DataFrame, value_suffix: str = "_mwh"
) -> EIA930Tables:
    """Melts a series data frame back into EIA-930 tables.

    Every series keeps a row for every hour in `series`, even where its value is blank,
    so that pivoting the tables again gives back the same series.

    Args:
        series (pd.DataFrame): series data frame, as returned by
            `pivot_eia930_tables_to_series`.
        value_suffix (str, optional): suffix of the value columns. Defaults to "_mwh".

    Returns:
        EIA930Tables: EIA-930 data.
    """
    series = series.rename_axis("datetime_utc")
    fields = series.columns.get_level_values("field")

    operations_fields = [DEMAND, NET_GENERATION, TOTAL_INTERCHANGE]
    operations_columns = {
        DEMAND: f"demand{value_suffix}",
        NET_GENERATION: f"net_generation{value_suffix}",
        TOTAL_INTERCHANGE: f"interchange{value_suffix}",
    }
    operations = (
        series.loc[:, fields.isin(operations_fields)]
        .droplevel("other", axis="columns")
        .stack("balancing_authority_code_eia", future_stack=True)
        .reindex(columns=operations_fields)
        .rename(columns=operations_columns)
        .reset_index()
    )
    operations.columns.name = None

    net_generation_by_energy_source = melt_series_with_other_key(
        series.loc[:, fields == NET_GENERATION_BY_ENERGY_SOURCE],
        other_name="generation_energy_source",
        value_name=f"net_generation{value_suffix}",
    )
    interchange = melt_series_with_other_key(
        series.loc[:, fields == INTERCHANGE],
        other_name="balancing_authority_code_adjacent_eia",
        value_name=f"interchange{value_suffix}",
    )

    return EIA930Tables(operations, net_generation_by_energy_source, interchange)


def melt_series_with_other_key(
    series: pd.DataFrame, other_name: str, value_name: str
) -> pd.DataFrame:
    """Melts series of a single field identified by a BA and an "other" key.

    Args:
        series (pd.DataFrame): series data frame with columns of a single field, e.g.
            interchange or generation by energy source. May have no columns.
        other_name (str): name of the column for the "other" key, e.g.
            "generation_energy_source".
        value_name (str): name of the value column.

    Returns:
        pd.DataFrame: long table with columns "datetime_utc",
            "balancing_authority_code_eia", `other_name` and `value_name`.
    """
    columns = ["datetime_utc", "balancing_authority_code_eia", other_name, value_name]
    # e.g. reconciled data before 2018-07-01 has no generation by energy source
    if series.shape[1] == 0:
        return pd.DataFrame(columns=columns)
    return (
        series.droplevel("field", axis="columns")
        .stack(["balancing_authority_code_eia", "other"], future_stack=True)
        .rename(value_name)
        .reset_index()
        .rename(columns={"other": other_name})
    )[columns]


def write_eia930_tables(tables: EIA930Tables, output_folder: str, stage: str):
    """Writes EIA-930 tables to parquet files named "eia930_{stage}_{table}.parquet".

    Args:
        tables (EIA930Tables): EIA-930 data.
        output_folder (str): folder to write the files to.
        stage (str): processing stage, used in the file names.
    """
    os.makedirs(output_folder, exist_ok=True)
    for table_name, table in tables._asdict().items():
        table.to_parquet(
            os.path.join(
                output_folder,
                f"eia930_{stage}_{TABLE_FILE_NAMES[table_name]}.parquet",
            ),
            index=False,
        )


def read_eia930_tables(folder: str, stage: str) -> EIA930Tables:
    """Reads EIA-930 tables written by `write_eia930_tables`.

    Args:
        folder (str): folder the files were written to.
        stage (str): processing stage, used in the file names.

    Returns:
        EIA930Tables: EIA-930 data.
    """
    return EIA930Tables(
        **{
            table_name: pd.read_parquet(
                os.path.join(folder, f"eia930_{stage}_{file_name}.parquet")
            )
            for table_name, file_name in TABLE_FILE_NAMES.items()
        }
    )


def eia930_tables_exist(folder: str, stage: str) -> bool:
    """Returns whether all EIA-930 tables of a stage were written to `folder`.

    Args:
        folder (str): folder the files were written to.
        stage (str): processing stage, used in the file names.

    Returns:
        bool: True if the files for all three tables exist.
    """
    return all(
        os.path.exists(os.path.join(folder, f"eia930_{stage}_{file_name}.parquet"))
        for file_name in TABLE_FILE_NAMES.values()
    )


def get_series_bas(series: pd.DataFrame, field: str | None = None) -> list[str]:
    """Returns the BAs with at least one series, or with a series of a given field.

    Args:
        series (pd.DataFrame): series data frame.
        field (str | None, optional): if given, only return BAs with a series of this
            field. Defaults to None.

    Returns:
        list[str]: BAs, in alphabetical order.
    """
    keys = series.columns
    if field is not None:
        keys = keys[keys.get_level_values("field") == field]
    return sorted(set(keys.get_level_values("balancing_authority_code_eia")))


def apply_basic_cleaning(series: pd.DataFrame) -> pd.DataFrame:
    """Drops unrealistic values and adds series for BAs outside of the US.

    The steps are:
        - Add interchange between AVRN and BPAT if AVRN has no interchange data. AVRN
          only exports to BPAT, so its exports are assumed to be its net generation.
        - Add demand, net generation, total interchange and generation by energy source
          for BAs outside of the US, which only appear as the adjacent BA of
          interchange series. Their total interchange is the opposite of the
          interchange reported by US BAs, exports count as net generation and imports
          as demand. Their generation is assumed to be hydro for Canadian BAs in
          `HYDRO_FOREIGN_BAS` and "other" for the rest.
        - Make sure that every interchange series has a series in the opposite
          direction, adding the opposite of the existing series if needed.
        - Drop values outside of realistic bounds (see `create_reject_bounds`), only
          for BAs in the original data.
        - Set negative solar generation to zero. There are many values of around -50
          MWh at night, from station use.
        - For BAs with no generation by energy source, assign all net generation to
          "other".

    Args:
        series (pd.DataFrame): series data frame.

    Returns:
        pd.DataFrame: cleaned series data frame.
    """
    start_time = time.time()
    series = series.copy()
    # the bounds and the energy source checks only apply to the BAs in the original
    # data, not to the BAs outside of the US added below
    regions = get_series_bas(series)
    interchange_bas = get_series_bas(series, INTERCHANGE)
    adjacent_bas = sorted(
        set(
            series.columns[
                series.columns.get_level_values("field") == INTERCHANGE
            ].get_level_values("other")
        )
    )
    net_generation_bas = get_series_bas(series, NET_GENERATION)

    # AVRN only exports to BPAT. Its interchange is sometimes missing.
    if "AVRN" not in interchange_bas:
        logger.info("Adding interchange for AVRN")
        series[(INTERCHANGE, "AVRN", "BPAT")] = (
            series[(NET_GENERATION, "AVRN", "")] - 1.0
        )
        series[(INTERCHANGE, "BPAT", "AVRN")] = (
            -series[(NET_GENERATION, "AVRN", "")] + 1.0
        )

    # add series for BAs outside of the US, which only appear as adjacent BAs
    foreign_bas = [ba for ba in adjacent_bas if ba not in net_generation_bas]
    logger.info(
        "Adding demand, generation and total interchange for "
        f"{len(foreign_bas)} BAs outside of the US"
    )
    for ba in foreign_bas:
        trade_keys = [
            key for key in series.columns if key[0] == INTERCHANGE and key[2] == ba
        ]
        total_interchange = -series.loc[:, trade_keys].sum(axis=1)
        series[(TOTAL_INTERCHANGE, ba, "")] = total_interchange
        # this keeps blank values blank, as in the original code
        exports = total_interchange.where(~(total_interchange < 0), 0)
        imports = total_interchange.where(~(total_interchange > 0), 0)
        series[(DEMAND, ba, "")] = -imports
        series[(NET_GENERATION, ba, "")] = exports
        energy_source = "hydro" if ba in HYDRO_FOREIGN_BAS else "other"
        series[(NET_GENERATION_BY_ENERGY_SOURCE, ba, energy_source)] = exports
        for key in trade_keys:
            series[(INTERCHANGE, ba, key[1])] = -series.loc[:, key]

    # make sure that interchange series exist in both directions
    for ba in regions:
        for adjacent_ba in regions:
            key = (INTERCHANGE, ba, adjacent_ba)
            opposite_key = (INTERCHANGE, adjacent_ba, ba)
            if key in series.columns and opposite_key not in series.columns:
                logger.info(f"Adding interchange from {adjacent_ba} to {ba}")
                series[opposite_key] = -series.loc[:, key]

    # drop unrealistic values
    reject_bounds = create_reject_bounds(series, regions)
    keys_to_filter = (
        [(DEMAND, ba, "") for ba in regions]
        + [(NET_GENERATION, ba, "") for ba in regions]
        + [(TOTAL_INTERCHANGE, ba, "") for ba in regions]
        + get_interchange_keys(series, regions)
    )
    for key in keys_to_filter:
        series[key] = filter_to_bounds(series.loc[:, key], reject_bounds[key])

    # do the same for generation by energy source, keeping values of solar down to
    # -1 GW and then setting negative values to zero
    for ba in regions:
        has_energy_sources = False
        for energy_source in ENERGY_SOURCES:
            key = (NET_GENERATION_BY_ENERGY_SOURCE, ba, energy_source)
            if key in series.columns:
                has_energy_sources = True
                if energy_source == "solar":
                    reject_bounds[key] = (-1e3, 200e3)
                series[key] = filter_to_bounds(series.loc[:, key], reject_bounds[key])
                if energy_source == "solar":
                    series[key] = series.loc[:, key].where(~(series.loc[:, key] < 0), 0)
        if not has_energy_sources:
            series[(NET_GENERATION_BY_ENERGY_SOURCE, ba, "other")] = series.loc[
                :, (NET_GENERATION, ba, "")
            ]

    series.columns.names = SERIES_LEVELS
    logger.info(f"Basic cleaning took {time.time() - start_time:.2f} seconds")
    return series


def get_interchange_keys(series: pd.DataFrame, bas: list[str]) -> list[tuple]:
    """Returns the keys of interchange series between pairs of BAs in `bas`.

    Args:
        series (pd.DataFrame): series data frame.
        bas (list[str]): BAs to include.

    Returns:
        list[tuple]: keys of the interchange series that exist, ordered by BA and then
            by adjacent BA, in the order of `bas`.
    """
    return [
        (INTERCHANGE, ba, adjacent_ba)
        for ba in bas
        for adjacent_ba in bas
        if (INTERCHANGE, ba, adjacent_ba) in series.columns
    ]


def create_reject_bounds(series: pd.DataFrame, regions: list[str]) -> defaultdict:
    """Creates the ranges outside of which values are considered unrealistic.

    The default range is -1 MWh to 200 GWh. Interchange and total interchange can be
    between -100 GWh and 100 GWh. A few series have their own ranges.

    Args:
        series (pd.DataFrame): series data frame.
        regions (list[str]): BAs in the original data.

    Returns:
        defaultdict: (lower bound, upper bound) in MWh, keyed by series key.
    """
    reject_bounds = defaultdict(lambda: (-1.0, 200e3))
    for ba in regions:
        reject_bounds[(TOTAL_INTERCHANGE, ba, "")] = (-100e3, 100e3)
    for key in get_interchange_keys(series, regions):
        reject_bounds[key] = (-100e3, 100e3)
    reject_bounds[(DEMAND, "AZPS", "")] = (1.0, 30e3)
    reject_bounds[(DEMAND, "BANC", "")] = (1.0, 6.5e3)
    reject_bounds[(TOTAL_INTERCHANGE, "BANC", "")] = (-5e3, 5e3)
    reject_bounds[(NET_GENERATION, "CISO", "")] = (5e3, 60e3)
    return reject_bounds


def filter_to_bounds(values: pd.Series, bounds: tuple[float, float]) -> pd.Series:
    """Blanks out values outside of `bounds`.

    Args:
        values (pd.Series): values to filter.
        bounds (tuple[float, float]): lower and upper bounds, inclusive.

    Returns:
        pd.Series: values, blank where outside of the bounds.
    """
    return values.where((values >= bounds[0]) & (values <= bounds[1]))


def apply_rolling_window_cleaning(
    series: pd.DataFrame, number_of_runs: int = 2
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drops anomalous values using a rolling window filter, then fills gaps.

    Large outliers can mask smaller ones, so the filter is applied `number_of_runs`
    times. Gaps are then filled with the average of the same hour on the previous and
    next days, up to four times, and then without limit. Series with no data at all are
    set to zero. Any remaining gaps are filled with the rolling median.

    Args:
        series (pd.DataFrame): series data frame.
        number_of_runs (int, optional): number of times to apply the rolling window
            filter. Defaults to 2.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: cleaned series data frame, and the weight of
            each value in the optimization, based on the rolling mean of each series.
    """
    start_time = time.time()
    series = series.copy()

    for _ in range(number_of_runs):
        series = filter_with_rolling_window(series, replace_nan_with_median=False)

    # first fill gaps with the average of the same hour on the previous and next days.
    # There are strong daily patterns, so this works well. Fill one day at a time,
    # which makes for a smoother transition across multi-day gaps. After four passes,
    # fill without a limit.
    for key in series.columns:
        number_of_passes = 0
        while (series.loc[:, key].isna().sum() > 0) and (number_of_passes < 4):
            number_of_passes += 1
            series[key] = pd.concat(
                [
                    series.loc[:, key].groupby(series.index.hour).ffill(limit=1),
                    series.loc[:, key].groupby(series.index.hour).bfill(limit=1),
                ],
                axis=1,
            ).mean(axis=1)
        if number_of_passes == 4:
            logger.debug(f"A lot of bad data for {key}")
            series[key] = pd.concat(
                [
                    series.loc[:, key].groupby(series.index.hour).ffill(),
                    series.loc[:, key].groupby(series.index.hour).bfill(),
                ],
                axis=1,
            ).mean(axis=1)

        # series with no data at all
        if series.loc[:, key].isna().sum() == len(series):
            series[key] = 0.0

    # fill any remaining gaps with the rolling median
    series, rolling_mean = filter_with_rolling_window(
        series, replace_nan_with_median=True, return_mean=True
    )

    if series.isna().sum().sum() > 0:
        logger.warning("There are still blank values after rolling window cleaning")

    # weights are inversely proportional to the size of each series, so that the
    # optimization adjusts large series more than small ones
    absolute_mean = rolling_mean.abs()
    weights = WEIGHT_NUMERATOR / absolute_mean.where(
        absolute_mean > WEIGHT_MINIMUM_DENOMINATOR, WEIGHT_MINIMUM_DENOMINATOR
    )

    logger.info(f"Rolling window cleaning took {time.time() - start_time:.2f} seconds")
    return series, weights


def filter_with_rolling_window(
    series: pd.DataFrame,
    window: int = 20 * 24,
    min_periods: int = 100,
    replace_nan_with_median: bool = True,
    return_mean: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Blanks out values far from the rolling median of each series.

    Values are rejected if they are more than 3 times the rolling 5th to 95th
    percentile range away from the rolling median, using a centered window. The range
    is at least 150 MWh, so that small series are not over-filtered.

    Note that the values of `series` are modified in place.

    Args:
        series (pd.DataFrame): series data frame.
        window (int, optional): size of the rolling window, in hours. Defaults to 480.
        min_periods (int, optional): minimum number of values in the window. Defaults
            to 100.
        replace_nan_with_median (bool, optional): whether to fill gaps after filtering,
            by linear interpolation for gaps of up to 3 hours and with the rolling
            median otherwise. Defaults to True.
        return_mean (bool, optional): whether to also return the rolling mean of each
            series. Defaults to False.

    Returns:
        pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]: filtered series data frame,
            and its rolling mean if `return_mean` is True.
    """
    for key in series.columns:
        rolling = series[key].rolling(window, min_periods=min_periods, center=True)
        median = rolling.median()
        # at least 150 MWh, including where the range can't be computed
        percentile_range = rolling.quantile(0.95) - rolling.quantile(0.05)
        percentile_range = percentile_range.where(percentile_range > 150, 150)
        upper_bound = median + 3 * percentile_range
        lower_bound = median - 3 * percentile_range
        to_reject = (series[key] >= upper_bound) | (series[key] <= lower_bound)
        series.loc[to_reject, key] = np.nan
        if replace_nan_with_median:
            series.loc[:, key] = series.loc[:, key].interpolate(limit=3)
            is_blank = series[key].isnull()
            series.loc[is_blank, key] = median.loc[is_blank]

    if return_mean:
        rolling_mean = series.rolling(
            window, min_periods=min_periods, center=True
        ).mean()
        return series, rolling_mean
    return series


def reconcile_energy_balance(
    series: pd.DataFrame, weights: pd.DataFrame, include_energy_sources: bool
) -> tuple[pd.DataFrame, pd.Series]:
    """Makes the data physically consistent in each hour.

    For each hour, finds the smallest weighted adjustments to the data such that:
        - demand + total interchange = net generation, for each BA
        - total interchange = sum of interchange with each adjacent BA, for each BA
        - interchange from A to B = -(interchange from B to A)
        - net generation = sum of generation by energy source, for each BA (only if
          `include_energy_sources`)
        - demand, net generation and generation by energy source are at least 1 MWh

    The problem is the same for every hour except for its data and weights, so it is
    built once (see `HourlyReconciliationProblem`) and solved for each hour in turn.
    Blank values and weights are treated as zero.

    Args:
        series (pd.DataFrame): series data frame, with no blank values except where a
            series has no data.
        weights (pd.DataFrame): weight of each value, with the same shape as `series`.
        include_energy_sources (bool): whether to reconcile generation by energy
            source.

    Returns:
        tuple[pd.DataFrame, pd.Series]: reconciled series data frame, and the value of
            the optimization objective for each hour.
    """
    start_time = time.time()
    logger.info(f"Running the reconciliation for {len(series)} hours")
    series = series.fillna(0)
    weights = weights.fillna(0)

    problem = HourlyReconciliationProblem(series, include_energy_sources)
    values = series.loc[:, problem.keys].to_numpy()
    value_weights = weights.loc[:, problem.keys].to_numpy()
    results = [
        problem.solve(values[hour], value_weights[hour]) for hour in range(len(series))
    ]

    reconciled = pd.DataFrame(
        np.array([reconciled_values for reconciled_values, _ in results]),
        index=series.index,
        columns=pd.MultiIndex.from_tuples(problem.keys, names=SERIES_LEVELS),
    )
    objective = pd.Series(
        [objective for _, objective in results], index=series.index, name="objective"
    )

    check_energy_balance(reconciled)
    logger.info(f"Reconciliation took {time.time() - start_time:.2f} seconds")
    return reconciled, objective


class HourlyReconciliationProblem:
    """The quadratic program that reconciles the data of one hour.

    The problem is built from the series in the data, with their values and weights as
    parameters, so cvxpy compiles it only once. `solve` then sets the parameters for an
    hour and solves it with Clarabel, which is much faster than building the problem
    again for each hour.

    The objective is the sum of the squared adjustments, each multiplied by its weight.
    See `reconcile_energy_balance` for the constraints.

    Attributes:
        keys (list[tuple]): keys of the series in the problem, in the order of the
            values passed to and returned by `solve`. Generation by energy source is
            only included if `include_energy_sources`, and only for the energy sources
            in `ENERGY_SOURCES`.
    """

    def __init__(self, series: pd.DataFrame, include_energy_sources: bool):
        """Builds the problem for the series of a series data frame.

        Args:
            series (pd.DataFrame): series data frame. Only its columns are used. Every
                BA must have demand, net generation and total interchange series, and
                every interchange series must have a series in the opposite direction.
            include_energy_sources (bool): whether to reconcile generation by energy
                source.
        """
        regions = get_series_bas(series)
        # (BA, adjacent BA) of each interchange series, and (energy source, BA) of each
        # generation by energy source series
        interchange_pairs = [
            (ba, adjacent_ba)
            for ba in regions
            for adjacent_ba in regions
            if (INTERCHANGE, ba, adjacent_ba) in series.columns
        ]
        energy_source_pairs = [
            (energy_source, ba)
            for energy_source in ENERGY_SOURCES
            for ba in regions
            if (NET_GENERATION_BY_ENERGY_SOURCE, ba, energy_source) in series.columns
        ]

        keys_by_field = {
            DEMAND: [(DEMAND, ba, "") for ba in regions],
            NET_GENERATION: [(NET_GENERATION, ba, "") for ba in regions],
            TOTAL_INTERCHANGE: [(TOTAL_INTERCHANGE, ba, "") for ba in regions],
            INTERCHANGE: [(INTERCHANGE, ba, adj) for ba, adj in interchange_pairs],
        }
        if include_energy_sources:
            keys_by_field[NET_GENERATION_BY_ENERGY_SOURCE] = [
                (NET_GENERATION_BY_ENERGY_SOURCE, ba, energy_source)
                for energy_source, ba in energy_source_pairs
            ]
        self.keys = [key for keys in keys_by_field.values() for key in keys]

        # the data and the square root of the weights are set for each hour by `solve`
        self.values = cp.Parameter(len(self.keys))
        self.sqrt_weights = cp.Parameter(len(self.keys), nonneg=True)
        self.adjustments = cp.Variable(len(self.keys))

        # reconciled values of each field, as slices of the reconciled data
        reconciled = self.values + self.adjustments
        reconciled_by_field = {}
        start = 0
        for field, keys in keys_by_field.items():
            reconciled_by_field[field] = reconciled[start : start + len(keys)]
            start += len(keys)
        demand = reconciled_by_field[DEMAND]
        net_generation = reconciled_by_field[NET_GENERATION]
        total_interchange = reconciled_by_field[TOTAL_INTERCHANGE]
        interchange = reconciled_by_field[INTERCHANGE]

        # sums the interchange series of each BA
        interchange_to_ba = create_incidence_matrix(
            regions, [ba for ba, _ in interchange_pairs]
        )
        # the position of the series in the opposite direction of each interchange
        # series
        pair_positions = {pair: i for i, pair in enumerate(interchange_pairs)}
        opposite = [pair_positions[(adj, ba)] for ba, adj in interchange_pairs]

        constraints = [
            demand >= OPTIMIZATION_LOWER_BOUND,
            net_generation >= OPTIMIZATION_LOWER_BOUND,
            # energy balance of each BA
            demand + total_interchange == net_generation,
            # total interchange is the sum of interchange with each adjacent BA
            total_interchange == interchange_to_ba @ interchange,
            # interchange from A to B is the opposite of interchange from B to A. Each
            # pair of BAs is constrained twice, which is redundant but harmless.
            interchange == -interchange[opposite],
        ]
        if include_energy_sources:
            generation = reconciled_by_field[NET_GENERATION_BY_ENERGY_SOURCE]
            # sums the generation by energy source series of each BA
            energy_sources_to_ba = create_incidence_matrix(
                regions, [ba for _, ba in energy_source_pairs]
            )
            constraints += [
                generation >= OPTIMIZATION_LOWER_BOUND,
                # net generation is the sum of generation by energy source
                net_generation == energy_sources_to_ba @ generation,
            ]

        # sum of weight * adjustment ** 2
        objective = cp.sum_squares(cp.multiply(self.sqrt_weights, self.adjustments))
        self.problem = cp.Problem(cp.Minimize(objective), constraints)

    def solve(
        self, values: np.ndarray, weights: np.ndarray
    ) -> tuple[np.ndarray, float]:
        """Reconciles the data of one hour.

        Args:
            values (np.ndarray): value of each series in `keys`, with no blank values.
            weights (np.ndarray): weight of each series in `keys`, with no blank
                values.

        Raises:
            ValueError: if there are blank values, or if the solver does not find a
                solution.

        Returns:
            tuple[np.ndarray, float]: reconciled value of each series in `keys`, and
                the value of the objective.
        """
        if np.isnan(values).any() or np.isnan(weights).any():
            raise ValueError("Cannot reconcile data with blank values")
        self.values.value = values
        self.sqrt_weights.value = np.sqrt(weights)
        self.problem.solve(solver=cp.CLARABEL, **CLARABEL_SETTINGS)
        if self.problem.status != cp.OPTIMAL:
            raise ValueError(
                f"The reconciliation solver stopped with status {self.problem.status}"
            )
        return values + self.adjustments.value, self.problem.value


def create_incidence_matrix(
    row_bas: list[str], column_bas: list[str]
) -> scipy.sparse.csr_array:
    """Creates a matrix that sums the series of each BA.

    Args:
        row_bas (list[str]): BA of each row.
        column_bas (list[str]): BA of each column, i.e. of each series to sum.

    Returns:
        scipy.sparse.csr_array: matrix with one row per BA in `row_bas` and one column
            per series, with 1 where the series belongs to the BA and 0 elsewhere.
    """
    row_positions = {ba: i for i, ba in enumerate(row_bas)}
    rows = [row_positions[ba] for ba in column_bas]
    columns = range(len(column_bas))
    return scipy.sparse.csr_array(
        (np.ones(len(column_bas)), (rows, columns)),
        shape=(len(row_bas), len(column_bas)),
    )


def check_energy_balance(series: pd.DataFrame, tolerance: float = 1e-2):
    """Logs an error for each BA whose reconciled data are not physically consistent.

    Checks for each BA that there are no blank values, that demand + total interchange
    = net generation, that total interchange = sum of interchange, that interchange is
    antisymmetric, and that demand and net generation are not negative.

    Args:
        series (pd.DataFrame): reconciled series data frame.
        tolerance (float, optional): largest difference allowed, in MWh. Defaults to
            1e-2.
    """
    for ba in get_series_bas(series):
        partners = [
            adjacent_ba
            for adjacent_ba in get_series_bas(series)
            if (INTERCHANGE, ba, adjacent_ba) in series.columns
            and (INTERCHANGE, adjacent_ba, ba) in series.columns
        ]
        demand = series[(DEMAND, ba, "")]
        net_generation = series[(NET_GENERATION, ba, "")]
        total_interchange = series[(TOTAL_INTERCHANGE, ba, "")]

        for field, values in [
            (DEMAND, demand),
            (NET_GENERATION, net_generation),
            (TOTAL_INTERCHANGE, total_interchange),
        ]:
            if values.isna().sum() != 0:
                logger.error(
                    f"There are still {values.isna().sum()} blanks for {ba} {field}"
                )
        for adjacent_ba in partners:
            blanks = series[(INTERCHANGE, ba, adjacent_ba)].isna().sum()
            if blanks != 0:
                logger.error(f"There are still {blanks} blanks for {ba}-{adjacent_ba}")

        if ((net_generation - (demand + total_interchange)).abs() > tolerance).sum():
            logger.error(f"{ba}: demand + total interchange == net generation violated")

        interchange_sum = series.loc[
            :, [(INTERCHANGE, ba, adjacent_ba) for adjacent_ba in partners]
        ].sum(axis=1)
        if ((total_interchange - interchange_sum).abs() > tolerance).sum():
            logger.error(f"{ba}: total interchange == sum of interchange violated")

        for adjacent_ba in partners:
            antisymmetry = (
                series[(INTERCHANGE, ba, adjacent_ba)]
                + series[(INTERCHANGE, adjacent_ba, ba)]
            )
            if (antisymmetry.abs() > tolerance).sum():
                logger.error(f"{ba}-{adjacent_ba}: interchange antisymmetry violated")

        for field, values in [(DEMAND, demand), (NET_GENERATION, net_generation)]:
            if (values < 0).sum():
                logger.error(
                    f"{ba}: there are {(values < 0).sum()} negative values for {field}"
                )
