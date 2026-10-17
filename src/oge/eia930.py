import pandas as pd
from datetime import timedelta

import oge.load_data as load_data
from oge.filepaths import (
    outputs_folder,
    reference_table_folder,
)
from oge.logging_util import get_logger
from oge.physics_reconciliation import (
    DEMAND,
    ENERGY_SOURCES,
    INTERCHANGE,
    NET_GENERATION,
    NET_GENERATION_BY_ENERGY_SOURCE,
    NEW_ENERGY_SOURCES,
    TOTAL_INTERCHANGE,
    EIA930Tables,
    melt_series_to_eia930_tables,
    pivot_eia930_tables_to_series,
    reconcile_eia930,
    write_eia930_tables,
)

logger = get_logger(__name__)


# PUDL renamed hydro, solar and wind on 2024-07-01, for every BA, when EIA-930
# introduced pumped storage and integrated battery storage as separate energy sources.
# The new names don't mean that a BA reports the new energy sources (BAs started to
# report those at different times), so each pair is combined into one series, under
# the original name.
renamed_energy_sources = {
    "hydro_excluding_pumped_storage": "hydro",
    "solar_wo_integrated_battery_storage": "solar",
    "wind_wo_integrated_battery_storage": "wind",
}

# Dictionary mapping EIA-930 energy sources to the fuel categories in
# reference_tables/energy_source_groups.csv. The energy sources that EIA-930
# introduced in 2024 are mapped to the category that included them before.
# "unknown" is not mapped.
energy_source_to_fuel_category_eia930 = {
    "coal": "coal",
    "gas": "natural_gas",
    "other": "other",
    "hydro": "hydro",
    "wind": "wind",
    "solar": "solar",
    "nuclear": "nuclear",
    "oil": "petroleum",
    "pumped_storage": "hydro",
    "solar_w_integrated_battery_storage": "solar",
    "wind_w_integrated_battery_storage": "wind",
    "battery_storage": "other",
    "unknown_energy_storage": "other",
    "other_energy_storage": "other",
    "geothermal": "other",
}


def clean_930(year: int, path_prefix: str = "") -> EIA930Tables | None:
    """This is the coordinating function for loading, adjusting and reconciling
    EIA-930 data.

    Writes the data at each step to `outputs/{path_prefix}/eia930/` as parquet files:
    "unadjusted" (as loaded from PUDL), "raw" (after timestamp adjustments) and the
    steps written by `physics_reconciliation.reconcile_eia930`.

    Args:
        year (int): a four-digit year.
        path_prefix (str): path to folder enclosing the eia930 folder with files needed
            for the physics-based data cleaning

    Returns:
        EIA930Tables | None: reconciled EIA-930 data. None if the data has no hours
            from 2018-07-01 onwards.
    """
    data_folder = outputs_folder(f"{path_prefix}/eia930/")

    tables = load_eia930_balance_data(year)
    write_eia930_tables(tables, data_folder, "unadjusted")

    logger.info("Adjusting EIA-930 time stamps")
    tables = manual_930_adjust(tables)
    write_eia930_tables(tables, data_folder, "raw")

    logger.info("Running physics-based data cleaning")
    # local timezone of each BA, so that the months in which generation by energy
    # source is zero match the months used for hourly shaping
    ba_reference = load_data.load_ba_reference().dropna(subset=["timezone_local"])
    ba_timezones = dict(zip(ba_reference["ba_code"], ba_reference["timezone_local"]))
    return reconcile_eia930(tables, data_folder, ba_timezones)


def load_eia930_balance_data(year: int) -> EIA930Tables:
    """Loads EIA-930 operations, generation by energy source and interchange from PUDL.

    Loads data from October of the previous year through January 2 of the next year,
    since the rolling window cleaning needs data before and after the year. The energy
    sources that EIA-930 introduced in 2024 (e.g. battery storage) are kept as their
    own series, and the energy sources that PUDL renamed on 2024-07-01 are combined
    with their original names (see `combine_renamed_energy_sources`).

    Args:
        year (int): a four-digit year.

    Returns:
        EIA930Tables: EIA-930 data, with adjusted net generation and interchange and
            demand imputed by PUDL.
    """
    start = pd.Timestamp(year - 1, 10, 1)
    end = pd.Timestamp(year + 1, 1, 2)

    interchange = load_data.load_pudl_table(
        "core_eia930__hourly_interchange",
        columns=[
            "datetime_utc",
            "balancing_authority_code_eia",
            "balancing_authority_code_adjacent_eia",
            "interchange_reported_mwh",
        ],
        dt=start,
        end_dt=end,
    ).rename(columns={"interchange_reported_mwh": "interchange_mwh"})

    net_generation_by_energy_source = combine_renamed_energy_sources(
        load_data.load_pudl_table(
            "core_eia930__hourly_net_generation_by_energy_source",
            columns=[
                "datetime_utc",
                "balancing_authority_code_eia",
                "generation_energy_source",
                "net_generation_adjusted_mwh",
            ],
            dt=start,
            end_dt=end,
        ).rename(columns={"net_generation_adjusted_mwh": "net_generation_mwh"})
    )

    operations = load_data.load_pudl_table(
        "out_eia930__hourly_operations",
        columns=[
            "datetime_utc",
            "balancing_authority_code_eia",
            "net_generation_adjusted_mwh",
            "interchange_adjusted_mwh",
            "demand_imputed_pudl_mwh",
        ],
        dt=start,
        end_dt=end,
    ).rename(
        columns={
            "net_generation_adjusted_mwh": "net_generation_mwh",
            "interchange_adjusted_mwh": "interchange_mwh",
            "demand_imputed_pudl_mwh": "demand_mwh",
        }
    )[
        [
            "datetime_utc",
            "balancing_authority_code_eia",
            "demand_mwh",
            "net_generation_mwh",
            "interchange_mwh",
        ]
    ]

    tables = add_missing_operations(
        EIA930Tables(operations, net_generation_by_energy_source, interchange)
    )
    check_bas_outside_us(tables)
    return tables


def combine_renamed_energy_sources(
    net_generation_by_energy_source: pd.DataFrame,
) -> pd.DataFrame:
    """Combines the energy sources that PUDL renamed on 2024-07-01 with their original
    names, and drops new energy sources that a BA never reports.

    PUDL's table has a row for every BA, hour and energy source, blank when nothing is
    reported. Each renamed energy source (see `renamed_energy_sources`) is combined
    with its original name, which is blank on the other side of the rename. The energy
    sources that EIA-930 introduced in 2024 (`NEW_ENERGY_SOURCES`) are kept as their
    own series, except for the BAs that don't report them in any hour, which would only
    add series of zeros to the reconciliation. The original energy sources are kept for
    every BA, as before the new ones were introduced.

    Args:
        net_generation_by_energy_source (pd.DataFrame): PUDL's generation by energy
            source, with columns "datetime_utc", "balancing_authority_code_eia",
            "generation_energy_source" and "net_generation_mwh".

    Returns:
        pd.DataFrame: generation by energy source with one row per BA, hour and energy
            source.
    """
    keys = ["datetime_utc", "balancing_authority_code_eia", "generation_energy_source"]
    data = net_generation_by_energy_source.astype(
        {"generation_energy_source": "str"}
    ).replace({"generation_energy_source": renamed_energy_sources})

    # a BA should report each pair under only one of the names in any hour
    reported_twice = (
        data.groupby(keys, dropna=True)["net_generation_mwh"].count() > 1
    )  # the keys are never blank in PUDL
    if reported_twice.any():
        logger.warning(
            "Some BAs report an energy source under both its original and its "
            "renamed PUDL name in the same hour. The values are added together:\n"
            + reported_twice[reported_twice]
            .reset_index()
            .groupby(keys[1:], dropna=False)
            .size()
            .rename("hours")
            .to_string()
        )

    # new energy sources that a BA never reports
    is_reported = (
        data.groupby(keys[1:], dropna=True)["net_generation_mwh"]
        .transform("count")
        .gt(0)
    )
    is_new = data["generation_energy_source"].isin(NEW_ENERGY_SOURCES)
    data = data[is_reported | ~is_new]

    # note that sum() turns groups with only blank values into zero, as in the
    # original gridemissions-format code
    return data.groupby(keys, dropna=True)["net_generation_mwh"].sum().reset_index()


def add_missing_operations(tables: EIA930Tables) -> EIA930Tables:
    """Adds operations data for BAs that report generation or interchange but have no
    operations data.

    PUDL's operations table can miss BAs that EIA-930 does report, e.g. SIKE, a
    generation-only BA that started reporting on 2025-06-01 with net generation and
    total interchange but no demand. Without operations data, the reconciliation would
    treat such a BA as outside of the US and estimate all of its data from its
    neighbors' interchange.

    For these BAs, net generation is the sum of their generation by energy source,
    total interchange is the sum of their reported interchange with each adjacent BA,
    and demand is net generation minus total interchange. Hours where a value can't be
    computed are left blank.

    Args:
        tables (EIA930Tables): EIA-930 data.

    Returns:
        EIA930Tables: EIA-930 data, with operations data added for BAs that were
            missing it.
    """
    operations_bas = set(tables.operations["balancing_authority_code_eia"])
    reporting_bas = set(
        tables.net_generation_by_energy_source["balancing_authority_code_eia"]
    ) | set(tables.interchange["balancing_authority_code_eia"])
    missing_bas = sorted(reporting_bas - operations_bas)
    if len(missing_bas) == 0:
        return tables
    logger.warning(
        f"Adding EIA-930 operations data for {missing_bas}, which report generation "
        "or interchange but have no operations data in PUDL"
    )

    # min_count=1 keeps hours with no data blank
    keys = ["datetime_utc", "balancing_authority_code_eia"]
    net_generation = (
        tables.net_generation_by_energy_source[
            tables.net_generation_by_energy_source["balancing_authority_code_eia"].isin(
                missing_bas
            )
        ]
        .groupby(keys, dropna=False)["net_generation_mwh"]
        .sum(min_count=1)
    )
    total_interchange = (
        tables.interchange[
            tables.interchange["balancing_authority_code_eia"].isin(missing_bas)
        ]
        .groupby(keys, dropna=False)["interchange_mwh"]
        .sum(min_count=1)
    )
    added_operations = pd.concat(
        [net_generation, total_interchange], axis="columns"
    ).reset_index()
    added_operations["demand_mwh"] = (
        added_operations["net_generation_mwh"] - added_operations["interchange_mwh"]
    )

    operations = pd.concat(
        [tables.operations, added_operations[tables.operations.columns]],
        ignore_index=True,
    )
    return tables._replace(operations=operations)


def check_bas_outside_us(tables: EIA930Tables):
    """Warns about BAs that will be treated as outside of the US but are listed as US
    BAs in `ba_reference.csv`.

    The reconciliation treats BAs that only appear as the adjacent BA of interchange
    data as being outside of the US, and estimates their data from their neighbors'
    interchange (see `physics_reconciliation.apply_basic_cleaning`).

    Args:
        tables (EIA930Tables): EIA-930 data.
    """
    adjacent_only_bas = set(
        tables.interchange["balancing_authority_code_adjacent_eia"]
    ) - set(tables.operations["balancing_authority_code_eia"])
    ba_ref = pd.read_csv(reference_table_folder("ba_reference.csv"))
    us_bas = set(ba_ref.loc[ba_ref["us_ba"] == "Yes", "ba_code"])
    us_adjacent_only_bas = sorted(adjacent_only_bas & us_bas)
    if len(us_adjacent_only_bas) > 0:
        logger.warning(
            f"{us_adjacent_only_bas} only appear as adjacent BAs in EIA-930 interchange "
            "data, so the reconciliation will treat them as outside of the US"
        )


def prepare_eia930_for_shaping(
    net_generation_by_energy_source: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Formats reconciled EIA-930 generation by energy source for hourly shaping.

    Keeps the year plus or minus one day, to allow for the conversion to local time
    later, drops BAs outside of the US, adds the energy sources that EIA-930 introduced
    in 2024 to the fuel category that included them before (e.g. battery storage to
    "other"), and sets negative generation to zero.

    Args:
        net_generation_by_energy_source (pd.DataFrame): reconciled EIA-930 generation
            by energy source.
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: EIA-930 generation with columns "ba_code",
            "fuel_category_eia930", "datetime_utc", "datetime_local", "report_date" and
            "net_generation_mwh_930".
    """
    data = net_generation_by_energy_source.rename(
        columns={
            "balancing_authority_code_eia": "ba_code",
            "net_generation_mwh": "net_generation_mwh_930",
        }
    )
    data["datetime_utc"] = data["datetime_utc"].astype("datetime64[ns, UTC]")

    # only keep data for the single year we are interested, plus or minus one
    # day, to account for conversion to local time later
    data = data.loc[
        (data["datetime_utc"] >= f"{year - 1}-12-31")
        & (data["datetime_utc"] < f"{year + 1}-01-02")
    ]

    # drop BAs not located in the United States
    ba_ref = pd.read_csv(reference_table_folder("ba_reference.csv"))
    foreign_bas = list(ba_ref.loc[ba_ref["us_ba"] == "No", "ba_code"])
    data = data[~data["ba_code"].isin(foreign_bas)]

    # rename the energy sources using the fuel categories in
    # reference_tables/energy_source_groups.csv, and add up the energy sources in each
    # category. The new energy sources are added before negative values are set to
    # zero below, so that e.g. battery charging still offsets other generation in the
    # "other" category, as it did before BAs reported batteries separately. Blank
    # categories ("unknown") are kept.
    data = (
        data.assign(
            fuel_category_eia930=data["generation_energy_source"].map(
                energy_source_to_fuel_category_eia930
            )
        )
        .groupby(["ba_code", "fuel_category_eia930", "datetime_utc"], dropna=False)[
            "net_generation_mwh_930"
        ]
        .sum(min_count=1)
        .reset_index()
    )

    data["datetime_local"] = ""
    for ba in list(data["ba_code"].unique()):
        data.loc[data.ba_code == ba, "datetime_local"] = (
            data.loc[data.ba_code == ba, "datetime_utc"]
            .dt.tz_convert(load_data.ba_timezone(ba=ba, type="local"))
            .astype(str)
        )

    # create a report date column
    data["report_date"] = data["datetime_local"].str[:7]
    data["report_date"] = pd.to_datetime(data["report_date"]).astype("datetime64[s]")

    # the hourly profiles describe when positive generation happens, so hours when an
    # energy source is a net load (e.g. station use at night, or charging storage) are
    # set to zero. Otherwise, a month whose total is close to zero or negative would
    # give profiles with very large or inverted shares of the monthly total.
    data["net_generation_mwh_930"] = data["net_generation_mwh_930"].clip(lower=0)

    data = data[
        [
            "ba_code",
            "fuel_category_eia930",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "net_generation_mwh_930",
        ]
    ].reset_index(drop=True)

    return data


def remove_months_with_zero_data(eia930_data: pd.DataFrame) -> pd.DataFrame:
    """Removes the BA-fuel-months in which EIA-930 generation is zero in every hour.

    An all-zero month has no shape to give, so it is removed from the data used to
    calculate hourly profiles. If EIA-923 reports generation for that BA, fuel and
    month, the month is then treated as missing and gets an imputed profile (see
    `impute_hourly_profiles.impute_missing_hourly_profiles`).

    Generation by energy source that is zero for a whole month before reconciliation
    stays exactly zero after it (see `physics_reconciliation.find_zero_months`). Months
    with only negative values (e.g. storage charging) are also all zero, since
    negative values are set to zero in `prepare_eia930_for_shaping`.

    Args:
        eia930_data (pd.DataFrame): EIA-930 generation, as returned by
            `prepare_eia930_for_shaping`.

    Returns:
        pd.DataFrame: `eia930_data` without the BA-fuel-months that are zero in every
            hour.
    """
    keys = ["ba_code", "fuel_category_eia930", "report_date"]
    is_zero_month = (
        eia930_data.groupby(keys, dropna=False)["net_generation_mwh_930"]
        .transform(lambda values: (values == 0).all())
        .astype(bool)
    )
    logger.info(
        f"Removing {eia930_data.loc[is_zero_month, keys].drop_duplicates().shape[0]} "
        "EIA-930 BA-fuel-months with zero generation in every hour"
    )
    return eia930_data[~is_zero_month].reset_index(drop=True)


###########################################################
# Code for adjusting 930 timestamps
#
###########################################################


def manual_930_adjust(tables: EIA930Tables) -> EIA930Tables:
    """Adjust timestamps in EIA-930 data.

    Adjustment Steps:
        - Make all end-of-hour
            - Generation
                - PJM: + 1 hour.
                - CISO: + 1 hour before June 16, 2022.
                        - 1 hour from November 1, 2023 to December 1, 2024.
                - TEPC: + 1 hour. Fixed after Nov 1, 2021.
                - SC: - 4 hours during daylight savings hours.
                      - 5 hours during standard hours (this happens to be
                      Eastern <-> UTC offset)
                      Fixed after Jan 1, 2021
                - AVA: + 1 hour during daylight saving hours in 2024.
                - BANC: + 1 hour during first month of saving hours in 2024.
            - Interchange
                - PJM: + 4 hours
                - TEPC: + 7 hours
                - IID:  + 4 hours
            - Interchange sign flipped
                - PJM-{CPLE, CPLW, DUK, LGEE, MISO, NYIS, TVA} before 10/31/2019,
                  4:00 UTC. This affects all interchange partners except OVEC, and
                  excluding total interchange
            - Other interchange
                - AZPS - SRP flips gradually in Nov 2019, then abruptly back in
                  06/2020 throughout, SRP - AZPS remains constant around 3000 MWh
                  imported to AZPS from SRP. We assume SRP - AZPS is correct, and
                  assign AZPS - SRP to be the inverse
        - Make all start-of-hour.

    Args:
        tables (EIA930Tables): EIA-930 data to adjust.

    Returns:
        EIA930Tables: adjusted EIA-930 data.
    """
    return melt_series_to_eia930_tables(
        adjust_series_timestamps(pivot_eia930_tables_to_series(tables))
    )


def adjust_series_timestamps(raw: pd.DataFrame) -> pd.DataFrame:
    """Applies the adjustments described in `manual_930_adjust` to a series data frame.

    Shifting a series by whole hours adds hours at the start or end of the data, where
    the other series are blank.

    Args:
        raw (pd.DataFrame): series data frame, as returned by
            `physics_reconciliation.pivot_eia930_tables_to_series`.

    Returns:
        pd.DataFrame: adjusted series data frame.
    """
    # SC offset = UTC <-> Eastern offset
    # This issue was corrected after Dec 31, 2020
    sc_offsets = (
        raw.index.tz_convert("US/Eastern").to_series().apply(lambda s: s.utcoffset())
    )
    sc_offsets["2020-12-31 00:00:00+00":] = timedelta(0)
    # make new data so we don't mess up other data indexing
    sc_dat = raw[get_columns("SC", raw.columns)].copy()
    sc_idx = pd.DatetimeIndex(sc_dat.index + sc_offsets)  # make shifted dates
    sc_dat.index = sc_idx  # use shifted dates
    sc_dat = sc_dat[~sc_dat.index.duplicated(keep="first")]
    raw = raw.drop(columns=sc_dat.columns)
    raw = pd.concat([raw, sc_dat], axis="columns")

    # PJM data reports start of hour instead of end of hour
    ba = "PJM"
    cols = get_columns(ba, raw.columns)
    new = raw[cols].shift(1, freq="h")
    raw = raw.drop(columns=cols)
    raw = pd.concat([raw, new], axis="columns")

    # TEPC data reports start of hour instead of end of hour. This issue was corrected
    # after 2021-11-01
    ba = "TEPC"
    cols = get_columns(ba, raw.columns)
    new = raw[cols].copy()
    new.loc[raw.index < "2021-11-01 00:00:00+00", cols] = new.loc[
        raw.index < "2021-11-01 00:00:00+00", cols
    ].shift(1, freq="h")
    raw = raw.drop(columns=cols)
    raw = pd.concat([raw, new], axis="columns")

    # CISO reported start of hour data through June 13, 2022
    # June 14, 2022: end-of-hour
    # June 15, 2022: start-of-hour
    # June 16, 2022: end-of-hour
    # November 1, 2023: end-of-hour + 1h
    # December 1, 2024: lag leading to highest correlation is 6h. Could not be
    # confirmed visually.
    ba = "CISO"
    cols = get_columns(ba, raw.columns)
    new = raw[cols].copy()
    new.loc[
        (raw.index < "2022-06-14 07:00:00+00")
        | (
            (raw.index >= "2022-06-15 07:00:00+00")
            & (raw.index < "2022-06-16 07:00:00+00")
        ),
        cols,
    ] = new.loc[
        (raw.index < "2022-06-14 07:00:00+00")
        | (
            (raw.index >= "2022-06-15 07:00:00+00")
            & (raw.index < "2022-06-16 07:00:00+00")
        ),
        cols,
    ].shift(1, freq="h")
    new.loc[
        (raw.index >= "2023-11-01 00:00:00+00"),
        cols,
    ] = new.loc[
        (raw.index >= "2023-11-01 00:00:00+00"),
        cols,
    ].shift(-1, freq="h")
    raw = raw.drop(columns=cols)
    raw = pd.concat([raw, new], axis="columns")

    # AVA
    # March 11, 2024: start-of-hour
    # November 3, 2024: end of hour
    ba = "AVA"
    cols = get_columns(ba, raw.columns)
    new = raw[cols].copy()
    new.loc[
        (
            (raw.index >= "2024-03-11 07:00:00+00")
            & (raw.index < "2024-11-03 00:00:00+00")
        ),
        cols,
    ] = new.loc[
        (
            (raw.index >= "2024-03-11 07:00:00+00")
            & (raw.index < "2024-11-03 00:00:00+00")
        ),
        cols,
    ].shift(1, freq="h")
    raw = raw.drop(columns=cols)
    raw = pd.concat([raw, new], axis="columns")

    # BANC
    # March 10, 2024: start of hour
    # April 1, 2024: end-of-hour
    ba = "BANC"
    cols = get_columns(ba, raw.columns)
    new = raw[cols].copy()
    new.loc[
        (
            (raw.index >= "2024-03-10 13:00:00+00")
            & (raw.index < "2024-04-01 08:00:00+00")
        ),
        cols,
    ] = new.loc[
        (
            (raw.index >= "2024-03-10 13:00:00+00")
            & (raw.index < "2024-04-01 08:00:00+00")
        ),
        cols,
    ].shift(1, freq="h")
    raw = raw.drop(columns=cols)
    raw = pd.concat([raw, new], axis="columns")

    # Interchange sign. Do before we change interchange time for PJM, because
    # identification of sign shift is based on raw data
    cols = get_int_columns(
        "PJM", raw.columns, ["CPLE", "CPLW", "DUK", "LGEE", "MISO", "NYIS", "TVA"]
    )
    raw.loc[raw.index < "2019-10-31 04:00:00+00", cols] = (
        raw.loc[raw.index < "2019-10-31 04:00:00+00", cols] * -1
    )

    # Interchange AZPS - SRP is wonky before 6/1/2020 7:00 UTC. Use SRP - AZPS (inverted)
    azps_srp = get_int_columns("AZPS", raw.columns, ["SRP"])
    srp_azps = get_int_columns("SRP", raw.columns, ["AZPS"])
    raw.loc[:"2020-06-01 07:00:00+00", azps_srp[0]] = raw.loc[
        :"2020-06-01 07:00:00+00", srp_azps[0]
    ] * (-1)
    # Update total interchange
    all_cols = [
        c for c in get_int_columns("AZPS", raw.columns) if c[0] != TOTAL_INTERCHANGE
    ]
    total_col = (TOTAL_INTERCHANGE, "AZPS", "")
    raw.loc[:"2020-06-01 07:00:00+00", total_col] = raw.loc[
        :"2020-06-01 07:00:00+00", all_cols
    ].sum(axis=1)

    # Interchange TEPC is uniformly lagged until 10-25-2021
    cols = get_int_columns("TEPC", raw.columns)
    new = raw[cols].copy()
    new.loc[raw.index < "2021-10-25 00:00:00+00", cols] = new.loc[
        raw.index < "2021-10-25 00:00:00+00", cols
    ].shift(-7, freq="h")
    raw = raw.drop(columns=cols)
    raw = pd.concat([raw, new], axis="columns")

    # Interchange IID-CISO is lagged by 4 hours in 2021
    col = get_int_columns("IID", raw.columns, ["CISO"])
    new = raw[col].copy()
    new.loc["2021-01-01 08:00:00+00:00":"2022-01-01 07:00:00+00:00", col] = new.loc[
        "2021-01-01 08:00:00+00:00":"2022-01-01 07:00:00+00:00", col
    ].shift(4, freq="h")
    raw = raw.drop(columns=col)
    raw = pd.concat([raw, new], axis="columns")

    # Interchange PJM is lagged differently across DST boundary
    is_dst = raw.index.tz_convert("US/Eastern").to_series().apply(
        lambda s: s.utcoffset()
    ) == timedelta(hours=-4)
    pjm_offset = [
        timedelta(hours=-3) if is_d else timedelta(hours=-4) for is_d in is_dst
    ]

    # make new data so we don't mess up other data indexing
    pjm_dat = raw[
        get_int_columns(
            "PJM",
            raw.columns,
            ["CPLE", "CPLW", "DUK", "LGEE", "MISO", "NYIS", "TVA", "ALL"],
        )
    ].copy()
    # make shifted dates
    pjm_idx = pd.DatetimeIndex(pjm_dat.index + pd.Series(pjm_offset))
    pjm_dat.index = pjm_idx  # use shifted dates
    # delete duplicates
    pjm_dat = pjm_dat[~pjm_dat.index.duplicated(keep="first")]
    # exchange old rows with new
    raw = raw.drop(columns=pjm_dat.columns)
    raw = pd.concat([raw, pjm_dat], axis="columns")

    # Shift all -1 hour to make start-of-hour
    raw = raw.shift(-1, freq="h")
    raw.index.name = "datetime_utc"
    return raw


def get_columns(ba: str, columns: pd.MultiIndex) -> list[tuple]:
    """Returns the keys of a BA's generation by energy source, net generation and
    demand series.

    Args:
        ba (str): BA code.
        columns (pd.MultiIndex): columns of a series data frame.

    Returns:
        list[tuple]: keys of the BA's generation by energy source series that exist,
            followed by the keys of its net generation and demand series.
    """
    cols = [
        (NET_GENERATION_BY_ENERGY_SOURCE, ba, energy_source)
        for energy_source in ENERGY_SOURCES
        if (NET_GENERATION_BY_ENERGY_SOURCE, ba, energy_source) in columns
    ]
    cols.append((NET_GENERATION, ba, ""))
    cols.append((DEMAND, ba, ""))
    return cols


def get_int_columns(
    ba1: str, columns: pd.MultiIndex, ba2: list[str] | None = None
) -> list[tuple]:
    """Returns the keys of a BA's interchange series.

    Args:
        ba1 (str): BA code.
        columns (pd.MultiIndex): columns of a series data frame.
        ba2 (list[str] | None, optional): adjacent BAs to include. "ALL" means total
            interchange. If None or empty, include all adjacent BAs and total
            interchange. Defaults to None.

    Returns:
        list[tuple]: keys of the interchange series that exist, followed by the key of
            total interchange if included.
    """
    # Looking for everyone, including ALL
    if not ba2:
        ba2 = [key[2] for key in columns if key[0] == INTERCHANGE and key[1] == ba1] + [
            "ALL"
        ]

    cols = [(INTERCHANGE, ba1, ba) for ba in ba2 if (INTERCHANGE, ba1, ba) in columns]
    if "ALL" in ba2:
        if (TOTAL_INTERCHANGE, ba1, "") in columns:  # CFE lacks "ALL" interchange
            cols.append((TOTAL_INTERCHANGE, ba1, ""))
    return cols
