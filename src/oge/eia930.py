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
    TOTAL_INTERCHANGE,
    EIA930Tables,
    melt_series_to_eia930_tables,
    pivot_eia930_tables_to_series,
    reconcile_eia930,
    write_eia930_tables,
)

logger = get_logger(__name__)


# Dictionary mapping the new energy sources in EIA-930 to existing ones
new_to_existing_energy_sources = {
    "coal": "coal",
    "gas": "gas",
    "hydro": "hydro",
    "hydro_excluding_pumped_storage": "hydro",
    "pumped_storage": "hydro",
    "nuclear": "nuclear",
    "oil": "oil",
    "battery_storage": "other",
    "geothermal": "other",
    "other": "other",
    "other_energy_storage": "other",
    "unknown_energy_storage": "other",
    "solar": "solar",
    "solar_w_integrated_battery_storage": "solar",
    "solar_wo_integrated_battery_storage": "solar",
    "unknown": "unknown",
    "wind": "wind",
    "wind_w_integrated_battery_storage": "wind",
    "wind_wo_integrated_battery_storage": "wind",
}

# Dictionary mapping EIA-930 energy sources to the fuel categories in
# reference_tables/energy_source_groups.csv. "unknown" is not mapped.
energy_source_to_fuel_category_eia930 = {
    "coal": "coal",
    "gas": "natural_gas",
    "other": "other",
    "hydro": "hydro",
    "wind": "wind",
    "solar": "solar",
    "nuclear": "nuclear",
    "oil": "petroleum",
}


def clean_930(year: int, path_prefix: str = "") -> EIA930Tables:
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
        EIA930Tables: reconciled EIA-930 data.
    """
    data_folder = outputs_folder(f"{path_prefix}/eia930/")

    tables = load_eia930_balance_data(year)
    write_eia930_tables(tables, data_folder, "unadjusted")

    logger.info("Adjusting EIA-930 time stamps")
    tables = manual_930_adjust(tables)
    write_eia930_tables(tables, data_folder, "raw")

    logger.info("Running physics-based data cleaning")
    return reconcile_eia930(tables, data_folder)


def load_eia930_balance_data(year: int) -> EIA930Tables:
    """Loads EIA-930 operations, generation by energy source and interchange from PUDL.

    Loads data from October of the previous year through January 2 of the next year,
    since the rolling window cleaning needs data before and after the year. Generation
    by energy source is mapped to the energy sources used before EIA-930 introduced new
    ones (see `new_to_existing_energy_sources`), since BAs started to use the new energy
    sources at different times.

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

    # TODO: consider individual energy sources (e.g. geothermal, battery storage,
    # etc.) and BAs reporting new energy sources at different times.
    # note that sum() turns groups with only blank values into zero, as in the
    # original gridemissions-format code
    net_generation_by_energy_source = (
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
        )
        .astype(
            {"generation_energy_source": "str"}
        )  # to avoid a FutureWarning when using replace on CategoricalDtype
        .replace({"generation_energy_source": new_to_existing_energy_sources})
        .groupby(
            [
                "datetime_utc",
                "balancing_authority_code_eia",
                "generation_energy_source",
            ],
            dropna=True,  # the keys are never blank in PUDL
        )["net_generation_adjusted_mwh"]
        .sum()
        .reset_index()
        .rename(columns={"net_generation_adjusted_mwh": "net_generation_mwh"})
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
    later, and drops BAs outside of the US.

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
    data = data[~data["ba_code"].isin(foreign_bas)].copy()

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

    # rename the energy sources using the fuel categories in
    # reference_tables/energy_source_groups.csv
    data["fuel_category_eia930"] = data["generation_energy_source"].map(
        energy_source_to_fuel_category_eia930
    )

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


def remove_imputed_ones(eia930_data: pd.DataFrame) -> pd.DataFrame:
    """Remove ones in data frame.

    Args:
        eia930_data (pd.DataFrame): EIA-930 table

    Returns:
        pd.DataFrame: input data frame with ones replace by zeros
    """
    filter = eia930_data["net_generation_mwh_930"].abs() < 1.5

    # replace all 1.0 values with zero
    logger.info(f"Replacing {sum(filter)} imputed 1 values with 0")
    eia930_data.loc[filter, "net_generation_mwh_930"] = 0

    return eia930_data


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
