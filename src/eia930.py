import pandas as pd
import re
from datetime import timedelta
import os
from os.path import join

import load_data
from column_checks import get_dtypes
from filepaths import top_folder, downloads_folder, outputs_folder, manual_folder

# Tell gridemissions where to find config before we load gridemissions
os.environ["GRIDEMISSIONS_CONFIG_FILE_PATH"] = top_folder("config/gridemissions.json")

from gridemissions.workflows import make_dataset


def convert_balance_file_to_gridemissions_format(year: int, small: bool = False):
    """Converts downloaded EIA-930 Balance files to gridemissions format."""
    files = [
        downloads_folder() + "eia930/EIA930_{}_{}_Jul_Dec.csv",
        downloads_folder() + "eia930/EIA930_{}_{}_Jan_Jun.csv",
        downloads_folder() + "eia930/EIA930_{}_{}_Jul_Dec.csv",
    ]

    years = [year - 1, year, year]

    if small:
        files = [downloads_folder() + "eia930/EIA930_{}_{}_Jan_Jun.csv"]
        years = [year]

    name_map = {
        "Total Interchange (MW)": "EBA.{}-ALL.TI.H",
        "Interchange (MW)": "EBA.{}-{}.ID.H",
        "Demand (MW) (Adjusted)": "EBA.{}-ALL.D.H",
        "Net Generation (MW) (Adjusted)": "EBA.{}-ALL.NG.H",
        "Net Generation (MW) from Coal": "EBA.{}-ALL.NG.COL.H",
        "Net Generation (MW) from Natural Gas": "EBA.{}-ALL.NG.NG.H",
        "Net Generation (MW) from Nuclear": "EBA.{}-ALL.NG.NUC.H",
        "Net Generation (MW) from All Petroleum Products": "EBA.{}-ALL.NG.OIL.H",
        "Net Generation (MW) from Hydropower and Pumped Storage": "EBA.{}-ALL.NG.WAT.H",
        "Net Generation (MW) from Solar": "EBA.{}-ALL.NG.SUN.H",
        "Net Generation (MW) from Wind": "EBA.{}-ALL.NG.WND.H",
        "Net Generation (MW) from Other Fuel Sources": "EBA.{}-ALL.NG.OTH.H",
        "Net Generation (MW) from Unknown Fuel Sources": "EBA.{}-ALL.NG.UNK.H",
    }

    out = pd.DataFrame()
    for i, file in enumerate(files):
        dat_file = file.format("BALANCE", years[i])
        int_file = file.format("INTERCHANGE", years[i])

        # Format balance files in series format (for gridemissions)
        dat = pd.read_csv(
            dat_file,
            usecols=[
                "Balancing Authority",
                "UTC Time at End of Hour",
                "Total Interchange (MW)",
                "Demand (MW) (Adjusted)",
                "Net Generation (MW) (Adjusted)",
                "Net Generation (MW) from Coal",
                "Net Generation (MW) from Natural Gas",
                "Net Generation (MW) from Nuclear",
                "Net Generation (MW) from All Petroleum Products",
                "Net Generation (MW) from Hydropower and Pumped Storage",
                "Net Generation (MW) from Solar",
                "Net Generation (MW) from Wind",
                "Net Generation (MW) from Other Fuel Sources",
                "Net Generation (MW) from Unknown Fuel Sources",
            ],
            parse_dates=["UTC Time at End of Hour"],
            thousands=",",
        )
        # Wide to long
        dat = dat.melt(id_vars=["Balancing Authority", "UTC Time at End of Hour"])
        # Find series name
        dat["column"] = dat.apply(
            lambda x: name_map[x.variable].format(x["Balancing Authority"]),
            axis="columns",
        )
        # Long to wide
        dat = dat[["UTC Time at End of Hour", "value", "column"]].pivot(
            index="UTC Time at End of Hour", columns="column", values="value"
        )

        # Now for interchange
        int = pd.read_csv(
            int_file,
            usecols=[
                "Balancing Authority",
                "Directly Interconnected Balancing Authority",
                "Interchange (MW)",
                "UTC Time at End of Hour",
            ],
            parse_dates=["UTC Time at End of Hour"],
            thousands=",",
        )
        int["column"] = int.apply(
            lambda x: name_map["Interchange (MW)"].format(
                x["Balancing Authority"],
                x["Directly Interconnected Balancing Authority"],
            ),
            axis="columns",
        )
        int = int[["UTC Time at End of Hour", "column", "Interchange (MW)"]].pivot(
            index="UTC Time at End of Hour", columns="column", values="Interchange (MW)"
        )

        # Combine
        dat = pd.concat([dat, int], axis="columns")
        out = pd.concat([out, dat], axis="index")

    out.index = out.index.tz_localize("UTC")
    # Balance files are all inclusive, so hours at boundaries (July 1, Jan 1) are duplicated.
    # Drop those duplicate rows
    out = out[~out.index.duplicated(keep="first")]

    return out


def clean_930(year: int, small: bool = False, path_prefix: str = ""):
    """
        Scrape and process EIA data.

    Arguments:
        `year`: Year to process. Prior years, downloaded from chalendar-hosted files, are used for rolling cleaning

    """

    data_folder = outputs_folder(f"{path_prefix}/eia930/")

    # Format raw file
    df = convert_balance_file_to_gridemissions_format(year, small=small)
    raw_file = data_folder + "eia930_unadjusted_raw.csv"
    df.to_csv(raw_file)

    # if not small, scrape 2 months before start of year for rolling window cleaning
    start = f"{year}0101T00Z" if small else f"{year-1}1001T00Z"
    # Scrape 1 week if small, else 1 year (plus one day for timezone flexibility)
    end = f"{year}0107T23Z" if small else f"{year+1}0101T23Z"
    if small:
        df = df.loc[start:end]  # Don't worry about processing everything

    # Adjust
    print("    Adjusting EIA-930 time stamps")
    df = manual_930_adjust(df)
    df.to_csv(
        join(data_folder, "eia930_raw.csv")
    )  # Will be read by gridemissions workflow

    # Run cleaning
    print("    Running physics-based data cleaning")
    make_dataset(
        start,
        end,
        file_name="eia930",
        tmp_folder=data_folder,
        folder_hist=data_folder,
        scrape=False,
        add_ca_fuels=False,
        calc_consumed=False,
    )


def reformat_chalendar(raw):
    """
        reformat_chalendar
    Reformat wide-format data (one row per time stamp) from Chalendar
    to long (one row per data point)
    Drop all columns that are not fuel-specific generation
    """
    # where we have variable (NG = net generation) and fuel type
    target_cols = [c for c in raw.columns if len(c.split(".")) == 5]
    print("Filtering")
    cleaned = (
        raw.loc[:, target_cols]
        .melt(ignore_index=False, value_name="generation", var_name="variable")
        .reset_index()
    )
    print("Expanding cols")
    cleaned[["dtype", "BA", "other BA", "var", "fuel", "interval"]] = cleaned[
        "variable"
    ].str.split(r"[.-]", expand=True, regex=True)
    print("Dropping and renaming")
    cleaned = cleaned.drop(columns=["dtype", "var", "interval", "other BA"])
    cleaned = cleaned.rename(columns={"index": "datetime_utc"})
    return cleaned


def load_chalendar(fname: str, year: int):
    raw = pd.read_csv(fname, index_col=0, parse_dates=True)
    raw = raw[raw.index.year == year]
    return reformat_chalendar(raw)


def load_chalendar_for_pipeline(cleaned_data_filepath, year):
    """
    Loads and formats cleaned hourly net generation data
    for use in the data pipeline
    """
    # read the data, only keeping net generation columns
    data = pd.read_csv(
        cleaned_data_filepath, index_col=0, parse_dates=True, dtype=get_dtypes()
    ).filter(like="-ALL.NG.")

    # name the index
    data.index = data.index.rename("datetime_utc")

    # only keep data for the single year we are interested, plus or minus one
    # day, to account for conversion to local time later
    data = data.loc[
        (data.index >= f"{year-1}-12-31") & (data.index < f"{year+1}-01-02"), :
    ]

    # remove columns for total net generation
    data = data.loc[
        :,
        [
            col
            for col in data.columns
            if col not in data.filter(like="-ALL.NG.H").columns
        ],
    ]

    # melt the data into long format
    data = data.reset_index().melt(
        id_vars="datetime_utc", value_name="net_generation_mwh_930"
    )

    # create columns for ba_code and fuel category
    data[["ba_code", "fuel_category_eia930"]] = data["variable"].str.split(
        r"[.-]", expand=True, regex=True
    )[[1, 4]]

    # drop BAs not located in the United States
    ba_ref = pd.read_csv(manual_folder("ba_reference.csv"))
    foreign_bas = list(ba_ref.loc[ba_ref["us_ba"] == "No", "ba_code"])
    data = data[~data["ba_code"].isin(foreign_bas)]

    data["datetime_local"] = ""
    for ba in list(data["ba_code"].unique()):
        data.loc[data.ba_code == ba, "datetime_local"] = (
            data.loc[data.ba_code == ba, "datetime_utc"]
            .dt.tz_convert(load_data.ba_timezone(ba=ba, type="local"))
            .astype(str)
        )

    # create a report date column
    data["report_date"] = data["datetime_local"].str[:7]
    data["report_date"] = pd.to_datetime(data["report_date"])

    # rename the fuel categories using format in
    # data/manual/energy_source_groups
    fuel_categories = {
        "COL": "coal",
        "NG": "natural_gas",
        "OTH": "other",
        "WAT": "hydro",
        "WND": "wind",
        "SUN": "solar",
        "NUC": "nuclear",
        "OIL": "petroleum",
        "BIO": "biomass",
        "GEO": "geothermal",
    }

    data["fuel_category_eia930"] = data["fuel_category_eia930"].map(fuel_categories)

    # reorder the columns and remove the variable column
    data = data[
        [
            "ba_code",
            "fuel_category_eia930",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "net_generation_mwh_930",
        ]
    ]

    return data


def remove_imputed_ones(eia930_data):

    filter = [eia930_data["net_generation_mwh_930"].abs() < 1.5]

    # replace all 1.0 values with zero
    print(f"  replacing {sum(filter)} imputed 1 values with 0")
    eia930_data.loc[filter, "net_generation_mwh_930"] = 0

    return eia930_data


def remove_months_with_zero_data(eia930_data):
    # remove data where the entire month is zero
    zero_data = (
        eia930_data.groupby(["ba_code", "fuel_category_eia930", "report_date"])
        .sum()
        .reset_index()
    )

    zero_data = zero_data[zero_data["net_generation_mwh_930"] == 0].drop(
        columns="net_generation_mwh_930"
    )

    # filter these ba-fuel-months out of the eia930 data
    eia930_data = eia930_data.merge(
        zero_data,
        how="outer",
        on=["ba_code", "fuel_category_eia930", "report_date"],
        indicator="zero_filter",
        validate="m:1",
    )
    eia930_data = eia930_data[eia930_data["zero_filter"] == "left_only"].drop(
        columns="zero_filter"
    )

    return eia930_data


###########################################################
# Code for adjusting 930 data in gridemissions format
#
###########################################################


def get_columns(ba: str, columns):
    GEN_ID = "EBA.{}-ALL.NG.H"
    GEN_TYPE_ID = "EBA.{}-ALL.NG.{}.H"
    DEM_ID = "EBA.{}-ALL.D.H"
    SRC = ["COL", "NG", "NUC", "OIL", "OTH", "SUN", "UNK", "WAT", "WND", "GEO", "BIO"]

    cols = [
        GEN_TYPE_ID.format(ba, f) for f in SRC if GEN_TYPE_ID.format(ba, f) in columns
    ]
    cols.append(GEN_ID.format(ba))
    cols.append(DEM_ID.format(ba))
    return cols


def get_int_columns(ba1: str, columns, ba2: list = []):
    INTER_ID = "EBA.{}-{}.ID.H"
    IT_ID = "EBA.{}-ALL.TI.H"

    # Looking for everyone, including ALL
    if ba2 == []:
        other_cols = [
            c
            for c in columns
            if re.split(r"[-.]", c)[1] == ba1 and re.split(r"[-.]", c)[2] != "ALL"
        ]
        ba2 = [re.split(r"[-.]", c)[2] for c in other_cols]
        ba2.append("ALL")

    cols = [
        INTER_ID.format(ba1, ba) for ba in ba2 if (INTER_ID.format(ba1, ba) in columns)
    ]
    if "ALL" in ba2:
        if IT_ID.format(ba1) in columns:  # CFE lacks "ALL" interchange
            cols.append(IT_ID.format(ba1))
    return cols


def manual_930_adjust(raw: pd.DataFrame):
    """
        manual_930_adjust

    Adjusts time stamps in 930 data. Assumes dataframe with timestamp index and
    one column per series,
    where column names correspond to EIA series IDs:
        "EBA.%s-ALL.D.H",  # Demand
        "EBA.%s-ALL.NG.H",  # Generation
        "EBA.%s-ALL.NG.%s.H",  # Generation by fuel type
        "EBA.%s-ALL.TI.H",  # Total Interchange
        "EBA.%s-%s.ID.H",  # Interchange

    Adjustment Steps:

    - Make all end-of-hour
        - Generation
            - PJM: + 1 hour
            - CISO: + 1 hour
            - TEPC: + 1 hour
            - SC: -4 hours during daylight savings hours; -5 hours during
                standard hours (this happens to = the Eastern <-> UTC offset)
                Fixed in BALANCE files starting Jan 1, 2021
        - Interchange
            - PJM: + 4 hours
            - TEPC: + 7 hours
            - CFE:  -11 hours
        - Interchange sign
            - PJM-{CPLE, CPLW, DUK, LGEE, MISO, NYIS, TVA} before
                    Oct 31, 2019, 4:00 UTC
                - this is all interchange partners except OVEC, and excluding
                    total interchange
        - Interchange mysterious
            - AZPS - SRP flips gradually in Nov 2019, then abruptly back in June 2020.
                throughout, SRP - AZPS remains constant around 3000 lb imported to AZPS from SRP
                We assume SRP - AZPS is correct, and assign AZPS - SRP to be the inverse
    - Make all start-of-hour
        - Generation
            - - 1 hour
    """
    # SC offset = UTC <-> Eastern offset
    sc_offsets = (
        raw.index.tz_convert("US/Eastern").to_series().apply(lambda s: s.utcoffset())
    )
    # After Dec 31, 2020, the offset is 0
    sc_offsets["2020-12-31T00:00":] = timedelta(0)
    # make new data so we don't mess up other data indexing
    sc_dat = raw[get_columns("SC", raw.columns)].copy()
    sc_idx = pd.DatetimeIndex(sc_dat.index + sc_offsets)  # make shifted dates
    sc_dat.index = sc_idx  # use shifted dates
    sc_dat = sc_dat[~sc_dat.index.duplicated(keep="first")]
    # exchange old rows with new
    raw = raw.drop(columns=sc_dat.columns)
    raw = pd.concat([raw, sc_dat], axis="columns")

    # PJM, CISO, TEPC: shift by one hour
    for ba in ["PJM", "CISO", "TEPC"]:
        cols = get_columns(ba, raw.columns)
        new = raw[cols].shift(1, freq="H")
        raw = raw.drop(columns=cols)
        raw = pd.concat([raw, new], axis="columns")

    # Interchange sign. Do before we change interchange time for PJM, because
    # identification of sign shift is based on raw data
    cols = get_int_columns(
        "PJM", raw.columns, ["CPLE", "CPLW", "DUK", "LGEE", "MISO", "NYIS", "TVA"]
    )
    raw.loc[raw.index < "2019-10-31T04", cols] = (
        raw.loc[raw.index < "2019-10-31T04", cols] * -1
    )

    # Interchange AZPS - SRP is wonky before 6/1/2020 7:00 UTC. Use SRP - AZPS (inverted)
    azps_srp = get_int_columns("AZPS", raw.columns, ["SRP"])
    srp_azps = get_int_columns("SRP", raw.columns, ["AZPS"])
    replacement = (raw.loc[:, srp_azps] * (-1)).rename(
        columns={srp_azps[0]: azps_srp[0]}  # rename so Pandas will do the right thing
    )
    raw.loc[:"2020-06-01T07:00+00", azps_srp] = replacement[:"2020-06-01T07:00+00"]
    # Update total interchange
    all_cols = [c for c in get_int_columns("AZPS", raw.columns) if "ALL" not in c]
    total_col = "EBA.AZPS-ALL.TI.H"
    raw.loc[:"2020-06-01T07:00+00", total_col] = raw.loc[
        :"2020-06-01T07:00+00", all_cols
    ].sum(axis=1)

    # Interchange TEPC is uniformly lagged
    cols = get_int_columns("TEPC", raw.columns)
    new = raw[cols].shift(-7, freq="H")
    raw = raw.drop(columns=cols)
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
    return raw.shift(-1, freq="H")
