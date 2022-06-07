import pandas as pd
import re
from datetime import timedelta
import src.data_cleaning as data_cleaning
import src.load_data as load_data

# Map from 923 fuel types (added to cems data in data_pipeline)
# to 930 fuel types
# Based on EIA 923 instructions:
# https://www.eia.gov/survey/form/eia_923/instructions.pdf
fuel_code_map = {
    "NG": "NG",
    "BIT": "COL",  # bituminous
    "SUB": "COL",  # subbituminous
    "BLQ": "OTH",  # black liquor
    "RC": "COL",  # refined coal
    "WDS": "OTH",  # wood / wood waste solids
    "SUN": "SUN",  # solar. TODO: why is this in cems???
    "KER": "OIL",  # kerosene
    "RFO": "OIL",  # residual fuel oil
    "JF": "OIL",  # jet fuel
    "DFO": "OIL",  # distillate fuel oil
    "SGC": "NG",  # synthetic gas derived from coal ???
    "LIG": "COL",  # lignite coal
    "OG": "NG",  # other gas
    "PC": "COL",  # petroleum coke
    "BFG": "NG",  # blast furnace gas
    "WC": "COL",
}  # waste coal


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
        .melt(ignore_index=False, value_name="generation")
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


def load_chalendar(fname: str, year: int = 2020):
    raw = pd.read_csv(fname, index_col=0, parse_dates=True)
    raw = raw[raw.index.year == year]
    return reformat_chalendar(raw)


def load_chalendar_for_pipeline(cleaned_data_filepath, year):
    """
    Loads and formats cleaned hourly net generation data
    for use in the data pipeline
    """
    # read the data, only keeping net generation columns
    data = pd.read_csv(cleaned_data_filepath, index_col=0, parse_dates=True).filter(
        like="-ALL.NG."
    )

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
    data[["ba_code", "fuel_category"]] = data["variable"].str.split(
        r"[.-]", expand=True, regex=True
    )[[1, 4]]

    # drop BAs not located in the United States
    foreign_bas = [
        "AESO",
        "BCHA",
        "CEN",
        "CFE",
        "HQT",
        "IESO",
        "MHEB",
        "SPC",
    ]
    data = data[~data["ba_code"].isin(foreign_bas)]

    # create a local datetime column
    # TODO: Convert to local prevailing time, not local standard time
    data["datetime_local"] = data["datetime_utc"]
    for ba in list(data["ba_code"].unique()):
        data.loc[data.ba_code == ba, "datetime_local"] = data.loc[
            data.ba_code == ba, "datetime_utc"
        ].dt.tz_convert(load_data.ba_timezone(ba=ba, type="local"))

    # create a report date column
    data["report_date"] = data["datetime_local"].astype(str).str[:7]
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

    data["fuel_category"] = data["fuel_category"].map(fuel_categories)

    # reorder the columns and remove the variable column
    data = data[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "net_generation_mwh_930",
        ]
    ]

    return data


###########################################################
# Code for adjusting 930 data in gridemissions format
#
###########################################################


def get_columns(ba: str, columns):
    GEN_ID = "EBA.{}-ALL.D.H"
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
        - Interchange
            - PJM: + 4 hours
            - TEPC: + 7 hours
            - CFE:  -11 hours
        - Interchange sign
            - PJM-{CPLE, CPLW, DUK, LGEE, MISO, NYIS, TVA} before
                    Oct 31, 2019, 4:00 UTC
                - this is all interchange partners except OVEC, and excluding
                    total interchange
            - PJM-OVEC, all time. Based on OVEC demand - generation, OVEC 
                should be a net exporter to PJM
                - Note: OVEC's data repeats daily starting in 2018...
    - Make all start-of-hour
        - Generation
            - - 1 hour
    """
    # SC offset = UTC <-> Eastern offset
    sc_offsets = (
        raw.index.tz_convert("US/Eastern").to_series().apply(lambda s: s.utcoffset())
    )
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
    # OVEC sign still appears to be wrong
    ovec_col = get_int_columns("PJM", raw.columns, ["OVEC"])
    raw.loc[:, ovec_col] = raw.loc[:, ovec_col] * -1

    # Interchange TEPC is uniformly lagged
    cols = get_int_columns("TEPC", raw.columns)
    new = raw[cols].shift(-7, freq="H")
    raw = raw.drop(columns=cols)
    raw = pd.concat([raw, new], axis="columns")

    # Interchange PJM->OVEC is uniformly lagged
    cols = get_int_columns("PJM", raw.columns, ["OVEC"])
    new = raw[cols].shift(-2, freq="H")
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
