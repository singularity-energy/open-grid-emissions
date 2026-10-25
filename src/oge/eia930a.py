"""Cleans EIA-930A, the annual generator inventory of each BA that reports to EIA-930.

Each BA lists the generators reflected in its EIA-930 net generation data as of the end
of the year (Schedule 2), the generators it plans to add in the following year
(Schedule 3), and the generators with which it has a pseudo-tie or dynamic scheduling
arrangement (Schedule 4), whether or not they are included in Schedule 2. EIA matches
each generator to EIA-860. Plant and generator IDs are published as validated EIA-860
IDs where they match, and as submitted by the BA otherwise.

EIA-930A is not yet in PUDL (catalyst-cooperative/pudl#5218), so it is downloaded from
the EIA website (see `download_data.download_raw_eia930a`).
"""

import os
import re
from typing import NamedTuple

import numpy as np
import pandas as pd

import oge.load_data as load_data
from oge.download_data import get_eia930a_year
from oge.filepaths import outputs_folder
from oge.logging_util import get_logger

logger = get_logger(__name__)


class EIA930ATables(NamedTuple):
    """Cleaned EIA-930A schedules, with one row per generator listed by each BA."""

    schedule_2: pd.DataFrame
    schedule_3: pd.DataFrame
    schedule_4: pd.DataFrame


# BA codes in the Schedule 4 partner column that are not EIA-930 BA codes
PARTNER_BA_ALIASES = {
    "APS": "AZPS",
    "BPA": "BPAT",
    "CAISO": "CISO",
    "CHLN": "CHPD",
    "LADWP": "LDWP",
    "NVE": "NEVP",
    "NVEP": "NEVP",
    "TWPR": "TPWR",
}
# values in the Schedule 4 partner column that are not BAs. These are kept in the
# `partner_note` column. WAPA operates several BAs (WACM, WALC and WAUW), so it can't be
# resolved to one of them.
PARTNER_NOTES = ["Not Generated", "WAPA"]
# CISO wrote the type of arrangement in the partner column instead of the partner BA,
# and left "reported in Schedule 2" blank
PARTNER_ARRANGEMENT_TYPES = {
    "Dynamic import": "dynamic_schedule",
    "Pseudo tie": "pseudo_tie",
}

# values that some BAs submitted in place of a generator ID
PLACEHOLDER_ID_PATTERN = r"^(?:DOES NOT EXIST|NOT GENERATED|TBD)"


def clean_eia930a(year: int, path_prefix: str) -> EIA930ATables:
    """Loads and cleans the EIA-930A data to use for a year, and writes the cleaned
    schedules to `outputs/{path_prefix}/eia930/`.

    Args:
        year (int): a four-digit year.
        path_prefix (str): path to the folder enclosing the eia930 folder.

    Returns:
        EIA930ATables: cleaned EIA-930A data, with the year of the EIA-930A data used in
            the "report_year" column (see `download_data.get_eia930a_year`).
    """
    report_year = get_eia930a_year(year)
    if report_year != year:
        logger.warning(
            f"EIA-930A data is not available for {year}. Using {report_year} instead."
        )
    tables = clean_raw_eia930a(load_data.load_raw_eia930a(report_year), report_year)
    write_eia930a_tables(tables, outputs_folder(f"{path_prefix}/eia930/"))
    return tables


def clean_raw_eia930a(
    schedules: dict[str, pd.DataFrame], report_year: int
) -> EIA930ATables:
    """Cleans the EIA-930A schedules of a year.

    Args:
        schedules (dict[str, pd.DataFrame]): schedules from `load_data.load_raw_eia930a`.
        report_year (int): the year of the EIA-930A data.

    Returns:
        EIA930ATables: cleaned EIA-930A data.
    """
    schedules = {
        name: clean_eia930a_generators(df, report_year)
        for name, df in schedules.items()
    }
    schedules["schedule_4"] = clean_eia930a_schedule_4(
        schedules["schedule_4"], schedules["schedule_2"]
    )
    return EIA930ATables(**schedules)


def clean_eia930a_generators(df: pd.DataFrame, report_year: int) -> pd.DataFrame:
    """Cleans the columns that all EIA-930A schedules have in common.

    Plant and generator IDs are the validated EIA-860 IDs where EIA matched the
    generator to EIA-860 (`validated_eia860`), and the IDs submitted by the BA
    otherwise. Submitted IDs that aren't IDs (e.g. "Does Not Exist" or "maybe 55086")
    are missing. Generator IDs are upper case, so they should be matched with other
    generator IDs in upper case. The capacity, state and county are from EIA-860 where
    available, and as submitted otherwise.

    Generators that a BA submitted as one entry (e.g. "1, 2, 3" or "GEN1/2") are kept as
    one row, since their generator ID can't be split reliably.

    Args:
        df (pd.DataFrame): an EIA-930A schedule, from `load_data.load_raw_eia930a`.
        report_year (int): the year of the EIA-930A data.

    Returns:
        pd.DataFrame: the cleaned schedule.
    """
    df = (
        df.dropna(how="all").apply(lambda column: column.str.strip()).replace("", pd.NA)
    )
    df.insert(0, "report_year", report_year)
    df["validated_eia860"] = df["plant_id_eia"].notna() & df["generator_id"].notna()
    df["plant_id_eia"] = (
        pd.to_numeric(df["plant_id_eia"], errors="coerce")
        .fillna(pd.to_numeric(df["plant_id_as_submitted"], errors="coerce"))
        .astype("Int64")
    )
    generator_id = (
        df["generator_id"]
        .fillna(df["generator_id_as_submitted"])
        .astype("string")
        .str.upper()
    )
    df["generator_id"] = generator_id.mask(
        generator_id.str.contains(PLACEHOLDER_ID_PATTERN, na=False)
    )
    df["capacity_mw"] = pd.to_numeric(df["capacity_mw"], errors="coerce").fillna(
        pd.to_numeric(df["capacity_mw_as_submitted"], errors="coerce")
    )
    df["state"] = df["state"].fillna(df["state_as_submitted"])
    df["county"] = df["county"].fillna(df["county_as_submitted"])
    for column in ["latitude", "longitude"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    for column in df.columns.intersection(
        [
            "operating_year",
            "operating_month",
            "planned_operating_year",
            "planned_operating_month",
        ]
    ):
        df[column] = pd.to_numeric(df[column], errors="coerce").astype("Int64")
    return df.drop(
        columns=[
            "capacity_mw_as_submitted",
            "state_as_submitted",
            "county_as_submitted",
        ]
    )


def clean_eia930a_schedule_4(
    schedule_4: pd.DataFrame, schedule_2: pd.DataFrame
) -> pd.DataFrame:
    """Cleans the Schedule 4 columns that describe each pseudo-tie or dynamic
    scheduling arrangement.

    The partner BAs are split into `partner_ba_codes` (comma-separated EIA-930 BA
    codes, with aliases such as "CAISO" replaced by their EIA-930 code),
    `arrangement_type` ("dynamic_schedule" or "pseudo_tie", where the BA wrote it in
    place of the partner) and `partner_note` (other values that aren't BAs, e.g.
    "WAPA").

    `reported_in_schedule_2` is True for "Yes" and False for "No". Where it's blank,
    it's True if the filing BA lists the generator in its Schedule 2, and False
    otherwise. In 2023 and 2024, CISO left it blank for the generators it imports
    through a dynamic schedule or pseudo-tie, none of which it lists in Schedule 2.

    Args:
        schedule_4 (pd.DataFrame): Schedule 4, from `clean_eia930a_generators`.
        schedule_2 (pd.DataFrame): Schedule 2 of the same year, from
            `clean_eia930a_generators`.

    Returns:
        pd.DataFrame: the cleaned Schedule 4.
    """
    schedule_4 = schedule_4.copy()
    ba_codes = set(load_data.load_ba_reference()["ba_code"])
    partners = schedule_4["partner_bas_as_submitted"].apply(
        parse_partner_bas, ba_codes=ba_codes
    )
    schedule_4["partner_ba_codes"] = partners.str[0]
    schedule_4["arrangement_type"] = partners.str[1]
    schedule_4["partner_note"] = partners.str[2]

    answer = schedule_4["reported_in_schedule_2"].str.lower()
    if not answer.dropna().isin(["yes", "no"]).all():
        raise ValueError(
            f"Unexpected values in Schedule 4 reported_in_schedule_2: {set(answer)}"
        )
    keys = ["balancing_authority_code_eia", "plant_id_eia", "generator_id"]
    listed = pd.MultiIndex.from_frame(schedule_4[keys]).isin(
        pd.MultiIndex.from_frame(schedule_2[keys].dropna())
    )
    schedule_4["reported_in_schedule_2"] = np.where(
        answer.isna(), listed, answer.eq("yes")
    )
    return schedule_4


def parse_partner_bas(
    partner_bas: str | None, ba_codes: set[str]
) -> tuple[str | None, str | None, str | None]:
    """Parses the "BAs having a Pseudo-Tie or Dynamic Relationship" column of EIA-930A
    Schedule 4.

    Args:
        partner_bas (str | None): the value as submitted, e.g. "CAISO, WALC".
        ba_codes (set[str]): valid BA codes, from `reference_tables/ba_reference.csv`.

    Returns:
        tuple[str | None, str | None, str | None]: the partner BA codes, sorted and
            separated by commas; the type of arrangement, if the BA wrote it in place of
            the partner; and any other values that aren't BA codes, separated by commas.
            Each is None if there is nothing to return.
    """
    if pd.isna(partner_bas):
        return None, None, None
    if partner_bas in PARTNER_ARRANGEMENT_TYPES:
        return None, PARTNER_ARRANGEMENT_TYPES[partner_bas], None
    values = [v.strip() for v in re.split(r"[,;]", partner_bas) if v.strip()]
    partner_ba_codes = {PARTNER_BA_ALIASES.get(v, v) for v in values} & ba_codes
    notes = [v for v in values if PARTNER_BA_ALIASES.get(v, v) not in ba_codes]
    for note in set(notes) - set(PARTNER_NOTES):
        logger.warning(f"Unrecognized BA '{note}' in EIA-930A Schedule 4 partner BAs")
    return ",".join(sorted(partner_ba_codes)) or None, None, ",".join(notes) or None


def write_eia930a_tables(tables: EIA930ATables, output_folder: str):
    """Writes EIA-930A tables to parquet files named "eia930a_{schedule}.parquet".

    Args:
        tables (EIA930ATables): cleaned EIA-930A data.
        output_folder (str): folder to write the files to.
    """
    os.makedirs(output_folder, exist_ok=True)
    for schedule, table in tables._asdict().items():
        table.to_parquet(
            os.path.join(output_folder, f"eia930a_{schedule}.parquet"), index=False
        )


def read_eia930a_tables(folder: str) -> EIA930ATables:
    """Reads EIA-930A tables written by `write_eia930a_tables`.

    Args:
        folder (str): folder the files were written to.

    Returns:
        EIA930ATables: cleaned EIA-930A data.
    """
    return EIA930ATables(
        **{
            schedule: pd.read_parquet(
                os.path.join(folder, f"eia930a_{schedule}.parquet")
            )
            for schedule in EIA930ATables._fields
        }
    )
