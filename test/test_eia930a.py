import logging

import openpyxl
import pandas as pd
import pytest

import oge.eia930a as eia930a
import oge.load_data as load_data
from oge.eia930a import (
    clean_eia930a,
    clean_raw_eia930a,
    parse_partner_bas,
)
from oge.load_data import load_raw_eia930a

SCHEDULE_2_HEADER = [
    "BA Code",
    "Plant Name",
    "EIA Plant ID",
    "EIA Generator ID",
    "Plant ID (EIA-930A)",
    "Generator ID (EIA-930A)",
    "Electric Generator Technology",
    "Nameplate Capacity (MW) (EIA-860)",
    "State (EIA-860)",
    "County (EIA-860)",
    "Latitude (EIA-860)",
    "Longitude (EIA-860)",
    "Operation Year (EIA-860)",
    "Opeartion Month (EIA-860)",
    "Nameplate Capacity (MW) (EIA-930A)",
    "State (EIA-930A)",
    "County (EIA-930A)",
]
SCHEDULE_3_HEADER = [
    (
        column.replace("Operation", "Planned Operation").replace(
            "Opeartion", "Planned Operation"
        )
    )
    for column in SCHEDULE_2_HEADER
]
SCHEDULE_4_HEADER = (
    SCHEDULE_2_HEADER[:6]
    + [
        "Generator Reported on EIA-930A Schedule 2?",
        "BAs having a Pseudo-Tie or Dynamic Relationship",
    ]
    + SCHEDULE_2_HEADER[6:]
)


def generator_row(
    ba, plant_name, plant_id, generator_id, capacity_mw=100, submitted=False
):
    """A row of Schedule 2 or 3. If `submitted`, the generator didn't validate against
    EIA-860, so only the IDs, capacity, state and county submitted by the BA are filled.
    """
    if submitted:
        return (
            [ba, plant_name, None, None, plant_id, generator_id, "Batteries"]
            + [None] * 7
            + [capacity_mw, "WA", "Asotin"]
        )
    return [
        ba,
        plant_name,
        plant_id,
        generator_id,
        None,
        None,
        "Natural Gas Fired Combined Cycle",
        capacity_mw,
        "WA",
        "Grant",
        47.6,
        -119.3,
        2001,
        6,
        None,
        None,
        None,
    ]


def schedule_4_row(ba, plant_name, plant_id, generator_id, reported, partners):
    row = generator_row(ba, plant_name, plant_id, generator_id)
    return row[:6] + [reported, partners] + row[6:]


def write_eia930a_file(
    folder, year, schedule_2, schedule_3, schedule_4, blank_ba_header=False
):
    """Writes the schedule sheets of an EIA-930A file, with the layout published by
    EIA."""
    workbook = openpyxl.Workbook()
    sheets = {
        "SCH 2": [[None if blank_ba_header else "BA Code"] + SCHEDULE_2_HEADER[1:]]
        + schedule_2,
        "SCH 3": [SCHEDULE_3_HEADER] + schedule_3,
        "SCH 4": [SCHEDULE_4_HEADER] + schedule_4,
    }
    workbook.remove(workbook.active)
    for name, rows in sheets.items():
        sheet = workbook.create_sheet(name)
        for row in rows:
            sheet.append(row)
    (folder / "eia930a").mkdir(exist_ok=True)
    workbook.save(folder / "eia930a" / f"EIA_930A_{year}.xlsx")


def load_eia930a(year):
    return clean_raw_eia930a(load_raw_eia930a(year), year)


@pytest.fixture
def eia930a_folder(tmp_path, monkeypatch):
    """A downloads folder with a 2023 EIA-930A file.

    - Schedule 2 has no header on its BA column, like the published 2023 file.
    - "GEN 1" at plant 1 is listed by both BAT1 and BAT2.
    - BAT1's battery didn't validate against EIA-860, so only its submitted IDs are
      available. BAT2 submitted a placeholder instead of the plant and generator ID.
    - In Schedule 4, BAT1 hosts plant 1, and AZPS and CISO take part of its output (as
      "APS, CAISO"). CISO imports plant 1's "GEN 2" and plant 3 through dynamic
      schedules, with "reported in Schedule 2" left blank. CISO lists plant 3 in its
      Schedule 2. LDWP takes part of Hoover (plant 4) from "WAPA".
    """
    monkeypatch.setattr(
        load_data, "downloads_folder", lambda path: str(tmp_path / path)
    )
    write_eia930a_file(
        tmp_path,
        2023,
        schedule_2=[
            generator_row("BAT1", "Plant One", 1, " gen 1 "),
            generator_row("BAT1", "Plant One", 1, "GEN 2"),
            generator_row("BAT2", "Plant One", 1, "GEN 1"),
            generator_row(
                "BAT1 ", "Battery", "2", " bess1 ", capacity_mw=25, submitted=True
            ),
            generator_row(
                "BAT2",
                "Unknown",
                "Does Not Exist",
                "Does Not Exist (Recommend AS1)",
                submitted=True,
            ),
            generator_row("CISO", "Plant Three", 3, "PV1"),
            generator_row("WALC", "Hoover", 4, "A1"),
        ],
        schedule_3=[
            generator_row("BAT2", "New Solar", "maybe 55086", "S1", submitted=True)
        ],
        schedule_4=[
            schedule_4_row("BAT1", "Plant One", 1, "GEN 1", "Yes", "APS, CAISO"),
            schedule_4_row("CISO", "Plant One", 1, "GEN 2", None, "Dynamic import"),
            schedule_4_row("CISO", "Plant Three", 3, "PV1", None, "Pseudo tie"),
            schedule_4_row("LDWP", "Hoover", 4, "A1", "No", "WAPA"),
        ],
        blank_ba_header=True,
    )
    return tmp_path


def test_load_eia930a_schedule_2(eia930a_folder):
    """Schedule 2 should be loaded without its BA header, using validated IDs where
    available and submitted IDs otherwise, and keeping generators that several BAs
    list."""
    schedule_2 = load_eia930a(2023).schedule_2.set_index(
        ["balancing_authority_code_eia", "plant_name", "generator_id"], drop=False
    )

    assert schedule_2["report_year"].eq(2023).all()
    assert list(schedule_2["balancing_authority_code_eia"]) == [
        "BAT1",
        "BAT1",
        "BAT2",
        "BAT1",
        "BAT2",
        "CISO",
        "WALC",
    ]
    # generator IDs are stripped and upper case
    assert ("BAT1", "Plant One", "GEN 1") in schedule_2.index
    # a generator listed by two BAs is kept for each
    assert len(schedule_2.loc[(slice(None), "Plant One", "GEN 1"), :]) == 2

    battery = schedule_2.loc[("BAT1", "Battery", "BESS1")]
    assert battery["plant_id_eia"] == 2
    assert not battery["validated_eia860"]
    assert battery["capacity_mw"] == 25
    assert battery["county"] == "Asotin"
    assert pd.isna(battery["latitude"])

    validated = schedule_2.loc[("BAT1", "Plant One", "GEN 2")]
    assert validated["validated_eia860"]
    assert validated["plant_id_eia"] == 1
    assert validated["operating_year"] == 2001

    # placeholders are not used as IDs
    placeholder = schedule_2[schedule_2["plant_name"] == "Unknown"].iloc[0]
    assert pd.isna(placeholder["plant_id_eia"])
    assert pd.isna(placeholder["generator_id"])
    assert placeholder["plant_id_as_submitted"] == "Does Not Exist"

    schedule_3 = load_eia930a(2023).schedule_3
    assert pd.isna(schedule_3["plant_id_eia"].iloc[0])
    assert schedule_3["generator_id"].iloc[0] == "S1"


def test_load_eia930a_schedule_4(eia930a_folder):
    """Schedule 4 partners should be mapped to EIA-930 BA codes, and CISO's blank
    "reported in Schedule 2" values filled from its Schedule 2."""
    schedule_4 = load_eia930a(2023).schedule_4.set_index(
        ["balancing_authority_code_eia", "plant_id_eia", "generator_id"]
    )

    host = schedule_4.loc[("BAT1", 1, "GEN 1")]
    assert host["reported_in_schedule_2"]
    assert host["partner_ba_codes"] == "AZPS,CISO"
    assert host["partner_bas_as_submitted"] == "APS, CAISO"
    assert pd.isna(host["arrangement_type"])

    dynamic_import = schedule_4.loc[("CISO", 1, "GEN 2")]
    assert not dynamic_import["reported_in_schedule_2"]
    assert dynamic_import["arrangement_type"] == "dynamic_schedule"
    assert pd.isna(dynamic_import["partner_ba_codes"])

    pseudo_tie = schedule_4.loc[("CISO", 3, "PV1")]
    assert pseudo_tie["reported_in_schedule_2"]
    assert pseudo_tie["arrangement_type"] == "pseudo_tie"

    hoover = schedule_4.loc[("LDWP", 4, "A1")]
    assert not hoover["reported_in_schedule_2"]
    assert pd.isna(hoover["partner_ba_codes"])
    assert hoover["partner_note"] == "WAPA"


def test_parse_partner_bas(caplog):
    ba_codes = {"BPAT", "CHPD", "NEVP", "PJM", "TPWR"}
    assert parse_partner_bas("NVEP; BPA,TWPR", ba_codes) == (
        "BPAT,NEVP,TPWR",
        None,
        None,
    )
    assert parse_partner_bas("CHLN, WAPA", ba_codes) == ("CHPD", None, "WAPA")
    assert parse_partner_bas("Not Generated", ba_codes) == (None, None, "Not Generated")
    assert parse_partner_bas("Pseudo tie", ba_codes) == (None, "pseudo_tie", None)
    assert parse_partner_bas(None, ba_codes) == (None, None, None)
    with caplog.at_level(logging.WARNING):
        assert parse_partner_bas("PJM, XYZ", ba_codes) == ("PJM", None, "XYZ")
    assert "Unrecognized BA 'XYZ'" in caplog.text


def test_load_raw_eia930a_checks_layout(tmp_path, monkeypatch):
    """A file whose layout changed should raise an error."""
    monkeypatch.setattr(
        load_data, "downloads_folder", lambda path: str(tmp_path / path)
    )
    rows = [generator_row("BAT1", "Plant One", 1, "GEN 1")]
    schedule_4 = [schedule_4_row("BAT1", "Plant One", 1, "GEN 1", "Yes", "PJM")]
    write_eia930a_file(tmp_path, 2024, rows, rows, schedule_4)
    workbook = openpyxl.load_workbook(tmp_path / "eia930a/EIA_930A_2024.xlsx")
    workbook["SCH 2"].cell(row=1, column=3).value = "Plant Code"
    workbook.save(tmp_path / "eia930a/EIA_930A_2024.xlsx")

    with pytest.raises(ValueError, match="Unexpected columns"):
        load_raw_eia930a(2024)


def test_clean_eia930a_logs_other_year(eia930a_folder, monkeypatch, caplog):
    """A data year without its own EIA-930A data should use and log another year."""
    monkeypatch.setattr(
        eia930a, "outputs_folder", lambda path: str(eia930a_folder / path)
    )
    monkeypatch.setattr(eia930a, "get_eia930a_year", lambda year: 2023)

    with caplog.at_level(logging.WARNING):
        tables = clean_eia930a(2025, path_prefix="2025/")

    assert "Using 2023 instead" in caplog.text
    assert tables.schedule_2["report_year"].eq(2023).all()
    assert (eia930a_folder / "2025/eia930/eia930a_schedule_4.parquet").exists()
