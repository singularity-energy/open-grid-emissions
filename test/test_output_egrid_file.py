import openpyxl
from openpyxl.utils import get_column_letter
import pandas as pd
import pytest

import oge.output_egrid_file as output_egrid_file
from oge.output_egrid_file import (
    EGRID_CATEGORIES,
    EGRID_SHEETS,
    get_egrid_columns,
    replace_units,
)

# number of columns in each sheet of the eGRID2023 U.S. units and metric files
US_COLUMN_COUNTS = {
    "UNT": 33,
    "GEN": 17,
    "PLNT": 150,
    "ST": 169,
    "BA": 169,
    "SRL": 169,
    "NRL": 169,
    "US": 167,
    "GGL": 6,
}
METRIC_COLUMN_COUNTS = {
    "UNT": 33,
    "GEN": 19,
    "PLNT": 187,
    "ST": 255,
    "BA": 255,
    "SRL": 255,
    "NRL": 255,
    "US": 253,
    "GGL": 9,
}


def headers(columns):
    """Returns the header of each column, keyed by its code."""
    return {code: header for code, header, _ in columns}


@pytest.mark.parametrize("sheet", EGRID_SHEETS)
@pytest.mark.parametrize("metric", [False, True])
def test_column_counts_and_unique_codes(sheet, metric):
    columns = get_egrid_columns(sheet, 2023, metric=metric)
    expected = (METRIC_COLUMN_COUNTS if metric else US_COLUMN_COUNTS)[sheet]
    assert len(columns) == expected
    codes = [code for code, _, _ in columns]
    assert len(set(codes)) == len(codes)


def test_sequence_number_code_uses_data_year():
    assert get_egrid_columns("PLNT", 2025)[0] == (
        "SEQPLT25",
        "Plant file sequence number",
        None,
    )
    assert get_egrid_columns("UNT", 2023)[0][0] == "SEQUNT23"


def test_regional_sheets_share_data_columns():
    us = get_egrid_columns("US", 2023)
    assert us[:3] == [
        ("YEAR", "Data year", None),
        ("USNAMEPCAP", "U.S. nameplate capacity (MW)", None),
        ("USHTIAN", "U.S. annual heat input from combustion (MMBtu)", 1),
    ]
    ba = headers(get_egrid_columns("BA", 2023))
    assert ba["BACODE"] == "Balancing authority code"
    assert ba["BANGENAN"] == "BA annual net generation (MWh)"
    srl = headers(get_egrid_columns("SRL", 2023))
    assert srl["SRCO2RTA"] == (
        "eGRID subregion annual CO2 total output emission rate (lb/MWh)"
    )


def test_metric_columns_convert_units_and_add_gj_columns():
    metric = get_egrid_columns("BA", 2023, metric=True)
    codes = [code for code, _, _ in metric]
    headers_by_code = headers(metric)
    # columns in MWh keep their units and are followed by a column in GJ
    assert codes[codes.index("BANGENAN") + 1] == "BANGENAN2"
    assert headers_by_code["BANGENAN"] == "BA annual net generation (MWh)"
    assert headers_by_code["BANGENAN2"] == "BA annual net generation (GJ)"
    # rates per MWh are converted to kg and followed by a rate per GJ
    assert (
        headers_by_code["BANOXRTA"]
        == "BA annual NOx total output emission rate (kg/MWh)"
    )
    assert (
        headers_by_code["BANOXRTA2"]
        == "BA annual NOx total output emission rate (kg/GJ)"
    )
    # other units are converted without adding a column
    assert headers_by_code["BACO2AN"] == "BA annual CO2 emissions (metric tons)"
    assert headers_by_code["BAHTIAN"] == "BA annual heat input from combustion (GJ)"
    assert headers_by_code["BANOXRA"] == "BA annual NOx input emission rate (kg/GJ)"
    assert headers_by_code["BAHGAN"] == "BA annual Hg emissions (kg)"
    assert "BAHTIAN2" not in headers_by_code
    # units without a metric equivalent are unchanged
    assert headers_by_code["BANAMEPCAP"] == "BA nameplate capacity (MW)"
    assert headers_by_code["BACLPR"] == "BA coal generation percent (resource mix)"


def test_headers_fix_egrid_inconsistencies():
    plant = headers(get_egrid_columns("PLNT", 2023))
    state = headers(get_egrid_columns("ST", 2023))
    # the other unknown/purchased total is labeled renewable in eGRID's plant sheet
    # and nonrenewable in its regional sheets
    assert plant["PLGENATO"] == (
        "Plant annual total nonrenewables other unknown/purchased net generation (MWh)"
    )
    assert state["STGENATO"] == (
        "State annual total nonrenewables other unknown/purchased net generation (MWh)"
    )
    assert state["STNBNOX"] == (
        "State annual NOx nonbaseload output emission rate (lb/MWh)"
    )
    assert plant["CHPCO2E"] == (
        "Plant annual CO2 equivalent emissions CHP adjustment value (tons)"
    )
    assert EGRID_CATEGORIES[14].endswith("per MMBtu)")


def test_replace_units_only_replaces_listed_units():
    units = {"MWh": "GJ"}
    assert replace_units("Direct use (MWh)", units) == "Direct use (GJ)"
    assert replace_units("Total disposition (MWh) without exports", units) == (
        "Total disposition (GJ) without exports"
    )
    assert replace_units("coal percent (resource mix)", units) == (
        "coal percent (resource mix)"
    )
    assert replace_units("nameplate capacity (MW)", units) == (
        "nameplate capacity (MW)"
    )


@pytest.mark.parametrize("sheet", EGRID_SHEETS)
@pytest.mark.parametrize("metric", [False, True])
def test_each_category_is_one_block_of_columns(sheet, metric):
    categories = [category for _, _, category in get_egrid_columns(sheet, 2023, metric)]
    assert set(categories) - {None} <= set(EGRID_CATEGORIES)
    blocks = [
        category
        for i, category in enumerate(categories)
        if category is not None and (i == 0 or categories[i - 1] != category)
    ]
    assert len(blocks) == len(set(blocks))


def test_categories_by_sheet():
    def category_counts(sheet):
        counts = {}
        for _, _, category in get_egrid_columns(sheet, 2023):
            counts[category] = counts.get(category, 0) + 1
        return counts

    assert category_counts("UNT") == {None: 19, 2: 14}
    assert category_counts("GEN") == {None: 15, 1: 2}
    assert category_counts("PLNT") == {
        None: 37,
        1: 15,
        2: 22,
        3: 16,
        4: 8,
        5: 8,
        6: 8,
        7: 11,
        8: 4,
        9: 3,
        10: 11,
        11: 4,
        12: 3,
    }
    assert category_counts("BA") == {
        None: 4,
        1: 15,
        4: 8,
        5: 8,
        6: 8,
        7: 11,
        8: 4,
        9: 3,
        10: 11,
        11: 4,
        12: 3,
        13: 30,
        14: 30,
        15: 8,
        16: 11,
        17: 11,
    }
    assert category_counts("GGL") == {None: 6}


def test_gj_columns_have_the_category_of_their_column():
    metric = get_egrid_columns("PLNT", 2023, metric=True)
    categories = {code: category for code, _, category in metric}
    assert categories["PLNGENAN"] == categories["PLNGENAN2"] == 1
    assert categories["PLGENACL"] == categories["PLGENACL2"] == 7


@pytest.fixture
def results(tmp_path, monkeypatch):
    """Writes small 2024 annual results and intermediate outputs for two plants in two
    BAs to temporary results and outputs folders.

    Plant 1 (BA1) burns gas and plant 2 (BA2) is solar. BA2 has no coal or other
    generation. Plant 1 has two CEMS units: unit 1 reports Hg in one of its two hours
    and unit 2 does not report Hg. Plant 1 has one generator (GT1) and plant 2 has one
    (PV1), which retired in 2024.
    """
    monkeypatch.setattr(
        output_egrid_file, "results_folder", lambda rel="": str(tmp_path / rel)
    )
    monkeypatch.setattr(
        output_egrid_file,
        "outputs_folder",
        lambda rel="": str(tmp_path / "outputs" / rel),
    )
    plant_folder = tmp_path / "2024/plant_data"
    (plant_folder / "annual/us_units").mkdir(parents=True)
    plant_data = pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "plant_name_eia": ["Gas Plant", "Solar Plant"],
            "ba_code": ["BA1", "BA2"],
            "state": ["TX", "AZ"],
            "county": ["Travis", "Pima"],
            "capacity_mw": [100.0, 50.0],
            "net_generation_mwh": [1000.0, 500.0],
            "fuel_consumed_mmbtu": [12000.0, 0.0],
            "fuel_consumed_for_electricity_mmbtu": [10000.0, 0.0],
            "co2_mass_lb": [2_000_000.0, 0.0],
            "co2_mass_lb_for_electricity_adjusted": [1_000_000.0, 0.0],
            "generated_co2_rate_lb_per_mwh_for_electricity_adjusted": [1000.0, 0.0],
        }
    )
    for column in output_egrid_file.PLANT_DATA_SOURCES.values():
        if column not in plant_data and column not in [
            "latitude",
            "longitude",
            "ba_name",
        ]:
            plant_data[column] = 0.0
    plant_data.to_csv(plant_folder / "annual/us_units/plant_data.csv", index=False)
    plant_data[["plant_id_eia", "plant_name_eia", "state"]].assign(
        latitude=[30.0, 32.0], longitude=[-97.0, -111.0]
    ).to_csv(plant_folder / "plant_static_attributes.csv", index=False)
    outputs_folder = tmp_path / "outputs/2024"
    outputs_folder.mkdir(parents=True)
    pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 1],
            "emissions_unit_id_epa": ["1", "1", "2", "2"],
            "energy_source_code": ["NG"] * 4,
            "operating_time_hours": [1.0, 0.5, 1.0, 0.0],
            "fuel_consumed_mmbtu": [100.0, 50.0, 80.0, None],
            "nox_mass_lb": [10.0, 10.0, 4.0, None],
            "so2_mass_lb": [None, None, 2.0, None],
            "co2_mass_lb": [4000.0, 2000.0, 1000.0, None],
            "hg_mass_lb": [0.25, None, None, None],
        }
    ).to_csv(outputs_folder / "cems_cleaned_2024.csv.zip", index=False)
    pd.DataFrame(
        {
            "report_date": ["2024-01-01", "2024-02-01", "2024-01-01"],
            "plant_id_eia": [1, 1, 2],
            "generator_id": ["GT1", "GT1", "PV1"],
            "prime_mover_code": ["GT", "GT", "PV"],
            "net_generation_mwh": [600.0, 400.0, 500.0],
        }
    ).to_csv(outputs_folder / "eia923_allocated_2024.csv.zip", index=False)
    pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 2],
            "emissions_unit_id_epa": ["1", "2", None],
            "generator_id": ["GT1", "GT1", "PV1"],
            "generator_operating_date": ["2001-06-01", "2001-06-01", "2015-03-01"],
            "generator_retirement_date": [None, None, "2024-09-30"],
        }
    ).to_csv(outputs_folder / "subplant_crosswalk_2024.csv.zip", index=False)
    fleet_folder = tmp_path / "2024/power_sector_data/annual/us_units"
    fleet_folder.mkdir(parents=True)
    columns = {column for _, column in output_egrid_file.REGION_DATA_SOURCES.values()}
    for ba, rows in {
        "BA1": {"total": 1000.0, "natural_gas": 900.0, "coal": 50.0, "other": 50.0},
        "BA2": {"total": 500.0, "solar": 500.0},
        "US": {"total": 1500.0, "natural_gas": 900.0, "coal": 50.0, "solar": 500.0},
    }.items():
        fleet = pd.DataFrame({"fuel_category": list(rows)})
        for column in columns:
            fleet[column] = 0.0
        fleet["net_generation_mwh"] = list(rows.values())
        fleet["co2_mass_lb_for_electricity_adjusted"] = [
            1_000_000.0 if fuel in ["total", "natural_gas"] else 0.0 for fuel in rows
        ]
        fleet.to_csv(fleet_folder / f"{ba}.csv", index=False)
    return tmp_path


def test_plant_table(results):
    plants = output_egrid_file.build_plant_table(2024, "2024/")
    # plants are sorted by state and numbered in that order
    assert list(plants["PNAME"]) == ["Solar Plant", "Gas Plant"]
    assert list(plants["SEQPLT24"]) == [1, 2]
    assert list(plants["YEAR"]) == [2024, 2024]
    gas = plants.set_index("ORISPL").loc[1]
    # masses in tons are converted from lb, and masses in lbs are not
    assert gas["PLCO2AN"] == 500.0
    assert gas["UNCO2"] == 1000.0
    assert gas["PLCO2RTA"] == 1000.0
    assert gas["PLHTIANT"] == 10000.0
    assert gas["UNHTIT"] == 12000.0
    assert gas["LAT"] == 30.0


def test_unit_table(results):
    units = output_egrid_file.build_unit_table(2024, "2024/").set_index("UNITID")
    assert list(units.index) == ["1", "2"]
    assert list(units["SEQUNT24"]) == [1, 2]
    assert list(units["CAPDFLAG"]) == ["Yes", "Yes"]
    assert list(units["PNAME"]) == ["Gas Plant", "Gas Plant"]
    assert units.loc["1", "ORISPL"] == 1
    assert units.loc["1", "FUELU1"] == "NG"
    # operating hours count the hours with any operation
    assert units.loc["1", "HRSOP"] == 2
    assert units.loc["2", "HRSOP"] == 1
    assert units.loc["1", "HTIAN"] == 150.0
    # masses are converted from lb to tons
    assert units.loc["1", "CO2AN"] == 3.0
    assert units.loc["1", "NOXAN"] == 0.01
    # Hg is the sum of the hours with Hg, and is blank for units without Hg
    assert units.loc["1", "HGAN"] == 0.25
    assert pd.isna(units.loc["2", "HGAN"])
    # pollutants that a unit does not report are blank
    assert pd.isna(units.loc["1", "SO2AN"])
    assert units.loc["2", "SO2AN"] == 0.001


def test_generator_table(results):
    generators = output_egrid_file.build_generator_table(2024, "2024/")
    # generators are sorted by state, so plant 2 (AZ) comes first
    assert list(generators["GENID"]) == ["PV1", "GT1"]
    assert list(generators["SEQGEN24"]) == [1, 2]
    generators = generators.set_index("GENID")
    assert generators.loc["GT1", "ORISPL"] == 1
    assert generators.loc["GT1", "PNAME"] == "Gas Plant"
    assert generators.loc["GT1", "PSTATABB"] == "TX"
    assert generators.loc["GT1", "PRMVR"] == "GT"
    assert generators.loc["GT1", "GENNTAN"] == 1000.0
    assert generators.loc["GT1", "GENYRONL"] == 2001
    assert pd.isna(generators.loc["GT1", "GENYRRET"])
    assert generators.loc["PV1", "GENYRRET"] == 2024


def test_region_tables(results):
    ba = output_egrid_file.build_region_table(2024, "2024/", "BA", ["BA1", "BA2"])
    ba = ba.set_index("BACODE")
    assert ba.loc["BA1", "BANGENAN"] == 1000.0
    assert ba.loc["BA1", "BACO2AN"] == 500.0
    assert ba.loc["BA1", "BAGENACL"] == 50.0
    # other generation is reported as other, nonrenewable and noncombustion
    assert ba.loc["BA1", ["BAGENAOP", "BAGENATO", "BAGENACO"]].tolist() == [50.0] * 3
    # BA2 has no coal or other generation
    assert ba.loc["BA2", "BAGENACL"] == 0
    assert ba.loc["BA2", "BAGENAOP"] == 0
    assert ba.loc["BA2", "BAGENASO"] == 500.0
    us = output_egrid_file.build_region_table(2024, "2024/", "US", ["US"])
    assert us.loc[0, "USNGENAN"] == 1500.0
    assert "BACODE" not in us


def test_egrid_tables_have_every_column(results):
    tables = output_egrid_file.build_egrid_tables(2024, "2024/")
    assert list(tables) == EGRID_SHEETS
    for sheet, table in tables.items():
        assert list(table.columns) == [
            code for code, _, _ in get_egrid_columns(sheet, 2024)
        ]
    assert list(tables["BA"]["BACODE"]) == ["BA1", "BA2"]
    assert len(tables["UNT"]) == 2
    assert len(tables["GEN"]) == 2
    assert tables["ST"].empty
    # columns that OGE does not calculate are blank
    assert tables["PLNT"]["PLGENACL"].isna().all()


def test_convert_to_metric():
    table = pd.DataFrame(
        {code: [1.0] for code, _, _ in get_egrid_columns("US", 2024)}
    ).assign(YEAR=2024)
    metric = output_egrid_file.convert_to_metric(table, "US", 2024)
    assert list(metric.columns) == [
        code for code, _, _ in get_egrid_columns("US", 2024, metric=True)
    ]
    assert metric.loc[0, "YEAR"] == 2024
    assert metric.loc[0, "USNGENAN"] == 1.0
    assert metric.loc[0, "USNGENAN2"] == pytest.approx(3.6)
    assert metric.loc[0, "USCO2AN"] == pytest.approx(0.907184)
    assert metric.loc[0, "USCH4AN"] == pytest.approx(0.453592)
    assert metric.loc[0, "USHTIANT"] == pytest.approx(1.055056)
    assert metric.loc[0, "USCO2RTA"] == pytest.approx(0.453592)
    assert metric.loc[0, "USCO2RTA2"] == pytest.approx(0.453592 / 3.6)
    assert metric.loc[0, "USCO2RA"] == pytest.approx(0.453592 / 1.055056)
    assert metric.loc[0, "USCLPR"] == 1.0


def test_write_egrid_files(results):
    output_egrid_file.write_egrid_files(2024, "2024/")
    for file, metric in [
        ("oge_egrid2024_data.xlsx", False),
        ("oge_egrid2024_data_metric.xlsx", True),
    ]:
        path = results / "2024/eGRID" / file
        workbook = openpyxl.load_workbook(path)
        assert workbook.sheetnames == ["Contents"] + [
            f"{sheet}24" for sheet in EGRID_SHEETS
        ]
        worksheet = workbook["PLNT24"]
        columns = get_egrid_columns("PLNT", 2024, metric)
        assert [cell.value for cell in worksheet[1]] == [h for _, h, _ in columns]
        assert [cell.value for cell in worksheet[2]] == [c for c, _, _ in columns]
        assert worksheet.max_row == 4
        codes = [code for code, _, _ in columns]
        # filled columns are colored by category and blank columns are grayed out
        filled = worksheet.cell(1, codes.index("PLNGENAN") + 1)
        blank = worksheet.cell(1, codes.index("PLGENACL") + 1)
        assert filled.fill.fgColor.rgb == "FFF2DCDB"
        assert blank.font.color.rgb == "FFA6A6A6"


def test_skip_outputs_writes_nothing(results):
    output_egrid_file.write_egrid_files(2024, "2024/", skip_outputs=True)
    assert not (results / "2024/eGRID").exists()


def test_number_formats_follow_egrid():
    plant = output_egrid_file.get_number_formats("PLNT", 2024)
    assert plant["PNAME"] == "General"
    assert plant["NAMEPCAP"] == "#,##0.0"
    assert plant["PLNGENAN"] == output_egrid_file.QUANTITY_FORMAT
    assert plant["PLCO2RTA"] == output_egrid_file.RATE_FORMAT
    assert plant["PLHGAN"] == output_egrid_file.HG_MASS_FORMAT
    assert plant["UNCO2"] == "#,##0"
    assert plant["PLCLPR"] == output_egrid_file.PERCENT_FORMAT
    state = output_egrid_file.get_number_formats("ST", 2024)
    assert state["STNAMEPCAP"] == output_egrid_file.QUANTITY_FORMAT
    assert state["STCCH4RT"] == output_egrid_file.SMALL_RATE_FORMAT
    assert state["STCCO2RT"] == output_egrid_file.RATE_FORMAT
    unit = output_egrid_file.get_number_formats("UNT", 2024)
    assert unit["NOXAN"] == output_egrid_file.UNIT_MASS_FORMAT
    assert unit["HTIAN"] == "#,##0"
    # GJ columns have the format of the column they follow
    metric = output_egrid_file.get_number_formats("US", 2024, metric=True)
    assert metric["USNGENAN2"] == metric["USNGENAN"]
    assert metric["USCO2RTA2"] == output_egrid_file.RATE_FORMAT
    assert set(metric) == {code for code, _, _ in get_egrid_columns("US", 2024, True)}


def test_contents_sheet_links(results):
    output_egrid_file.write_egrid_files(2024, "2024/")
    worksheet = openpyxl.load_workbook(results / "2024/eGRID/oge_egrid2024_data.xlsx")[
        "Contents"
    ]
    values = [cell.value for row in worksheet.iter_rows() for cell in row]
    assert "eGRID2024_Data.xlsx" in values
    assert "Balancing authority area year 2024 data" in values
    assert "Color Coding Legend" in values
    links = {
        cell.coordinate: cell.hyperlink.location
        for row in worksheet.iter_rows()
        for cell in row
        if cell.hyperlink is not None and cell.hyperlink.location
    }
    assert "'PLNT24'!A1" in links.values()
    # the legend links to the first column of each category on each sheet, e.g. the
    # plant annual heat input
    plant_heat_input = [code for code, _, _ in get_egrid_columns("PLNT", 2024)].index(
        "PLHTIAN"
    )
    assert f"'PLNT24'!{get_column_letter(plant_heat_input + 1)}1" in links.values()


def test_describe_sheets():
    assert output_egrid_file.describe_sheets() == (
        "Unit, Generator, Plant, State, Balancing Authority Area, eGRID Subregion, "
        "NERC Region, U.S., and Grid Gross Loss (%)"
    )


def test_sheets_filter_on_column_codes(results):
    output_egrid_file.write_egrid_files(2024, "2024/")
    workbook = openpyxl.load_workbook(results / "2024/eGRID/oge_egrid2024_data.xlsx")
    columns = get_egrid_columns("PLNT", 2024)
    assert workbook["PLNT24"].auto_filter.ref == (
        f"A2:{get_column_letter(len(columns))}4"
    )
    # sheets with no rows filter on the codes alone
    assert workbook["ST24"].auto_filter.ref == "A2:FM2"
