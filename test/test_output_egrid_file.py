from pathlib import Path

import numpy as np
import openpyxl
from openpyxl.utils import get_column_letter
import pandas as pd
import pytest

import oge.output_egrid_file as output_egrid_file
from oge.column_checks import DATA_COLUMNS, HG_DATA_COLUMNS
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
    (PV1), which retired in 2024. The CAMPD facility file and EIA-860 tables are
    replaced by small tables: plant 2 is missing its county from OGE's attributes, so it
    comes from EIA-860.
    """
    monkeypatch.setattr(
        output_egrid_file, "results_folder", lambda rel="": str(tmp_path / rel)
    )
    monkeypatch.setattr(
        output_egrid_file,
        "outputs_folder",
        lambda rel="": str(tmp_path / "outputs" / rel),
    )
    monkeypatch.setattr(
        output_egrid_file,
        "load_campd_facility",
        lambda year: pd.DataFrame(
            {
                "plant_id_eia": [1, 1],
                "emissions_unit_id_epa": ["1", "2"],
                "operating_status_epa": ["Operating", "Operating (Retired 06/30/2024)"],
                "program_code_epa": ["ARP, CSOSG2", "ARP"],
                "unit_type_epa": ["Combined cycle", "Combustion turbine"],
                "so2_controls_epa": [pd.NA, pd.NA],
                "nox_controls_epa": [
                    "Water Injection (Retired Jun 01, 2024)|"
                    "Selective Catalytic Reduction",
                    "Low NOx Burner Technology (Dry Bottom only)",
                ],
                "hg_controls_epa": [pd.NA, "Halogenated PAC Sorbent Injection"],
                "commercial_operation_date_epa": pd.to_datetime(
                    ["2001-06-01", "1975-01-01"]
                ),
            }
        ),
    )
    pudl_tables = {
        "core_eia860__scd_generators": pd.DataFrame(
            {
                "plant_id_eia": [1, 2],
                "generator_id": ["GT1", "PV1"],
                "operational_status_code": ["op", "RE"],
                "energy_source_code_1": ["NG", "SUN"],
                "capacity_mw": [100.0, 50.0],
            }
        ),
        "core_eia__entity_plants": pd.DataFrame(
            {
                "plant_id_eia": [1, 2],
                "plant_name_eia": ["Gas Plant EIA", "Solar Plant EIA"],
                "state": ["TX", "AZ"],
                "county": ["Travis EIA", "Pima"],
                "latitude": [30.5, 32.5],
                "longitude": [-97.5, -111.5],
                "state_id_fips": ["48", "04"],
                "county_id_fips": ["48453", "04019"],
            }
        ),
        "core_eia860__scd_plants": pd.DataFrame(
            {
                "plant_id_eia": [1, 2],
                "transmission_distribution_owner_name": ["Wires Co", "Wires Co"],
                "transmission_distribution_owner_id": [10, 10],
                "utility_id_eia": [20, 21],
                "sector_id_eia": [2, 1],
                "nerc_region": ["TRE", "WECC"],
            }
        ),
        "core_eia__entity_utilities": pd.DataFrame(
            {"utility_id_eia": [20, 21], "utility_name_eia": ["Gas Co", "Solar Co"]}
        ),
    }
    monkeypatch.setattr(
        output_egrid_file,
        "load_pudl_table",
        lambda table, year=None, columns=None: pudl_tables[table],
    )
    plant_folder = tmp_path / "2024/plant_data"
    (plant_folder / "annual/us_units").mkdir(parents=True)
    plant_data = pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "plant_name_eia": ["Gas Plant", "Solar Plant"],
            "ba_code": ["BA1", "BA2"],
            "state": ["TX", "AZ"],
            "county": ["Travis", None],
            "capacity_mw": [100.0, 50.0],
            "net_generation_mwh": [1000.0, 500.0],
            "fuel_consumed_mmbtu": [12000.0, 0.0],
            "fuel_consumed_for_electricity_mmbtu": [10000.0, 0.0],
            "co2_mass_lb": [2_000_000.0, 0.0],
            "co2_mass_lb_adjusted": [1_500_000.0, 0.0],
            "co2_mass_lb_for_electricity_adjusted": [1_000_000.0, 0.0],
            "co2e_mass_lb_for_electricity_adjusted": [1_100_000.0, 0.0],
            "generated_co2_rate_lb_per_mwh_for_electricity_adjusted": [1000.0, 0.0],
        }
    )
    # OGE's other data columns are zero, and Hg is incomplete for the year
    for column in DATA_COLUMNS:
        if column not in plant_data:
            plant_data[column] = 0.0
    for column in HG_DATA_COLUMNS:
        plant_data[column] = np.nan
    plant_data.to_csv(plant_folder / "annual/us_units/plant_data.csv", index=False)
    plant_data[["plant_id_eia", "plant_name_eia", "state"]].assign(
        latitude=[30.0, 32.0],
        longitude=[-97.0, -111.0],
        plant_primary_fuel=["NG", "SUN"],
    ).to_csv(plant_folder / "plant_static_attributes.csv", index=False)
    # monthly subplant data: plant 1 has a gas subplant that runs in January and June
    # and a solar subplant that runs in June, and plant 2 is solar
    (plant_folder / "monthly/us_units").mkdir(parents=True)
    pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 2],
            "subplant_id": [1, 1, 2, 1],
            "report_date": ["2024-01-01", "2024-06-01", "2024-06-01", "2024-06-01"],
            "fuel_consumed_mmbtu": [7000.0, 5000.0, 100.0, 1000.0],
            "fuel_consumed_for_electricity_mmbtu": [6000.0, 4000.0, 100.0, 1000.0],
            "net_generation_mwh": [600.0, 400.0, 30.0, 500.0],
            "nox_mass_lb": [20.0, 10.0, 0.0, 0.0],
            "nox_mass_lb_adjusted": [20.0, 10.0, 0.0, 0.0],
            "nox_mass_lb_for_electricity_adjusted": [16.0, 8.0, 0.0, 0.0],
            "so2_mass_lb": [0.0] * 4,
            "co2_mass_lb": [8000.0, 6000.0, 0.0, 0.0],
            "ch4_mass_lb": [1.0, 1.0, 0.0, 0.0],
            "n2o_mass_lb": [0.1, 0.1, 0.0, 0.0],
            "co2e_mass_lb": [8100.0, 6100.0, 0.0, 0.0],
            "hg_mass_lb": [0.5, None, None, None],
        }
    ).to_csv(plant_folder / "monthly/us_units/subplant_data.csv", index=False)
    # the gas subplant's January data is from CEMS and its June data from EIA
    pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 2],
            "subplant_id": [1, 1, 2, 1],
            "report_date": ["2024-01-01", "2024-06-01", "2024-06-01", "2024-06-01"],
            "data_source": ["CEMS", "EIA", "EIA", "EIA"],
        }
    ).to_csv(plant_folder / "plant_metadata.csv", index=False)
    pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 2],
            "subplant_id": [1, 2, 1],
            "fuel_category": ["natural_gas", "solar", "solar"],
        }
    ).to_csv(plant_folder / "subplant_attributes.csv", index=False)
    outputs_folder = tmp_path / "outputs/2024"
    outputs_folder.mkdir(parents=True)
    pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 1, 1],
            "emissions_unit_id_epa": ["1", "1", "2", "2"],
            "report_date": ["2024-01-01", "2024-06-01", "2024-06-01", "2024-10-01"],
            "energy_source_code": ["NG"] * 4,
            "operating_time_hours": [1.0, 0.5, 1.0, 0.0],
            "fuel_consumed_mmbtu": [100.0, 50.0, 80.0, None],
            "nox_mass_lb": [10.0, 10.0, 4.0, None],
            "so2_mass_lb": [None, None, 2.0, None],
            "co2_mass_lb": [4000.0, 2000.0, 1000.0, None],
            "hg_mass_lb": [0.25, None, None, None],
            "nox_mass_measurement_code": ["Measured", "Substitute", "Measured", None],
            "so2_mass_measurement_code": ["Measured"] * 4,
            "co2_mass_measurement_code": ["Measured", "Measured", "Imputed", None],
            "hg_mass_measurement_code": ["Measured", "Unavailable", None, None],
        }
    ).to_csv(outputs_folder / "cems_cleaned_2024.csv.zip", index=False)
    pd.DataFrame(
        {
            "report_date": ["2024-01-01", "2024-06-01", "2024-01-01"],
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
    # power sector results, with half of each annual value in January and half in June.
    # Heat input is ten times net generation
    fleet_folder = tmp_path / "2024/power_sector_data/{resolution}/us_units"
    for resolution in ["annual", "monthly"]:
        Path(str(fleet_folder).format(resolution=resolution)).mkdir(parents=True)
    columns = {
        column for _, column, _ in output_egrid_file.REGION_DATA_SOURCES.values()
    }
    for ba, rows in {
        "BA1": {"total": 1000.0, "natural_gas": 900.0, "coal": 50.0, "other": 50.0},
        "BA2": {"total": 500.0, "solar": 500.0},
        "US": {"total": 1500.0, "natural_gas": 900.0, "coal": 50.0, "solar": 500.0},
    }.items():
        fleet = pd.DataFrame({"fuel_category": list(rows)})
        for column in columns:
            fleet[column] = 0.0
        fleet["net_generation_mwh"] = list(rows.values())
        fleet["fuel_consumed_for_electricity_mmbtu"] = fleet["net_generation_mwh"] * 10
        fleet["co2_mass_lb_for_electricity_adjusted"] = [
            1_000_000.0 if fuel in ["total", "natural_gas"] else 0.0 for fuel in rows
        ]
        fleet.to_csv(
            Path(str(fleet_folder).format(resolution="annual")) / f"{ba}.csv",
            index=False,
        )
        monthly = pd.concat(
            [fleet.assign(report_date=date) for date in ["2024-01-01", "2024-06-01"]]
        )
        monthly[list(columns)] = monthly[list(columns)] / 2
        monthly.to_csv(
            Path(str(fleet_folder).format(resolution="monthly")) / f"{ba}.csv",
            index=False,
        )
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
    assert gas["PLPRMFL"] == "NG"
    assert gas["PLCO2EQA"] == 550.0


def test_plant_ozone_season_and_combustion_heat_input(results):
    plants = output_egrid_file.build_plant_table(2024, "2024/").set_index("ORISPL")
    # combustion heat input leaves out plant 1's solar subplant
    assert plants.loc[1, "PLHTIAN"] == 10000.0
    assert plants.loc[1, "UNHTI"] == 12000.0
    # the ozone season is May to September
    assert plants.loc[1, "PLHTIOZ"] == 4000.0
    assert plants.loc[1, "PLHTIOZT"] == 4100.0
    assert plants.loc[1, "UNHTIOZ"] == 5000.0
    assert plants.loc[1, "UNHTIOZT"] == 5100.0
    assert plants.loc[1, "PLNGENOZ"] == 430.0
    assert plants.loc[1, "PLNOXOZ"] == 0.004
    assert plants.loc[1, "UNNOXOZ"] == 0.005
    # plants without combustion subplants have no combustion heat input
    assert plants.loc[2, "PLHTIAN"] == 0
    assert plants.loc[2, "PLHTIOZT"] == 1000.0


def test_plant_generation_by_fuel(results):
    plants = output_egrid_file.build_plant_table(2024, "2024/").set_index("ORISPL")
    # each subplant's generation counts toward the groups of its OGE fuel category
    assert plants.loc[1, "PLGENAGS"] == 1000.0
    assert plants.loc[1, "PLGENASO"] == 30.0
    assert plants.loc[1, "PLGENACL"] == 0
    assert plants.loc[1, "PLGENATN"] == 1000.0
    assert plants.loc[1, "PLGENATR"] == 30.0
    assert plants.loc[1, "PLGENATH"] == 30.0
    assert plants.loc[1, "PLGENACY"] == 1000.0
    assert plants.loc[1, "PLGENACN"] == 30.0
    assert plants.loc[2, "PLGENASO"] == 500.0
    assert plants.loc[2, "PLGENACY"] == 0


def test_fuel_groups_use_oge_fuel_categories():
    groups = output_egrid_file.EGRID_FUEL_GROUPS
    assert groups["gas"] == ["natural_gas"]
    assert groups["other_fossil"] == ["waste"]
    assert groups["other"] == ["other", "storage"]
    # every OGE fuel category is in exactly one of eGRID's fuels
    fuels = ["coal", "oil", "gas", "nuclear", "hydro", "biomass", "wind", "solar"]
    fuels += ["geothermal", "other_fossil", "other"]
    categories = [category for fuel in fuels for category in groups[fuel]]
    assert len(categories) == len(set(categories)) == 12
    # nonrenewable and renewable, and combustion and noncombustion, together cover
    # every category except other unknown/purchased fuel
    assert sorted(groups["nonrenewable"] + groups["renewable"] + groups["other"]) == (
        sorted(categories)
    )
    assert sorted(
        groups["combustion"] + groups["noncombustion"] + groups["other"]
    ) == sorted(categories)


def test_plant_attributes_use_oge_then_eia860(results):
    plants = output_egrid_file.build_plant_table(2024, "2024/").set_index("ORISPL")
    # OGE's attributes are used where OGE has them, and EIA-860's otherwise
    assert plants.loc[1, "PNAME"] == "Gas Plant"
    assert plants.loc[1, "CNTYNAME"] == "Travis"
    assert plants.loc[2, "CNTYNAME"] == "Pima"
    assert plants.loc[1, "OPRNAME"] == "Wires Co"
    assert plants.loc[1, "UTLSRVID"] == 20
    assert plants.loc[1, "UTLSRVNM"] == "Gas Co"
    assert plants.loc[1, "SECTOR"] == "IPP Non-CHP"
    assert plants.loc[2, "SECTOR"] == "Electric Utility"
    assert plants.loc[1, "NERC"] == "TRE"
    assert plants.loc[1, "FIPSST"] == "48"
    assert plants.loc[1, "FIPSCNTY"] == "453"


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
    # the ozone season is May to September
    assert units.loc["1", "HTIOZ"] == 50.0
    assert units.loc["1", "NOXOZ"] == 0.005
    assert units.loc["2", "HTIOZ"] == 80.0
    # masses are converted from lb to tons
    assert units.loc["1", "CO2AN"] == 3.0
    assert units.loc["1", "NOXAN"] == 0.01
    # Hg is the sum of the hours with Hg, and is blank for units without Hg
    assert units.loc["1", "HGAN"] == 0.25
    assert pd.isna(units.loc["2", "HGAN"])
    # pollutants that a unit does not report are blank
    assert pd.isna(units.loc["1", "SO2AN"])
    assert units.loc["2", "SO2AN"] == 0.001
    # unit attributes come from the CAMPD facility file, as eGRID codes
    assert units.loc["1", "PRGCODE"] == "ARP, CSOSG2"
    assert units.loc["1", "UNTOPST"] == "OP"
    assert units.loc["2", "UNTOPST"] == "OP"
    assert units.loc["1", "BOTFIRTY"] == "CC"
    # the prime mover follows from the CAMPD unit type, as in eGRID
    assert units.loc["1", "PRMVR"] == "CT"
    assert units.loc["2", "PRMVR"] == "GT"
    assert units.loc["1", "NOXCTLDV"] == "H2O,SCR"
    assert units.loc["2", "NOXCTLDV"] == "LNB"
    assert pd.isna(units.loc["1", "SO2CTLDV"])
    assert pd.isna(units.loc["1", "HGCTLDV"])
    assert units.loc["2", "HGCTLDV"] == "Yes"
    assert units.loc["1", "UNTYRONL"] == 2001
    # sources are the measurement codes, from the one covering the most emissions
    assert units.loc["1", "CO2SRC"] == "Measured"
    assert units.loc["2", "CO2SRC"] == "Imputed"
    assert units.loc["1", "NOXANSRC"] == "Measured; Substitute"
    assert units.loc["1", "NOXOZSRC"] == "Substitute"
    assert units.loc["1", "HGSRC"] == "Measured"
    assert pd.isna(units.loc["2", "HGSRC"])
    assert units.loc["1", "HTIANSRC"] == "EPA/CAPD"
    assert units.loc["2", "SO2SRC"] == "Measured"


def test_translate_campd_values():
    values = pd.Series(
        ["Wet Lime FGD|Wet Limestone (Began Jan 01, 2024)", None, "New control"]
    )
    translated = output_egrid_file.translate_campd_values(
        values, output_egrid_file.CAMPD_SO2_CONTROL_CODES
    )
    assert translated[0] == "WL,WLS"
    assert pd.isna(translated[1])
    # values without a code are kept as they are
    assert translated[2] == "New control"


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
    assert generators.loc["GT1", "GENNTOZ"] == 400.0
    assert generators.loc["PV1", "GENNTOZ"] == 0
    assert generators.loc["GT1", "GENYRONL"] == 2001
    assert pd.isna(generators.loc["GT1", "GENYRRET"])
    assert generators.loc["PV1", "GENYRRET"] == 2024
    # generator attributes come from EIA-860
    assert generators.loc["GT1", "GENSTAT"] == "OP"
    assert generators.loc["PV1", "GENSTAT"] == "RE"
    assert generators.loc["GT1", "FUELG1"] == "NG"
    assert generators.loc["PV1", "NAMEPCAP"] == 50.0


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
    # combustion heat input sums the combustion fuel categories, and BA2 has none
    assert ba.loc["BA1", "BAHTIAN"] == 9500.0
    assert ba.loc["BA1", "BAHTIANT"] == 10000.0
    assert ba.loc["BA2", "BAHTIAN"] == 0
    # the ozone season sums May to September of the monthly results
    assert ba.loc["BA1", "BANGENOZ"] == 500.0
    assert ba.loc["BA1", "BAHTIOZT"] == 5000.0
    assert ba.loc["BA1", "BAHTIOZ"] == 4750.0
    # generation by fuel uses OGE's fuel categories
    assert ba.loc["BA1", "BAGENAGS"] == 900.0
    assert ba.loc["BA1", "BAGENATN"] == 950.0
    assert ba.loc["BA1", "BAGENACY"] == 950.0
    assert ba.loc["BA1", "BAGENATR"] == 0
    assert ba.loc["BA2", "BAGENATR"] == 500.0
    assert ba.loc["BA2", "BAGENACN"] == 500.0
    # rates by fuel type divide each fuel group's emissions by its generation or heat
    # input, and the fossil group includes coal
    assert ba.loc["BA1", "BAGCO2RT"] == pytest.approx(1_000_000 / 900)
    assert ba.loc["BA1", "BAGCO2R"] == pytest.approx(1_000_000 / 9000)
    assert ba.loc["BA1", "BAFSC2RT"] == pytest.approx(1_000_000 / 950)
    assert ba.loc["BA1", "BACCO2RT"] == 0
    assert ba.loc["BA1", "BAGNXORT"] == 0
    # rates are blank where a BA has none of the fuel group
    assert pd.isna(ba.loc["BA1", "BAOCO2RT"])
    assert pd.isna(ba.loc["BA2", "BAFSC2RT"])
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
    assert tables["PLNT"]["CAPFAC"].isna().all()


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
        blank = worksheet.cell(1, codes.index("CAPFAC") + 1)
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


def test_emission_rates_and_resource_mix():
    columns = [code for code, _, _ in get_egrid_columns("PLNT", 2024)]
    table = pd.DataFrame({code: [0.0, 0.0] for code in columns}).assign(
        PLNOXAN=[1.0, 1.0],
        PLNOXOZ=[0.5, 0.5],
        PLCH4AN=[30.0, 30.0],
        PLCO2EQA=[600.0, 0.0],
        PLNGENAN=[1000.0, 0.0],
        PLNGENOZ=[400.0, 0.0],
        PLHTIAN=[10000.0, 500.0],
        PLHTIOZ=[4000.0, 200.0],
        PLGENACY=[800.0, -10.0],
        PLGENAGS=[800.0, -10.0],
        PLGENASO=[200.0, 0.0],
        PLGENAOP=[-50.0, 0.0],
        PLGENATN=[800.0, -10.0],
        PLGENATR=[200.0, 0.0],
    )
    table[f"PL{output_egrid_file.OZONE_COMBUSTION_GENERATION}"] = [320.0, 0.0]
    rates = output_egrid_file.add_emission_rates_and_resource_mix(table, "PLNT", 2024)
    gas, idle = rates.iloc[0], rates.iloc[1]
    # tons are converted to lb, and masses in lb are used as they are
    assert gas["PLNOXRTA"] == 2.0
    assert gas["PLCH4RTA"] == 0.03
    assert gas["PLC2ERTA"] == 1200.0
    assert gas["PLNOXRA"] == 0.2
    assert gas["PLNOXCRT"] == 2.5
    # ozone season rates divide ozone season NOx by ozone season denominators
    assert gas["PLNOXRTO"] == 2.5
    assert gas["PLNOXRO"] == 0.25
    assert gas["PLNOXCRO"] == 3.125
    # rates are blank where the denominator is zero or negative
    assert pd.isna(idle["PLNOXRTA"])
    assert pd.isna(idle["PLNOXCRT"])
    assert idle["PLNOXRA"] == 4.0
    # shares are of positive generation, and negative generation has a 0% share
    assert gas["PLGSPR"] == 0.8
    assert gas["PLSOPR"] == 0.2
    assert gas["PLOPPR"] == 0
    assert gas["PLTNPR"] == 0.8
    assert gas["PLCYPR"] == 0.8
    assert pd.isna(idle["PLGSPR"])
    assert f"PL{output_egrid_file.OZONE_COMBUSTION_GENERATION}" not in rates


def test_fuel_type_rates_include_waste_in_fossil():
    annual = pd.DataFrame(
        {
            "net_generation_mwh": [100.0, 50.0, 0.0],
            "fuel_consumed_for_electricity_mmbtu": [1000.0, 800.0, 10.0],
            "co2_mass_lb_for_electricity_adjusted": [100_000.0, 50_000.0, 1000.0],
        },
        index=pd.Index(["natural_gas", "waste", "petroleum"], name="fuel_category"),
    )
    for column in output_egrid_file.FUEL_TYPE_RATE_MASSES.values():
        if column not in annual:
            annual[column] = 0.0
    rates = output_egrid_file.calculate_fuel_type_rates(annual, annual)
    assert rates["FSC2RT"] == 151_000 / 150
    assert rates["FSC2R"] == 151_000 / 1810
    # oil has emissions but no generation, so its output rate is blank
    assert pd.isna(rates["OCO2RT"])
    assert rates["OCO2R"] == 100.0
    assert pd.isna(rates["CCO2RT"])


def test_plant_adjustment_values_and_sources(results):
    plants = output_egrid_file.build_plant_table(2024, "2024/").set_index("ORISPL")
    gas = plants.loc[1]
    # adjustments are differences between OGE's unadjusted and adjusted values
    assert gas["BIOCO2"] == 250.0
    assert gas["CHPCO2"] == 250.0
    assert gas["UNCO2"] - gas["BIOCO2"] - gas["CHPCO2"] == gas["PLCO2AN"]
    assert gas["BIONOX"] == 0
    assert gas["CHPNOXOZ"] == 0.001
    assert gas["CHPCHTI"] == 2000.0
    assert gas["CHPCHTIOZ"] == 1000.0
    assert gas["RMBMFLAG"] == "Yes"
    assert gas["CHPFLAG"] == "Yes"
    assert pd.isna(plants.loc[2, "CHPFLAG"])
    # sources list the data sources of the subplant-months with a positive value
    assert gas["UNCO2SRC"] == "EPA/CAPD; EIA"
    assert gas["UNHTISRC"] == "EPA/CAPD; EIA"
    assert gas["UNNOZSRC"] == "EIA"
    assert gas["UNHOZSRC"] == "EIA"
    # sources are blank where the value is blank, like the gas plant's Hg
    assert pd.isna(gas["UNHG"])
    assert pd.isna(gas["UNHGSRC"])
    assert plants.loc[2, "UNHTISRC"] == "EIA"
    assert pd.isna(plants.loc[2, "UNCO2SRC"])
