import pandas as pd

import oge.load_data as load_data


def load_emissions_controls_helper(years, expect_empty=False):
    for year in years:
        print(f"-- Loading emission controls data from EIA-923 for {year}")
        df_emission_controls = load_data.load_emissions_controls_eia923(year)
        print("Columns:\n", df_emission_controls.columns)
        if expect_empty:
            assert len(df_emission_controls) == 0
        else:
            assert len(df_emission_controls) > 0
        print("Looks good!")


def test_load_emissions_controls_eia923_post_2015():
    load_emissions_controls_helper(list(reversed(range(2016, 2021))))


def test_load_emissions_controls_eia923_2012_to_2015():
    load_emissions_controls_helper(list(reversed(range(2012, 2016))))


def test_load_emissions_controls_eia923_pre_2012():
    """Not available before 2008 because EIA publishes 906/920 instead of 923."""
    load_emissions_controls_helper(list(reversed(range(2008, 2012))), expect_empty=True)


def test_load_diba_data():
    """Make sure that the DIBA data from EIA-930 is loaded without errors."""
    load_data.load_diba_data(2021)
    load_data.load_diba_data(2020)
    load_data.load_diba_data(2019)


def test_load_campd_facility_cleans_unit_ids_and_types(monkeypatch, tmp_path):
    """Leading zeros should be removed from unit IDs, each unit should keep the unit
    type in use at the end of the year, boilers should get their firing type and
    wet/dry bottom type, primary fuels should get an EIA energy source code, and plant
    IDs remapped in the manual crosswalk should get their EIA plant ID."""
    (tmp_path / "epa").mkdir()
    (tmp_path / "epa/facility-2025.csv").write_text(
        '"State","Facility Name","Facility ID","Unit ID","Source Category",'
        '"Unit Type","Primary Fuel Type","Secondary Fuel Type",'
        '"Associated Generators & Nameplate Capacity (MWe)"\n'
        '"AL","Plant A",1,"01","Electric Utility","Tangentially-fired","Coal",,'
        '"1 (100)"\n'
        '"CA","AES Alamitos",315,"CT1","Electric Utility",'
        '"Combustion turbine (Started May 19, 2025)","Pipeline Natural Gas",'
        '"Diesel Oil","CT1 (50), ST1 (25)"\n'
        '"FL","Plant B",2,"4","Electric Utility",'
        '"Integrated gasification combined cycle (Ended Jun 07, 2025), '
        'Combustion turbine (Started Jun 07, 2025)","Natural Gas",,\n'
        '"FL","Plant B",2,"5","Electric Utility","Dry bottom wall-fired boiler",'
        '"Coal, Pipeline Natural Gas",,"5 (300)"\n'
    )
    monkeypatch.setattr(
        load_data, "downloads_folder", lambda path: str(tmp_path / path)
    )

    facility = load_data.load_campd_facility(2025)

    assert facility["emissions_unit_id_epa"].tolist() == ["1", "CT1", "4", "5"]
    assert facility["unit_type_epa"].tolist() == [
        "Tangentially-fired",
        "Combustion turbine",
        "Combustion turbine",
        "Dry bottom wall-fired boiler",
    ]
    assert facility["boiler_firing_type"].tolist() == [
        "tangential_firing",
        pd.NA,
        pd.NA,
        "wall_fired",
    ]
    assert facility["wet_dry_bottom"].tolist() == [pd.NA, pd.NA, pd.NA, "dry"]
    assert facility["energy_source_code_epa"].tolist() == ["SUB", "NG", "NG", "SUB"]
    assert facility["plant_id_eia"].tolist() == [1, 62115, 2, 2]
    assert facility["associated_generators_and_capacity_epa"].tolist()[:2] == [
        "1 (100)",
        "CT1 (50), ST1 (25)",
    ]


def _write_excel_sheets(path, sheets: dict) -> None:
    """Writes each table to its own sheet, below a number of blank title rows."""
    with pd.ExcelWriter(path) as writer:
        for sheet_name, (title_rows, table) in sheets.items():
            table.to_excel(
                writer, sheet_name=sheet_name, startrow=title_rows, index=False
            )


def test_load_raw_eia860m_puerto_rico_generators(monkeypatch, tmp_path):
    """Operating and planned generators should be kept, along with generators that
    retired during the year, with the status code taken from the status description,
    dates built from their month and year columns, and the notes below each table
    dropped."""
    notes = {"Entity ID": "NOTES: capacity under 1 MW is excluded"}
    operating = pd.DataFrame(
        [
            {
                "Entity ID": 1,
                "Entity Name": "Utility A",
                "Plant ID": 61034,
                "Plant Name": "Plant A",
                "Plant State": "PR",
                "Generator ID": "1",
                "Unit Code": "CC1",
                "Nameplate Capacity (MW)": 100.0,
                "Energy Source Code": "NG",
                "Prime Mover Code": "CT",
                "Operating Month": 3,
                "Operating Year": 2000,
                "Planned Retirement Month": " ",
                "Planned Retirement Year": " ",
                "Status": "(SB) Standby/Backup: available for service but not normally used",
            },
            notes,
        ]
    )
    planned = pd.DataFrame(
        [
            {
                "Entity ID": 1,
                "Entity Name": "Utility A",
                "Plant ID": 61034,
                "Plant Name": "Plant A",
                "Plant State": "PR",
                "Generator ID": "BESS",
                "Unit Code": None,
                "Nameplate Capacity (MW)": 50.0,
                "Energy Source Code": "MWH",
                "Prime Mover Code": "BA",
                "Planned Operation Month": 6,
                "Planned Operation Year": 2025,
                "Status": "(V) Under construction, more than 50 percent complete",
            },
            notes,
        ]
    )
    retired_generator = {
        "Entity ID": 2,
        "Entity Name": "Utility B",
        "Plant ID": 61147,
        "Plant Name": "Plant B",
        "Plant State": "PR",
        "Generator ID": "2",
        "Unit Code": None,
        "Nameplate Capacity (MW)": 50.0,
        "Energy Source Code": "RFO",
        "Prime Mover Code": "ST",
        "Operating Month": 2,
        "Operating Year": 1959,
    }
    retired = pd.DataFrame(
        [
            retired_generator | {"Retirement Month": 7, "Retirement Year": 2024},
            retired_generator
            | {"Generator ID": "1", "Retirement Month": 7, "Retirement Year": 2009},
            notes,
        ]
    )
    (tmp_path / "eia860m").mkdir()
    _write_excel_sheets(
        tmp_path / "eia860m/december_generator2024.xlsx",
        {
            "Operating_PR": (2, operating),
            "Planned_PR": (1, planned),
            "Retired_PR": (2, retired),
        },
    )
    monkeypatch.setattr(
        load_data, "downloads_folder", lambda path: str(tmp_path / path)
    )

    generators = load_data.load_raw_eia860m_puerto_rico_generators(2024)

    assert generators["generator_id"].tolist() == ["1", "BESS", "2"]
    assert generators["unit_id_eia"].tolist() == ["CC1", pd.NA, pd.NA]
    assert generators["operational_status_code"].tolist() == ["SB", "V", "RE"]
    assert generators["operational_status"].tolist() == [
        "existing",
        "proposed",
        "retired",
    ]
    assert generators["generator_operating_date"].tolist()[0] == pd.Timestamp(
        "2000-03-01"
    )
    assert pd.isna(generators["planned_generator_retirement_date"].iloc[0])
    assert generators["current_planned_generator_operating_date"].tolist()[
        1
    ] == pd.Timestamp("2025-06-01")
    assert generators["generator_retirement_date"].tolist()[2] == pd.Timestamp(
        "2024-07-01"
    )
    assert (generators["report_date"] == pd.Timestamp("2024-01-01")).all()


def test_load_raw_eia923_puerto_rico_generation_fuel(monkeypatch, tmp_path):
    """Each month should become its own row with PUDL's column names, state increment
    rows should be dropped, and rows reporting the same prime mover and fuel for a
    plant should be combined, with their heat content recalculated."""
    months = pd.date_range("2024-01-01", periods=12, freq="MS").strftime("%B")

    def plant_row(plant_id, fuel, quantity, mmbtu_per_unit, netgen):
        row = {
            "Plant Id": plant_id,
            "Combined Heat And\nPower Plant": "N",
            "Plant Name": f"Plant {plant_id}",
            "Operator Name": "Utility A",
            "Operator Id": 1,
            "Reported\nPrime Mover": "IC",
            "Reported\nFuel Type Code": fuel,
        }
        for month in months:
            row |= {
                f"Quantity\n{month}": quantity,
                f"Elec_Quantity\n{month}": quantity,
                f"MMBtuPer_Unit\n{month}": mmbtu_per_unit,
                f"Tot_MMBtu\n{month}": quantity * mmbtu_per_unit,
                f"Elec_MMBtu\n{month}": quantity * mmbtu_per_unit,
                f"Netgen\n{month}": netgen,
            }
        return row

    page_1 = pd.DataFrame(
        [
            plant_row(61034, "DFO", 10, 5.8, 3),
            plant_row(62410, "NG", 100, 1.0, 10),
            plant_row(62410, "NG", 300, 1.2, 30),
            plant_row(99999, "DFO", 0, 0, 0),
        ]
    )
    (tmp_path / "eia923/f923_2024").mkdir(parents=True)
    filepath = (
        tmp_path / "eia923/f923_2024/EIA923_Schedules_2_3_4_5_M_12_2024_Final.xlsx"
    )
    _write_excel_sheets(filepath, {"Page 1 Puerto Rico": (6, page_1)})
    monkeypatch.setattr(
        load_data, "find_downloaded_file", lambda directory, pattern: str(filepath)
    )
    monkeypatch.setattr(
        load_data,
        "load_pudl_table",
        lambda table_name, columns=None: pd.DataFrame(
            {"code": ["DFO", "NG"], "fuel_type_code_pudl": ["oil", "gas"]}
        ),
    )

    generation_fuel = load_data.load_raw_eia923_puerto_rico_generation_fuel(2024)

    assert len(generation_fuel) == 24
    assert set(generation_fuel["plant_id_eia"]) == {61034, 62410}
    assert generation_fuel["report_date"].min() == pd.Timestamp("2024-01-01")
    assert generation_fuel["report_date"].max() == pd.Timestamp("2024-12-01")
    january = generation_fuel[
        generation_fuel["report_date"] == pd.Timestamp("2024-01-01")
    ].set_index("plant_id_eia")
    assert january.loc[61034, "energy_source_code"] == "DFO"
    assert january.loc[61034, "fuel_consumed_mmbtu"] == 58
    assert january.loc[61034, "fuel_mmbtu_per_unit"] == 5.8
    assert january.loc[62410, "fuel_consumed_units"] == 400
    assert january.loc[62410, "fuel_consumed_mmbtu"] == 460
    assert january.loc[62410, "fuel_mmbtu_per_unit"] == 1.15
    assert january.loc[62410, "net_generation_mwh"] == 40
    assert january["fuel_type_code_pudl"].tolist() == ["oil", "gas"]
    assert (generation_fuel["data_maturity"] == "final").all()


def test_replace_puerto_rico_data():
    """PUDL's Puerto Rico rows should be replaced by the data of each year loaded from
    2018 onward, limited to the columns loaded."""

    def load_puerto_rico_year(year):
        return pd.DataFrame(
            {
                "plant_id_eia": [61034],
                "generator_id": ["CT1"],
                "report_date": [pd.Timestamp(year=year, month=1, day=1)],
                "unit_id_eia": ["CC1"],
            }
        )

    pudl_table = pd.DataFrame(
        {
            "plant_id_eia": [3, 61034],
            "generator_id": ["1", "CT1"],
            "report_date": pd.to_datetime(["2017-01-01", "2026-01-01"]),
        }
    )

    all_years = load_data.replace_puerto_rico_data(
        pudl_table, load_puerto_rico_year, year=None, end_year=None, columns=None
    )
    assert all_years["report_date"].dt.year.tolist() == [2017] + list(
        range(2018, load_data.latest_validated_year + 1)
    )

    some_years = load_data.replace_puerto_rico_data(
        pudl_table[["plant_id_eia", "generator_id"]],
        load_puerto_rico_year,
        year=2016,
        end_year=2019,
        columns=["plant_id_eia", "generator_id", "unit_id_pudl"],
    )
    assert some_years.columns.tolist() == [
        "plant_id_eia",
        "generator_id",
        "unit_id_pudl",
    ]
    assert some_years["plant_id_eia"].tolist() == [3, 61034, 61034]

    before_2018 = load_data.replace_puerto_rico_data(
        pudl_table.iloc[[0]],
        load_puerto_rico_year,
        year=2010,
        end_year=2015,
        columns=None,
    )
    assert before_2018["plant_id_eia"].tolist() == [3]


def test_load_puerto_rico_generators_fills_late_reporters(monkeypatch):
    """Plants that report generation but are missing from the December EIA-860M should
    get the generators from PUDL's snapshot that operated during the year."""
    monkeypatch.setattr(
        load_data,
        "load_raw_eia860m_puerto_rico_generators",
        lambda year: pd.DataFrame(
            {
                "plant_id_eia": [61034],
                "generator_id": ["CT1"],
                "energy_source_code_1": ["NG"],
                "operational_status_code": ["OP"],
                "operational_status": ["existing"],
                "generator_operating_date": pd.to_datetime(["2000-01-01"]),
                "generator_retirement_date": pd.to_datetime([pd.NaT]),
                "unit_id_eia": ["CC1"],
                "report_date": pd.to_datetime(["2022-01-01"]),
            }
        ),
    )
    monkeypatch.setattr(
        load_data,
        "load_raw_eia923_puerto_rico_generation_fuel",
        lambda year: pd.DataFrame({"plant_id_eia": [61034, 66319]}),
    )
    snapshot = pd.DataFrame(
        {
            "plant_id_eia": [66319, 66319, 66319, 66319],
            "generator_id": ["OLD", "RETIRED_LATER", "NEW", "RETIRED_EARLIER"],
            "energy_source_code_1": ["DFO", "DFO", "SUN", "DFO"],
            "operational_status_code": ["SB", "RE", "OP", "RE"],
            "operational_status": ["existing", "retired", "existing", "retired"],
            "generator_operating_date": pd.to_datetime(
                ["1991-01-01", "2004-06-01", "2023-06-01", "1990-01-01"]
            ),
            "generator_retirement_date": pd.to_datetime(
                [pd.NaT, "2025-04-01", pd.NaT, "2020-01-01"]
            ),
        }
    )
    fuel_codes = pd.DataFrame(
        {"code": ["DFO", "NG", "SUN"], "fuel_type_code_pudl": ["oil", "gas", "solar"]}
    )
    for loader in ["load_pudl_table_from_local", "load_pudl_table_from_s3"]:
        monkeypatch.setattr(
            load_data, loader, lambda table_name, columns=None: snapshot
        )
    monkeypatch.setattr(
        load_data, "load_pudl_table", lambda table_name, columns=None: fuel_codes
    )
    load_data.load_puerto_rico_generators.cache_clear()

    generators = load_data.load_puerto_rico_generators(2022).set_index("generator_id")
    load_data.load_puerto_rico_generators.cache_clear()

    assert generators.index.tolist() == ["CT1", "OLD", "RETIRED_LATER"]
    assert generators.loc["RETIRED_LATER", "operational_status_code"] == "OP"
    assert pd.isna(generators.loc["RETIRED_LATER", "generator_retirement_date"])
    assert (generators["report_date"] == pd.Timestamp("2022-01-01")).all()
    assert generators["fuel_type_code_pudl"].tolist() == ["gas", "oil", "oil"]
    assert generators["fuel_type_count"].tolist() == [1, 1, 1]


def test_load_mats_data_before_mats():
    """There is no MATS data before 2015, so an empty table is returned."""
    mats = load_data.load_mats_data(2014)

    assert mats.empty
    assert mats.columns.tolist() == [
        "plant_id_epa",
        "emissions_unit_id_epa",
        "datetime_utc",
        "hg_mass_lb",
        "hg_mass_measurement_code",
    ]
    assert mats["hg_mass_lb"].dtype == "float64"
