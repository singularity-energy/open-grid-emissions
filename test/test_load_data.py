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
