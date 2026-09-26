from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

import oge.energy_storage as energy_storage


@pytest.fixture
def energy_storage_table():
    """A fake `core_eia923__monthly_energy_storage` table.

    Includes a battery, a solar array at the same plant (which should be excluded), a
    compressed air storage plant that reports natural gas consumption, and a pumped
    storage plant with one month that is missing discharge data.
    """
    return pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 2, 3, 3, 3],
            "report_date": pd.to_datetime(
                [
                    "2024-01-01",
                    "2024-01-01",
                    "2024-01-01",
                    "2024-01-01",
                    "2024-02-01",
                    "2024-03-01",
                ]
            ),
            "prime_mover_code": ["BA", "PV", "CE", "PS", "PS", "PS"],
            "energy_source_code": ["MWH", "SUN", "NG", "WAT", "WAT", "WAT"],
            "fuel_units": ["mwh", None, "mcf", "mwh", "mwh", "mwh"],
            "fuel_consumed_for_electricity_units": [100.0, np.nan, 3000.0, 50, 60, 70],
            "gross_generation_mwh": [85.0, 500.0, 200.0, 40, 48, 0],
            "net_generation_mwh": [-15.0, 490.0, -300.0, -10, -12, -14],
        }
    )


@pytest.fixture
def storage_dispatch():
    """Monthly charging and discharging data, as returned by
    `load_energy_storage_dispatch()`."""
    return pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 4, 4, 5, 5],
            "report_date": pd.to_datetime(
                [
                    "2024-01-01",
                    "2024-02-01",
                    "2024-01-01",
                    "2024-02-01",
                    "2024-01-01",
                    "2024-02-01",
                ]
            ),
            "prime_mover_code": ["BA", "BA", "PS", "PS", "PS", "PS"],
            "energy_source_code": ["MWH", "MWH", "WAT", "WAT", "WAT", "WAT"],
            "net_generation_mwh": [-15.0, -20.0, 50.0, 30.0, -20.0, -25.0],
            "discharge_filled": [False] * 6,
            "storage_charge_mwh": [100.0, 120.0, 50.0, 70.0, 100.0, 125.0],
            "storage_discharge_mwh": [85.0, 100.0, 100.0, 100.0, 80.0, 100.0],
        }
    )


def test_load_energy_storage_dispatch_excludes_non_storage(energy_storage_table):
    with patch.object(
        energy_storage.load_data,
        "load_pudl_table",
        return_value=energy_storage_table,
    ):
        result = energy_storage.load_energy_storage_dispatch(2024)

    assert "PV" not in result["prime_mover_code"].values
    battery = result[result["prime_mover_code"] == "BA"].iloc[0]
    assert battery["storage_charge_mwh"] == 100
    assert battery["storage_discharge_mwh"] == 85


def test_load_energy_storage_dispatch_compressed_air(energy_storage_table):
    with patch.object(
        energy_storage.load_data,
        "load_pudl_table",
        return_value=energy_storage_table,
    ):
        result = energy_storage.load_energy_storage_dispatch(2024)

    # charge is gross generation minus net generation, rather than the natural gas
    # consumed
    compressed_air = result[result["prime_mover_code"] == "CE"].iloc[0]
    assert compressed_air["storage_discharge_mwh"] == 200
    assert compressed_air["storage_charge_mwh"] == 500


def test_load_energy_storage_dispatch_fills_missing_discharge(energy_storage_table):
    with patch.object(
        energy_storage.load_data,
        "load_pudl_table",
        return_value=energy_storage_table,
    ):
        result = energy_storage.load_energy_storage_dispatch(2024)

    pumped_storage = result[result["prime_mover_code"] == "PS"].set_index("report_date")
    # March reports zero discharge, but all other months are consistent, so the
    # discharge is filled with the net generation plus the charge
    assert pumped_storage.loc["2024-03-01", "discharge_filled"]
    assert pumped_storage.loc["2024-03-01", "storage_discharge_mwh"] == 56
    assert not pumped_storage.loc["2024-01-01", "discharge_filled"]


def test_load_energy_storage_dispatch_does_not_fill_if_other_months_inconsistent(
    energy_storage_table,
):
    # make February inconsistent for the pumped storage plant
    energy_storage_table.loc[4, "net_generation_mwh"] = -100
    with patch.object(
        energy_storage.load_data,
        "load_pudl_table",
        return_value=energy_storage_table,
    ):
        result = energy_storage.load_energy_storage_dispatch(2024)

    pumped_storage = result[result["prime_mover_code"] == "PS"].set_index("report_date")
    assert not pumped_storage["discharge_filled"].any()
    assert pumped_storage.loc["2024-03-01", "storage_discharge_mwh"] == 0


def test_identify_pumped_storage_with_inflow(storage_dispatch):
    # add a pumped storage plant that did not charge or discharge
    storage_dispatch = pd.concat(
        [
            storage_dispatch,
            pd.DataFrame(
                {
                    "plant_id_eia": [6],
                    "prime_mover_code": ["PS"],
                    "storage_charge_mwh": [0.0],
                    "storage_discharge_mwh": [0.0],
                }
            ),
        ],
        ignore_index=True,
    )

    result = energy_storage.identify_pumped_storage_with_inflow(storage_dispatch)

    # plant 4 discharges more than it charges, plant 5 discharges less than it charges,
    # and plant 6 does not discharge. Plant 1 is not pumped storage
    assert result == [4]


def test_allocate_energy_storage_dispatch_to_subplants_by_capacity(storage_dispatch):
    primary_fuel_table = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 4, 5, 5],
            "subplant_id": [1, 2, 1, 1, 2],
            "generator_id": ["B1", "B2", "P1", "P1", "P2"],
        }
    )
    generators = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 4, 5, 5],
            "generator_id": ["B1", "B2", "P1", "P1", "P2"],
            "prime_mover_code": ["BA", "BA", "PS", "PS", "PS"],
            "capacity_mw": [75.0, 25.0, 100.0, 0.0, 0.0],
        }
    )
    with patch.object(
        energy_storage.load_data, "load_pudl_table", return_value=generators
    ):
        result = energy_storage.allocate_energy_storage_dispatch_to_subplants(
            storage_dispatch, primary_fuel_table, 2024
        ).set_index(["plant_id_eia", "subplant_id", "report_date"])

    # plant 1 is allocated based on capacity
    assert result.loc[(1, 1, "2024-01-01"), "storage_charge_mwh"] == 75
    assert result.loc[(1, 2, "2024-01-01"), "storage_charge_mwh"] == 25
    # plant 5 does not report capacity, so is allocated equally
    assert result.loc[(5, 1, "2024-01-01"), "storage_charge_mwh"] == 50
    assert result.loc[(5, 2, "2024-01-01"), "storage_charge_mwh"] == 50
    # the total data is unchanged
    assert (
        result["storage_charge_mwh"].sum()
        == storage_dispatch["storage_charge_mwh"].sum()
    )


def test_create_monthly_energy_storage_data_blanks_pumped_storage_with_inflow(
    storage_dispatch,
):
    primary_fuel_table = pd.DataFrame(
        {
            "plant_id_eia": [1, 4, 5],
            "subplant_id": [1, 1, 1],
            "generator_id": ["B1", "P1", "P1"],
            "subplant_storage_category_method": [
                "no_evidence",
                "pumped_storage_with_inflow",
                "same_plant",
            ],
        }
    )
    generators = pd.DataFrame(
        {
            "plant_id_eia": [1, 4, 5],
            "generator_id": ["B1", "P1", "P1"],
            "prime_mover_code": ["BA", "PS", "PS"],
            "capacity_mw": [100.0, 100.0, 100.0],
        }
    )
    with (
        patch.object(
            energy_storage,
            "load_energy_storage_dispatch",
            return_value=storage_dispatch,
        ),
        patch.object(
            energy_storage.load_data, "load_pudl_table", return_value=generators
        ),
    ):
        result = energy_storage.create_monthly_energy_storage_data(
            primary_fuel_table, 2024
        )

    pumped_storage_with_inflow = result[result["plant_id_eia"] == 4]
    assert pumped_storage_with_inflow["storage_discharge_mwh"].isna().all()
    assert (pumped_storage_with_inflow["storage_charge_mwh"] == [50, 70]).all()
    other_storage = result[result["plant_id_eia"] != 4]
    assert other_storage["storage_discharge_mwh"].notna().all()


def test_add_monthly_energy_storage_data_keeps_non_storage_blank():
    monthly_subplant_data = pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "subplant_id": [1, 1],
            "report_date": pd.to_datetime(["2024-01-01", "2024-01-01"]),
            "net_generation_mwh": [-15.0, 500.0],
        }
    )
    monthly_storage_data = pd.DataFrame(
        {
            "plant_id_eia": [1],
            "subplant_id": [1],
            "report_date": pd.to_datetime(["2024-01-01"]),
            "storage_charge_mwh": [100.0],
            "storage_discharge_mwh": [85.0],
        }
    )

    result = energy_storage.add_monthly_energy_storage_data(
        monthly_subplant_data, monthly_storage_data
    ).set_index("plant_id_eia")

    assert len(result) == 2
    assert result.loc[1, "storage_charge_mwh"] == 100
    assert np.isnan(result.loc[2, "storage_charge_mwh"])


def test_add_monthly_energy_storage_data_adds_missing_months():
    monthly_subplant_data = pd.DataFrame(
        {
            "plant_id_eia": [1],
            "subplant_id": [1],
            "report_date": pd.to_datetime(["2024-01-01"]),
            "net_generation_mwh": [-15.0],
        }
    )
    monthly_storage_data = pd.DataFrame(
        {
            "plant_id_eia": [1, 1, 2],
            "subplant_id": [1, 1, 1],
            "report_date": pd.to_datetime(["2024-01-01", "2024-02-01", "2024-01-01"]),
            "storage_charge_mwh": [100.0, 50.0, 20.0],
            "storage_discharge_mwh": [85.0, 40.0, 15.0],
        }
    )

    result = energy_storage.add_monthly_energy_storage_data(
        monthly_subplant_data, monthly_storage_data
    ).set_index(["plant_id_eia", "report_date"])

    # February is added for plant 1, since plant 1 is in the subplant data, but plant 2
    # is not added
    assert len(result) == 2
    assert result.loc[(1, "2024-02-01"), "storage_charge_mwh"] == 50
    assert np.isnan(result.loc[(1, "2024-02-01"), "net_generation_mwh"])
