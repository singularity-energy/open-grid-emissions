import pytest
from pandas import DataFrame

from src.hourly_egrid import PUDL


@pytest.fixture
def two_plants()->DataFrame:
    df = DataFrame({'plant_id_eia':[589,544], 'plant_name_eia':["J C McNeil","CT plant with data"], 'state':["VT", "CT"]})
    return df

def test_plants_from_ba():
    P = PUDL()
    assert "ISNE" in P.ba_list
    plants = P.plants_from_ba("ISNE")
    assert "J C McNeil" in plants["plant_name_eia"].tolist()
    # J C McNeil plant is in VT
    assert plants[plants["plant_name_eia"] == "J C McNeil"]["state"].array == "VT"

def test_hourly_from_plants(two_plants):
    P = PUDL()
    hourly = P.hourly_from_plants(two_plants, 2018)
    assert set(hourly.plant_id_eia.unique()) == set(two_plants.plant_id_eia)
