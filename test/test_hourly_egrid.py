import pytest

from src.hourly_egrid import PUDL

## TODO create small example DB to live in test directory
@pytest.fixture
def db_path() -> str:
    return "sqlite://///Users/gailin.pease/singularity/hourly-egrid/data/pudl-v0.5.0-2021-11-14/pudl_data/sqlite/pudl.sqlite"

def test_plants_from_ba(db_path):
    P = PUDL(db_path)
    assert "ISNE" in P.ba_list
    plants = P.plants_from_ba("ISNE")
    assert "J C McNeil" in plants["plant_name_eia"].tolist()
    # J C McNeil plant is in VT
    assert plants[plants["plant_name_eia"] == "J C McNeil"]["state"].array == "VT"
