import pytest

from src.hourly_egrid import PUDL

## TODO create small example DB to live in test directory
@pytest.fixture
def db_path() -> str:
    return "sqlite://///Users/gailin.pease/singularity/hourly-egrid/data/pudl-v0.5.0-2021-11-14/pudl_data/sqlite/pudl.sqlite"

def test_plants_from_ba(db_path):
    P = PUDL(db_path)
    assert "ISNE" in P.ba_list
    plants, names = P.plants_from_ba("ISNE")
    assert len(plants) == len(names)
    assert "J C McNeil" in names
