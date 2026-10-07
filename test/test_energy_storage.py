import numpy as np
import pandas as pd

import oge.energy_storage as energy_storage


def test_assign_storage_category_ignores_supported_plants_without_generators():
    # battery 1 supports a plant with an operating solar generator, while battery 2
    # supports a plant that only has a proposed solar generator and a battery
    storage_generators = pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "generator_id": ["BESS1", "BESS2"],
            "prime_mover_code": ["BA", "BA"],
            "energy_source_code_1": ["MWH", "MWH"],
            "is_dc_coupled_tightly": [False, False],
            "is_direct_support": [True, True],
            "served_co_located_renewable_firming": [False, False],
            "is_independent": [False, False],
            "plant_id_eia_direct_support_1": [10, 20],
            "plant_id_eia_direct_support_2": [np.nan, np.nan],
            "plant_id_eia_direct_support_3": [np.nan, np.nan],
        }
    )
    generators = pd.DataFrame(
        {
            "plant_id_eia": [1, 2, 10, 20, 20],
            "generator_id": ["BESS1", "BESS2", "PV1", "PV2", "BESS3"],
            "prime_mover_code": ["BA", "BA", "PV", "PV", "BA"],
            "energy_source_code_1": ["MWH", "MWH", "SUN", "SUN", "MWH"],
            "is_operating": [True, True, True, False, True],
            "latitude": [30.0, 31.0, 32.0, 33.0, 33.0],
            "longitude": [-100.0, -101.0, -102.0, -103.0, -103.0],
        }
    )

    result = energy_storage.assign_storage_category_to_generators(
        storage_generators, generators
    ).set_index("plant_id_eia")

    # battery 1 is co-located with the solar plant it supports
    assert result.loc[1, "storage_category_method"] == "direct_support_other_plant"
    assert result.loc[1, "co_located_plant_id_list"] == [10]
    # battery 2 is still reported as providing direct support, but the plant it
    # supports does not have an operating non-storage generator, so it is not recorded
    # as the co-located plant
    assert result.loc[2, "storage_category_method"] == "eia860_flag"
    assert result.loc[2, "storage_category"] == "co_located"
    assert result.loc[2, "co_located_plant_id_list"] == []
