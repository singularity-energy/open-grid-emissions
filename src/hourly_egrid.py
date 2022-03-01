import pandas as pd
from pandas import DataFrame
import sqlalchemy as sa
import os
from typing import List
import dask.dataframe as dd
import numpy as np

"""
    PUDL
Class for interacting with PUDL EIA data
(excluding CEMS data, for which we will use parquet)
"""
class PUDL(object):
    # TODO allow configuration of data sources
    PARQUET_SOURCE = "data/pudl-*/pudl_data/parquet/epacems/year={}/state={}/*.parquet"
    PUDL_SOURCE = "data/pudl-v0.5.0-2021-11-14/pudl_data/sqlite/pudl.sqlite"

    def __init__(self):
        # Note: pudl_db needs to be full path to pudl .sqlite file
        pudl_db = "sqlite:////"+os.getcwd()+"/"+self.PUDL_SOURCE
        self.pudl_engine = sa.create_engine(pudl_db)

        # Get distinct BAs so we can test against valid options later
        sql = """
            SELECT DISTINCT
                plants.balancing_authority_code_eia
            FROM
                plants_entity_eia AS plants
        """
        ba_list = pd.read_sql(sql, self.pudl_engine)
        self.ba_list = set(ba_list["balancing_authority_code_eia"])

    """
        PUDL.plants_from_ba(ba: str)

    According to eGrid documentation (2020, pg 30):
    "Balancing authority ID codes are assigned to a plant based on the EIA-860
    plant-level data and the balancing authority names are assigned to the
    corresponding balancing authority ID codes based on the EIA-861"

    EIA-860 data is availible via PUDL
    """
    def plants_from_ba(self, ba:str)->DataFrame:
        # Ensure that BA is valid.
        assert ba in self.ba_list

        sql = f"""
        SELECT
           plants.plant_id_eia,
           plants.plant_name_eia,
           plants.state
        FROM
           plants_entity_eia AS plants
        WHERE
           plants.balancing_authority_code_eia == \"{ba}\"
        """
        return pd.read_sql(sql, self.pudl_engine)

    """
        PUDL.hourly_from_plants(plants: pandas, year: int)

    Given a pandas dataset describing EIA ID, name, and state of plants
    (see plants_from_ba(.)), return CEMS data for all availible plants in year
    """
    def hourly_from_plants(self, plants:DataFrame, year:int):
        collected = DataFrame()
        # One state at a time
        for state in plants.state.unique():
            state_dat = dd.read_parquet(self.PARQUET_SOURCE.format(year,state))
            # Only select data that actually has either co2 or gross load
            state_dat = state_dat[np.logical_or(~np.isnan(state_dat["co2_mass_tons"]), \
            (state_dat["gross_load_mw"] > 0))]
            # Only select data with target IDs
            state_dat = state_dat[state_dat["plant_id_eia"].isin(plants.plant_id_eia)]
            collected = pd.concat([collected, state_dat.compute()])

        return collected

"""
   HourlyEGRID

Class for holding higher-level logic around hourly calculations
"""
class HourlyEGrid(object):
    def __init__(self, year: int, BA: str):
        self.year = year
        self.ba = BA

    def calc_hourly() -> List[float]:
        # Get plants

        # Get hourly data for every unit associated with each plant
        # Sum over all units reported in each hour
        return []
