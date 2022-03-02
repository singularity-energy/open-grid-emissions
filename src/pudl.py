import pandas as pd
from pandas import DataFrame, Series
import sqlalchemy as sa
from sqlalchemy.engine import Engine
import os
from typing import List, Set, NamedTuple
import dask.dataframe as dd
import numpy as np

"""
   PUDLConfig
parquet_source: Path (full or relative) to partitioned parquet pudl db
    (formattable with year and state)
pudl_source: full path to sqlite pudl db, prefixed by sqlite://// for creation
    of sqlalchemy engine connection
"""
class PUDLConfig(NamedTuple):
    parquet_source: str
    pudl_source:str

"""
   default_config()->PUDLConfig
Default config. Assumes working directory is top level hourly_egrid directory
For other working directories, explicitly specify config when constructing PUDL()
"""
def default_config()->PUDLConfig:
    PARQUET_SOURCE = "data/pudl-*/pudl_data/parquet/epacems/year={}/state={}/*.parquet"
    PUDL_SOURCE = "data/pudl-v0.5.0-2021-11-14/pudl_data/sqlite/pudl.sqlite"
    pudl_db = "sqlite:////"+os.getcwd()+"/"+PUDL_SOURCE
    return PUDLConfig(PARQUET_SOURCE, pudl_db)

"""
    PUDL
Class for interacting with PUDL EIA data
(excluding CEMS data, for which we will use parquet)
"""
class PUDL(object):
    # TODO allow configuration of data sources
    config:PUDLConfig
    pudl_engine:Engine
    ba_list:Set[str]


    def __init__(self, config:PUDLConfig=default_config()):
        self.config = config

        self.pudl_engine = sa.create_engine(self.config.pudl_source)

        # Get distinct BAs so we can test against valid options later
        sql = """
            SELECT DISTINCT
                plants.balancing_authority_code_eia
            FROM
                plants_entity_eia AS plants
            WHERE
                plants.balancing_authority_code_eia IS NOT NULL
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
        AND
           plants.plant_id_eia IS NOT NULL
        AND
           plants.plant_name_eia IS NOT NULL
        AND
           plants.state IS NOT NULL
        """
        return pd.read_sql(sql, self.pudl_engine)

    """
        PUDL.hourly_from_plants(plants: pandas, year: int)

    Given a pandas dataset describing EIA ID, name, and state of plants
    (see plants_from_ba(.)), return CEMS data for all availible plants in year
    """
    def hourly_from_plants(self, plants:DataFrame, year:int)->DataFrame:
        collected = DataFrame()
        # One state at a time
        for state in plants.state.unique():
            state_dat = dd.read_parquet(self.config.parquet_source.format(year,state))
            # Only select data that actually has either co2 or gross load
            state_dat = state_dat[np.logical_or(~np.isnan(state_dat["co2_mass_tons"]), \
            (state_dat["gross_load_mw"] > 0))]
            # Only select data with target IDs
            state_dat = state_dat[state_dat["plant_id_eia"].isin(plants.plant_id_eia)]
            collected = pd.concat([collected, state_dat.compute()])

        return collected
