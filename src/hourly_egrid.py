import pandas as pd
import sqlalchemy as sa

"""
    PUDL
Class for interacting with PUDL EIA data
(excluding CEMS data, for which we will use parquet)
"""
class PUDL(object):
    def __init__(self, pudl_db:str):
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
    def plants_from_ba(self, ba:str):
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

class HourlyEGrid(object):
    def __init__(self, year: int, BA: str):
        self.year = year
        self.ba = BA
