import pandas as pd
from pandas import DataFrame

"""
   calc_hourly()

Return a dataframe containing time, emissions, power, emission rate for every hour in year
Only includes units that report to CEMS during the hours they report to CEMS

"""
def calc_hourly(pudl:PUDL, ba:str, year:int) -> DataFrame:
    target_cols = ["co2_mass_tons", "gross_load_mw"] # CEMS data of interest, PUDL names

    # Get plants
    plants = pudl.plants_from_ba(ba)

    # Get hourly data of interest for every unit associated with each plant
    hourly = pudl.hourly_from_plants(plants, year)
    hourly = hourly.loc[:,target_cols+["operating_datetime_utc"]]
    hourly = hourly.dropna() # To get cems-only rate, need both Co2 and load

    # Sum over all units reported in each hour
    hours = DataFrame(data=None, columns=target_cols+["emission_rate_co2"],\
    index=sorted(hourly.operating_datetime_utc.unique()), \
    dtype=float)
    for hour in hours.index:
        for col in target_cols:
            hours.loc[hour,col] = sum(hourly.loc[hourly.operating_datetime_utc==hour,col])
        hours.loc[hour,"emission_rate_co2"] = hours.loc[hour,"co2_mass_tons"]*2000/hours.loc[hour,"gross_load_mw"]

    return hours
