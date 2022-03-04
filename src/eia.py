"""
Connection management for EIA endpoints
"""

import requests
import os
from os.path import exists
from typing import Dict
import pandas as pd
from pandas import Series
from dateutil.parser import parse as parse_dt
import csv

"""
    EIA
Class for interfacing with EIA API.
Caches results in eia930.

Function naming patterns: Each data access function x has a corresponding _x
x functions check cache, and call _x if needed data is not cached.
_x functions actually fetch from EIA API
"""
class EIA:
    # Cache
    CACHE_PATH = "hourly-egrid/data/eia930/" # created in setup.py
    REGION_F = "regions.csv"
    # API
    BASE_URL_CATEGORY = "https://api.eia.gov/category/?api_key={}&category_id={}"
    BASE_URL_SERIES = "https://api.eia.gov/series/?api_key={}&series_id={}"
    key:str
    # Data
    regions:Dict[str,str]

    ## TODO pull out API key
    def __init__(self):
        self.key = os.getenv('EIA_API_KEY')
        # Use cwd to build path so running from any dir (examples or top level dir) works
        self.cache = os.getcwd().split("hourly-egrid")[0]+"/"+self.CACHE_PATH
        self.get_regions()

    """
        get_regions()
    Get valid regions and BAs, and EIA series ids for each.
    Should only call once on set-up, after that access via self.regions
    """
    def get_regions(self)->None:
        region_f = self.cache + self.REGION_F
        if not exists(region_f):
            self.regions = self._get_regions()
            with open(region_f, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames = self.regions.keys)
                writer.writeheader()
                writer.writerow(self.regions)
        else:
            with open(region_f, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader: # should only be one
                    self.regions = row

        return

    """
        _get_regions
    Get valid regions and BAs from the EIA API, and EIA series ids for each
    """
    def _get_regions(self)->Dict[str,str]:
        # URL where net generation regions are listed
        url = self.BASE_URL_CATEGORY.format(self.key, "2122629")
        r = requests.get(url)

        # Save resulting regions for later use
        map = dict[str, int]()
        jsonres = r.json()
        for cat in jsonres["category"]["childcategories"]:
            name = cat["name"][cat["name"].find("(")+1:cat["name"].find(")")]

            # Category ID for each region contains two series IDs: UTC and local
            # look up UTC series id.
            url = self.BASE_URL_CATEGORY.format(self.key, cat["category_id"])
            r = requests.get(url)
            for option in r.json()["category"]["childseries"]:
                if "UTC" in option["name"]:
                    series_id = option["series_id"]
                    break
            if series_id is None:
                print(f"Could not find series for region {eia_region_code}")
            else:
                map[name] = series_id
        return map

    """
        get_net_generation
    Fetch a year of net generation
    Time window is always one year for ease of interfacing with PUDL data analysis.
    (could easily change this)
    """
    def get_net_generation(self, eia_region_code:str, year:int)->Series:
        assert eia_region_code in self.regions.keys()
        series_f = self.cache + self.regions[eia_region_code]
        if not exists(series_f):
            series = self._get_net_generation(eia_region_code)
            series.to_csv(series_f)
        else:
            series = pd.read_csv(series_f, index_col=0, parse_dates=True)
            series = series[series.columns[0]] # read_csv returns a DF, want series

        return series[f"{year}-01-01T00:00:00Z":f"{year}-12-30T23:59:59Z"]

    """
        _get_net_generation
    Fetch a non-cached year of net generation for a region or BA.
    """
    def _get_net_generation(self, eia_region_code:str)->Series:
        series_id = self.regions[eia_region_code]

        # Fetch time series data
        url = self.BASE_URL_SERIES.format(self.key, series_id)
        r = requests.get(url)

        # Turn dates and generation numbers into a time series
        datajson = r.json()
        index = [parse_dt(datajson["series"][0]["data"][i][0]) \
        for i in range(len(datajson["series"][0]["data"]))]
        data = [datajson["series"][0]["data"][i][1] for i in range(len(datajson["series"][0]["data"]))]

        toreturn = Series(index=index, data=data)
        toreturn.sort_index(inplace=True)
        return toreturn
