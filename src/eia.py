"""
Connection management for EIA endpoints
"""

import requests
import os
from typing import Dict
import pandas
from pandas import Series
from dateutil.parser import parse as parse_dt

class EIA:
    BASE_URL_CATEGORY = "https://api.eia.gov/category/?api_key={}&category_id={}"
    BASE_URL_SERIES = "https://api.eia.gov/series/?api_key={}&series_id={}&start={}&end={}"
    key:str
    regions:Dict[str,str]

    ## TODO pull out API key
    def __init__(self):
        self.key = "xixW0K8zI3O28PO9IXn8IdULiBZjLLIYezWrcfqk"
        self._get_regions()

    """
        _get_regions
    Get valid regions and BAs ("categories") in EIA terminology. Performed during
    initialization so we can be ready for call requests.
    """
    def _get_regions(self)->None:
        # URL where net generation regions are listed
        url = self.BASE_URL_CATEGORY.format(self.key, "2122629")
        r = requests.get(url)

        # Save resulting regions for later use
        map = dict[str, int]()
        jsonres = r.json()
        for cat in jsonres["category"]["childcategories"]:
            name = cat["name"][cat["name"].find("(")+1:cat["name"].find(")")]
            map[name] = cat["category_id"]
        self.regions = map

    """
        get_net_generation
    Fetch a year of net generation for a region or BA.
    Time window is always one year for ease of interfacing with PUDL data analysis.
    """
    def get_net_generation(self, eia_region_code:str, year:int)->Series:
        assert eia_region_code in self.regions.keys()
        # Look up UTC series IDs
        url = self.BASE_URL_CATEGORY.format(self.key, self.regions[eia_region_code])
        r = requests.get(url)
        for option in r.json()["category"]["childseries"]:
            if "UTC" in option["name"]:
                series_id = option["series_id"]
                break
        if series_id is None:
            throw(AssertionError(f"Could not find series for region {eia_region_code}"))

        # Fetch time series data
        url = self.BASE_URL_SERIES.format(self.key, series_id, f"{year}-01-01", f"{year}-12-30")
        r = requests.get(url)

        # Turn dates and generation numbers into a time series
        datajson = r.json()
        index = [parse_dt(datajson["series"][0]["data"][i][0]) \
        for i in range(len(datajson["series"][0]["data"]))]
        data = [datajson["series"][0]["data"][i][1] for i in range(len(datajson["series"][0]["data"]))]
        return Series(index=index, data=data)
