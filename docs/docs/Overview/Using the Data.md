---
stoplight-id: using_the_data
---

# Summary of the Data

## Data Availability

The latest release includes data for year 2005-2022 covering the contiguous United States, Alaska, and Hawaii. In future releases, we plan to expand the geographic coverage to additional U.S. territories (track progress [here](https://github.com/singularity-energy/open-grid-emissions/issues/79)).

Currently, 2005 is the earliest available year. Note that only monthly and annual data are available for year 2005-2018 because 2019 is the first year for which complete EIA-930 data is available (which is used to assign an hourly profile to non-CEMS data).

## Data Release Schedule

In general, annual data releases will be published in Q4 of the following year (i.e. 2023 data should be published in Q4 of 2024)

Parts of the input data used for the Open Grid Emissions dataset is released by the U.S. Energy Information Administration in the Autumn following the end of each year (2022 data was published in September 2023). Each release will include the most recent year of available data as well as updates of all previous available years based on any updates to the OGE methodology. All previous versions of the data will be archived on Zenodo.

Updated datasets will also be published whenever a new minor version of the open-grid-emissions code is released, usually representing methdological enhancements. These releases may happen mutliple times throughout a year.

## Files

The data is organized into three categories based on the primary use case

* **Carbon accounting data**: BA-level, _consumption-_based emission rates and fuel mix data.
    * Intended user: This data is intended for those wishing to understand the source and emission intensity of the electricity they are consuming.
    * File structure: data for each BA will be located in its own separate file.
* **Power sector data**: BA-level, _production_-based data on generation, emissions, and fuel consumption split out by fuel type.
    * Intended user: This data is intended for analysts who seek to understand the operating characteristics of the power generation fleets in each region.
    * File structure: data for each BA will be located in its own separate file
* **Power plant data**: Plant-level data on hourly power plant operations and static attributes of each plant.
    * Intended user: This data is intended for researchers who seek to understand individual power plants, or who want to perform their own aggregations of the data.
    * File structure: include all dynamic data in a single large csv (perhaps a tar.gz file) and all static attributes (e.g. fuel type, BA, etc) in a separate table, with the plant_id_eia as the primary key

Each of these data categories includes data at different temporal resolutions (hourly, monthly, annual), and in different units (U.S. and Metric). U.S.-unit files uses lb for mass and MMBtu for fuel consumption, while the metric files use kg for mass and GJ for fuel consumption.

<!-- theme: info -->

> #### Carbon accounting data in Hawaii and Alaska
>
> Our carbon accounting data does not include balancing authorities in Hawaii and Alaska, because those BAs do not report their electricity interchange to the EIA. Users interested in carbon accounting in these BAs should use the *power sector data* for their BA of interest.



## Uses and Users

![use_cases.png](https://stoplight.io/api/v1/projects/cHJqOjE1MzAxNA/images/Ax2MofwLkJE)


## Data Dictionary

Generally all column names in the OGE dataset are self-descriptive including the units associated with the data, if relevant. Thus, we do not include a full data dictionary for all files here.

A few notes on column names:
- Data column names generally follow the format `{value_name}_{adjustment_type}_{unit_of_measure}`
- Column names that include the `_for_electricity` adjustment type have been adjusted for CHP (see [methodology](../Methodology/Emissions%20Calculations/Adjusting%20Emissions%20for%20CHP.md)) and represent the portion of the total value that is related to electricity generation.
- Column names that include the `_adjusted` adjustment type have been adjusted for biomass (see [methodology](../Methodology/Emissions%20Calculations/Adjusting%20Emissions%20for%20Biomass.md)).
- Column names that include the `_for_electricity_adjusted` adjustment type have been adjusted for both CHP and biomass.
- The `report_date` column is in the format `MM/01/YYYY` and represents the month `MM` with which the data is associated. So data with a report date of "03/01/2020" is associated with March 2020.

### Plant Static Attributes Table
This table summarizes static attributes of each plant in the dataset.

Field Name | Type | Description | Source
---------|----------|---------|---------
plant_id_eia | integer | The unique six-digit facility identification number, also called an ORISPL, assigned by the Energy Information Administration | EIA-860
plant_name_eia | string | Plant name | EIA-860
capacity_mw | float | Total installed (nameplate) capacity, in megawatts | EIA-860
plant_primary_fuel | string | two- or three-character energy_source_code of the primary fuel consumed by the plant in a year | [Calculated](../Methodology/Data%20Aggregation/Plant%20Primary%20Fuel.md)
fuel_category | string | named fuel category associated with the identified `plant_primary_fuel` | [Mapping table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/energy_source_groups.csv)
fuel_category_eia930 | string | named fuel category used by EIA-930 associated with the identified `plant_primary_fuel` | [Mapping table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/energy_source_groups.csv)
state | string | Two-letter state abbreviation where the plant is physically located | EIA-860
county | string | County name | EIA-860
city | string | City name | EIA-860
ba_code | string | Three- or four-character code representing the commercial balancing authority with which the plant is associated | EIA-860
ba_code_physical | string | three- or four-character code representing the commercial balancing authority with which the plant is associated | [Inferred from EIA-860](../Methodology/Data%20Aggregation/Aggregating%20Data%20to%20Balancing%20Authority.md)
latitude | float | latitude of the plant (in deg.) with Equator as zero point | EIA-860
longitude | float | longitude of the plant (in deg.) measured eastward from Greenwich, UK | EIA-860
plant_operating_date | date | Date the plant began commercial operation | EIA-860
plant_retirement_date | date | Date of the scheduled or effected retirement of the plant | EIA-860
distribution_flag | boolean | Indicates whether a plant is connected to the distribution grid (<= 60kV interconnection) | EIA-860
timezone | string | IANA timezone name for the local timezone where the plant is physically located | EIA-860
data_availability | string | identifies whether this plant only reports to EIA, only reports to CEMS, or reports to both EIA and CEMS | Derived from EIA and CEMS data
shaped_plant_id | integer | A synthetic plant id code used when aggregating plants to the ba-fuel level. If blank, this plant is not aggregated | [See fleet shaping method](../Methodology/Assigning%20Hourly%20Profiles%20to%20Monthly%20Data/Shaping%20Using%20Fleet-Specific%20Profiles.md)

### Plant Metadata Table

This table sumarizes the methods used to assign hourly profiles to each subplant in the OGE dataset, including the aggregated "synthetic plants" representing all non-CEMS-reporting generation in each BA.