---
stoplight-id: cleaning_cems
---

## Loading data from PUDL
All of the CEMS data used by OGE is loaded from pre-cleaned versions of these files created by the PUDL project. The PUDL transformation process cleans the input data so that it is adjusted for uniformity, corrected for errors, and ready for bulk programmatic use.

A comprehensive list of transformations made by by pudl can be found [here](https://catalystcoop-pudl.readthedocs.io/en/latest/data_sources/epacems.html#pudl-data-transformations) , but the most notable data cleaning that pudl performs is:
 - converts all datetimes to UTC
 - Corrects anomalous `gross_load_mw` values (if a gross load value is greater than 2,000 MW, they assume that they have accidentally reported kWh, so the gross load value is divided by 1,000). This may cause errors if a single unit greater than 2GW of nameplate capacity is ever built.
 - Harmonizes plant ID codes between the EPA version and the EIA version (For some plants, the EPA uses different plant codes from the ones used by the EI.), based on the EPA-EIA power sector data crosswalk (Note: The PSDC does not contain the most up to date plant ID crosswalk. We supplement this matching using the `epa_eia_crosswalk_manual.csv` reference table in OGE)

After loading the transformed cems data into the pipeline, we perform several other data cleaning steps:
- Convert `co2_mass_tons` to `co2_mass_lb` so that all emissions mass units are standardized.
- Add a "report date" to each data point to enable crosswalking with monthly EIA-923 data. This report date reflects the month of the local prevailing time for each plant's reported hourly data.

## Data removed from the CEMS data

Certain data reported to CEMS are removed, including:
- Data for plants that are not connected to the electrical grid (These include all plants with a `plant_id_eia` of `88xxxx` and the plants in [this table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/plants_not_connected_to_grid.csv))
- Data for plants in Puerto Rico (These data will be added in a future release. Track progress [here](https://github.com/singularity-energy/open-grid-emissions/issues/79))
- Data for any unit-months where there is less than one day of non-missing data are dropped. This is generally months that have a single hour of data reported.


## Imputing missing emissions data
Although CEMS reports measured hourly CO2, NOx, and SO2 emissions, in some cases these data are missing. In addition, there is no data reported for CH4 and N2O.

In addition to actual missing values, CO2 data is considered to be missing if reported CO2 is zero and fuel consumption is positive. Currently, only missing CO2 values from CEMS are imputed; imputation of missing NOx and SO2 values will be included in a future release (track progress [here](https://github.com/singularity-energy/open-grid-emissions/issues/153)).

Before imputing the missing data, we must assign an `energy_source_code` to each observation so that the appropriate emission factor can be used. See "[Assigning energy source codes](../Emissions%20Calculations/Assigning%20Energy%20Source%20Codes.md)" for a description of this methodology.

When imputing missing CO2 data, if a unit has a unit-specific fuel type identified by the EPA-EIA crosswalk, we calculate co2 using a fuel-specific emission factor. If not (often in the cases when a unit burns multiple fuels), we fill missing data by using a plant-month weighted average emission factor of all fuels burned in that plant-month.

Hourly CH4 and N2O emissions are then calculated based on hourly fuel consumption and fuel type (see [here](../Emissions%20Calculations/GHG%20Emissions.md) for more details)

## Additional adjustments

Biomass adjusted emissions are calculating using the methodology described [here](../Emissions%20Calculations/Adjusting%20Emissions%20for%20Biomass.md).

Because our [methodology for calculating CHP-adjusted emissions](../Emissions%20Calculations/Adjusting%20Emissions%20for%20CHP.md) relies on net generation data, this adjustment is not calculated until after hourly gross generation has been [converted to hourly net generation](../Converting%20Gross%20to%20Net%20Generation.md).

## Validation Checks

While cleaning the CEMS data, we run several validation checks to ensure that the cleaned data does not contain any unexpected or anomalous values. These checks include:
 - Check that there are no missing energy source codes associated with records with non-zero fuel consumption (which would result in missing emissions data).
 - Check that there are no negative fuel consumption or emissions values
 - Check that all records have been assigned a `subplant_id`
 - Check that adjusted emissions values are less than or equal to total emissions values
 - Check that there are no duplicate timestamps for each subplant

