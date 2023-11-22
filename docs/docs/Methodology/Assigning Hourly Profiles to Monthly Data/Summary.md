---
stoplight-id: hourly_shaping_summary
---

## Assigning an hourly profile to lower-resolution data

One of the novel contributions of the Open Grid Emissions Initiative is a method for assigning hourly profiles to data that are reported at a lower resolution.


### Summary of approaches

The hourly data provided by the OGEI is a combination of measured hourly data reported in CEMS and monthly data reported in EIA-923 which is assigned an imputed hourly profile.

The first step in ensuring complete data coverage is identifying the subplant-months for which complete data is available from either source. If complete CEMS data is available for all hours in a month and for all units in a subplant, the reported CEMS data is used for that subplant month. If the CEMS data for a subplant-month is incomplete or not available, the reported EIA-923 is used, and an hourly profile must be assigned to it.

Prior to identifying the data source for each subplant-month, the pipeline drops any CEMS data for unit-months where zero operation is reported for the entire month. Currently, the pudl data pipeline fills all missing gross load and heat input values with zeros, so we have no way of identifying which data are measured zeros or missing data (although this issue [will be fixed](https://github.com/catalyst-cooperative/pudl/pull/1692) in a future pudl release). Removing these values gives the pipeline the opportunity to fill this missing data using reported EIA-923 values. In theory, if the reported zeros were correct, the reported EIA-923 should also be zero.



The method used to impute the hourly shape of monthly EIA-923 data depends on the availability of any partial data from CEMS. Whenever complete hourly CEMS data is not available, we always try to use plant-specific imputation methods if possible (methods 1 and 2), otherwise we use region- and fuel-specific imputation methods.



1. If complete hourly data is available for only a subset of units that make up a subplant, we scale this partial data to match the reported values in EIA-923 for the entire subplant.
2. If CEMS data for an entire subplant is missing, but there is complete or partial CEMS data available for other subplants at the same plant, we use the hourly profile of these other subplant(s) to shape the EIA-923 data.
3. If no CEMS data is available for any subplants at a plant, we use regional fuel-specific imputation methods to shape the reported EIA-923 data.


### Annually-reported EIA data

Although most generation and fuel data is reported to EIA-923 as monthly totals, certain plants only report annual totals (in 2020, this made up about 10% of the total EIA data). In these cases, EIA imputes a monthly profile for this data based on the monthly patterns of similar plants.

In general, the only implication of this is that the hourly profiles for these plants may be less accurate since they are relying on imputed rather than reported monthly data. However, in the case when one of these annually-reporting plants also only reports to CEMS for part of the year, and thus we use CEMS data for some months and EIA data for others, there is a potential that some emissions may be double-counted or under-counted depending on the accuracy of the EIA’s monthly imputation method. For 2020, this situation applied to less than 1% of the final data, but over time we will work to minimize any potential errors due to annually-reported data. To track progress on this issue see [https://github.com/singularity-energy/open-grid-emissions/issues/170](https://github.com/singularity-energy/open-grid-emissions/issues/170).