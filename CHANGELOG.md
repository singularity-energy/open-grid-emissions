CHANGELOG

-------------------------------------------------------------------------------
# Branch partial_cems (PR 2022-06-01) 
-------------------------------------------------------------------------------
Closes #22 and #30

## Partial CEMS calculations
Adds `data_cleaning.scale_partial_cems_data()` as step 7 in the pipeline. This function scales hourly subplant-level CEMS data to match the EIA monthly data totals when only certain CEMS units that make up a subplant report data in a month, and the CEMS totals are less than the reported EIA totals (indicating missing data from CEMS). For any subplants that were flagged as partial CEMS, but for which the CEMS totals are greater than the EIA totals, the function changes the `hourly_data_source` flag back to `cems` instead of `partial_cems`

## Fuel Categories
 - Updates `energy_source_groups.csv` and related functions
 - adds a new column to identigy the fuel categories used in egrid
 - renames columns from "fuel_group_" to "fuel_category_"
 - Both our custom fuel category (which is just labeled `fuel_category`) and `fuel_category_eia930` are now assigned to the EIA-923 and CEMS data during denormalization (step 8)

## Other
 - Implements --small argument for both EIA and CEMS data
 - updates column data checks

-------------------------------------------------------------------------------
# Branch grgmiller/car-342-clean-cems-and-eia923 (PR 2022-05-28) 
-------------------------------------------------------------------------------

Closes #20, partially addresses #17

## Repository Organization
- Moves all validation functions from the `data_pipeline.ipynb` to a separate notebook `data_validation.ipynb`. `data_pipeline` now exports all relevant intermediate files as csvs, that are then read into `data_validation`
- Reorganizes the `data` directory:
    - `data/downloads` contains all files that are downloaded by functions in `load_data`
    - `data/manual` contains all manually-created files, including the egrid static tables
    - `data/outputs` contains intermediate outputs from the data pipeline... any files created by our code that are not final results
    - `data/results` contains all final output files that will be published
- NOTE: results have not yet been generated/tested

## EIA-923 cleaning
- Many generators were missing prime mover codes, which is one of the primary keys used to allocate generation and fuel. These are missing because the PUDL 860 table used to treat a prime mover as a static attribute, so if it changed over time, it became NA. This has been fixed in the PUDL dev branch (https://github.com/catalyst-cooperative/pudl/issues/1585) and in Datasette, but not in the zenodo archive. Thus, we are temporarily adding the function `load_data.download_updated_pudl_database(download=True)` which replaces the zenodo version of `pudl.sqlite` with the Datasette version of the file.

## Emissions
- Adds calculations for N2O, CH4, NOx, and SO2 emissions in addition to CO2 emissions
- There is still some data cleaning that needs to be implemented for these new emissions, but the pipeline for these is substantially complete
- We are not yet implementing calculations of CO2e emissions, because there are a variety of GWPs that can be used to convert CO2, CH4, and N2O to CO2e
- there were certain plant months for which co2 emissions were not calculated. Most of these have a fuel code of "OTH" which has no default emisisons rate. We manually updated these plant fuel codes to OG, since they are refinery plants
- Fixed the function for identifying geothermal emissions factor. Now defaults to using a generator-specific emissions factor

## Emissions adjustments
- Adjust emissions for biomass. There are now three different columns related to emissions:
    - `co2_mass_tons`: total co2 emissions
    - `co2_mass_tons_for_electricity`: `co2_mass_tons` adjusted for CHP
    - `co2_mass_tons_adjusted`: `co2_mass_tons_for_electricity` adjusted for biomass
- Fixes issue where all MSW generation was assigned the MSW fuel code, instead of the MSB and MSN codes, which have different emission rates
- Updates the CHP allocation method to align with the method used in eGRID


## Crosswalking EIA and CEMS at the subplant level
Being able to match data reported in CEMS to data reported in EIA-923 is important for two reasons:
 - it allows us to identify generation and emissions data that is missing from CEMS and needs to be allocated to an hourly profile
 - it allows us to compare gross generation reported in CEMS to net generation reported in EIA-923, and then calculate the relationship between the two so that hourly gross generation in CEMS can be converted to hourly net generation.

However, performing this matchup is challenging for two reasons:
 - EIA and EPA do not always use the same plant ID code for the same plant
 - EPA reports data at the "unit" level, while EIA reports data at the "generator" level, which can be related to each other in 1:1, 1:m, m:1, or m:m relationships.

These challenges can be partially overcome using the EPA's power system data crosswalk, which identifies most of these relationships between plant IDs, generators, and units. When comparing data between the two datasets, we want to ensure that we are comparing the same sources of data. When units and generators have a 1:1 relationship, this is straightforward. However, when there are m:m relationships, this becomes more challenging. Thus, we create new subplant IDs, which identify distinct clusters of these unit-generator relationships.

Besides identifying these subplant clusters, we also want to identify if we have complete data for each cluster when we are comparing it, for example when calculating a gross to net generation ratio. For example, if we have a plant with units A, B, and C that are related to generators 1 and 2 in a m:m relationship, we want to ensure that we have data for all three units and both generators before comparing the data. For example, if we are missing data from unit C in a certain month, we would want to not use this month's data for calculating a gross to net ratio since the data is incomplete. However, sometimes certain units or generators in a subplant cluster retire over time. So if the reason that we were missing data from unit C in that month was because that unit retired in the previous month, we would now consider the remaining units and generators to be the complete subplant, so we should use that data. 

This process has now been integrated into the main data pipeline, using five years of historical data.

We have replaced `data_cleaning.identify_emissions_data_source()` with `data_cleaning.identify_hourly_data_source()` which now does this matching at the subplant level.

## Gross to net generation conversion
- Integrates the multi-year regressions of gross to net generation into the main data pipeline. 
- Improves the speed and memory use requirements of loading multiple years of data by aggregating the hourly data to monthly upon loading.
- If the regression leads to an intecept > 0 (implying house loads are adding generation), re-calculate the regression, forcing the intercept through zero.
- Implements a hierarchical approach to converting hourly gross generation in CEMS to net generation, which uses the following approaches, in order:
    - subplant ratio
    - subplant regression
    - plant ratio
    - plant regression

## Other
- Respond to validation checks that are failing.
- Clean up the `test_distribute_923.ipynb` notebook to remove cells that we no longer are using for testing. Keeping this notebook for now becuase it is still useful for loading and exploring EIA tables from PUDL
- Improve power sector data crosswalk with additional manual matches

-------------------------------------------------------------------------------
PR 2022-05-14
-------------------------------------------------------------------------------
### `calculate_residual_net_generation.ipynb`
- NOTE: This notebook modifies the code from `clean_930_compare_residual.ipynb` (I created a separate notebook for now because I did not want to delete your work), but this new notebook is what I've been integrating into the data pipeline.

Changes:
- Load chalendar data using `load_chalendar_for_pipeline()` instead of `load_chalendar()`
- Load cleaned CEMS net generation data data that was exported from `data_pipeline.ipynb` to `cems_2020_for_residual.csv`
- For each BA-fuel, aggregate the hourly CEMS data four different ways to test which aggregation works best:
    - All generation aggregated by commercial BA
    - All generation aggregated by physical BA
    - Transmission-connected generation aggregated by commercial BA
    - Transmission-connected generation aggregated by physical BA
- Calculates and visualizes residual based on raw CEMS data
- Whenever the CEMS netgen > 930 netgen, scale the entire timeseries so that in each hour the CEMS value is always <= the 930 netgen
- This scaling produces a residual that should always be positive
- Export the scaled residual for use in the data pipeline


### `data_pipeline.ipynb`
- Adds documentation and organizes pipeline into 12 steps
- Adds test and validation functions in `src.validation` to check cleaned CEMS and EIA-923 data for anomolous values (e.g. fuel consumption should always be positive, identify where there is still missing data, etc)
- Move more of the code into defined functions in `load_data` and `data_cleaning` 
- Create final output csvs for hourly BA data and plant-level data, and export to `data/final_outputs/`
- Previously, when converting UTC data to local time, we were using a UTC offset, meaning that we were converting to local standard time instead of local prevailing time. This meant when groupbing data by month, the data during DST would be off by an hour compared to the local prevailing month. Instead of using UTC offsets to convert to local time, we are now using tz database names (e.g. "US/Eastern"), which should automatically convert to prevailing time. At the BA level, we used the BA reference spreadsheet downloaded from the EIA-930 about page to create a new static table `data/manual/ba_reference.csv`

### Gross to Net Generation
- In `test_gross_to_net.ipynb` and `gross_to_net_generation.py`, investigate effectiveness of multi-year regression of monthly gross and net generation data
- We use `pudl.analysis.epa_crosswalk` to identify new "subplant_id"s for each plant. These subplant_id represent distinct clusters of EPA unitid and EIA generator_id, essentially the smallest subplant grouping at which it makes sense to aggregate data between CEMS and EIA data. The pudl module uses graph analysis to cluster units and generators based on their relationships (e.g. 1:1, 1:m, m:m) 
- We aggregate CEMS gross generation and EIA net generation by subplant and month, and run a simple linear regression. In theory, the slope should represent losses that scale with the amount of generation, and the intercept should represent static "house loads" that are used to power the plant, regardless of the amount of gross generation in any hour. For these regressions, we first convert total monthly MWh to average hourly MW for each month, based on the number of hours in each month. This helps control for differences in the number of days in each month, and it means that the regression coefficients are more directly applicable for converting hourly gross generation values.
- We export these results to `data/output/gross_to_net_regression.csv`
- These regressions have not yet been integrated into the main pipeline.
- The idea that the gross to net regressions would be run prior to running the main pipeline in order to 1) generate these regression values and 2) generate a list of subplant ids that can be used in the pipeline

### Other 
- Implement code formatting using `black`
- Add `.flake8` to ensure compatability with black formatter. See: https://github.com/psf/black/blob/06ccb88bf2bd35a4dc5d591bb296b5b299d07323/docs/guides/using_black_with_other_tools.md#flake8
- Add `ISSUES.md` to keep note of outstanding issues to address in the code. Eventually this should be moved to Github Issues
- Add `CHANGELOG.md` to keep note of changes made for each PR/commit


-------------------------------------------------------------------------------
Previous Work
-------------------------------------------------------------------------------
May 5
- [x] Adjust CEMS emissions for CHP (calculate fuel consumed for electricity?)
- [x] Update data_cleaning.get_epa_unit_fuel_types() to load manual data and look up fuel code of matched generator

May 4
 - [x] flag distribution-connected plants for 930 matching
 - [x] Assign physical BA for 930 matching
 - [x] Assign primary fuel code for 930 matching
 - [x] Explore why fraction > 1 in some cases
 - [x] Identified where eGRID data doesn't match EIA-923 data


April 29
 - [x] update fuel allocation
 - [x] add fuel for electricity
 - [x] test/update net generation allocation in pudl code
 - [x] remove non-grid connected plants
 - [x] add calculation of emissions
 - [x] remove non-grid connected plants
 - [x] remove retired or proposed plants unless they report data
 - [x] remove PR plants
 - [x] Integrate 930 shapes into profile
 - [x] assign plant primary fuel
 - [x] Point conda env to git+git://github.com/catalyst-cooperative/pudl@allocate_fuel (needs to be editable)
 - [x] Delete explore_data_cleaning.ipynb
 - [x] Test creating a env using the git branch

April 27
- [x] install pudl 
- [x] Integrate PUDL and heat allocation 

April 26
 - Added a new manual folder in the data directory for manually-created tables
 - Moved all of the data cleaning/investigation functions from `data_pipeline.ipynb` to `explore_data_cleaning.ipynb`
 - [x] Create hourly output timeseries showing data from CEMS and EIA
 - [x] Continue manually fixing prime mover codes
 - [x] Removed leading zeros from generator and unitids to help with matching

April 16, 2022
- [x] Figure out why output metrics don't have the same number of plants
- [x] Downloaded updated version of pudl database (v6.0.0)
- [x] Manually fix incorrect BA assignments
- [x] Fix BA assignments
- [x] Add plant count to BA metric

April 15, 2022
 - [x] move GTN functions to script
 - [x] in allocated EIA data, change "missing" co2 values to 0 for clean energy plants
 - [x] update geothermal emission factors


April 13, 2022
 - [x] Assign each plant to a balancing authority
 - [x] Aggregate all calculated generation and emissions data from CEMS and EIA to the annual level for each plant
 - [x] Developed functions to compare calculated totals to the reported eGRID totals
 - [x] Added nuclear generation to EIA data


April 12, 2022
 - [x] Fix rate issue in CEMS: the `co2_mass_tons` and `heat_content_mmbtu` columns are actually hourly rates, rather than measurements (see https://github.com/catalyst-cooperative/pudl/issues/1581). Thus, in `load_data.load_cems_data`, these columns are corrected by multiplying them by `operating_time_hours`
 - [x] Remove non-grid connected plants from EIA