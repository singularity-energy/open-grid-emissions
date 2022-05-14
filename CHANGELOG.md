
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