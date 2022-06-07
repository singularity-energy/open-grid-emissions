# HIGH PRIORITY

## Outputs
- [ ] Separate data by hourly data source
- [-] Output final files

## Emissions
- [ ] replace all OTH fuel type
- [ ] Add functions to clean_cems for filling missing NOx and SOx
- [ ] Use plant-season specific Nox emissions factors from EIA-923 schedule 8c
- [ ] Calculate Nox emissions using factors specific to boiler firing type
- [ ] Adjust SO2 emissions for control efficiencies reported in EIA-923
- see: https://github.com/catalyst-cooperative/pudl/issues/889

## General data cleaning
- [-] When subplant ID is NA in CEMS, make sure it is correctly merging with the EIA data
- [ ] Check MSW codes - make sure this was done correctly
- [ ] Identify and remove hourly values in CEMS that appear to be outliers

## Gross to net generation conversion
- [ ] Test gross to net emissions ordering using residual 
- [ ] Create annual residual metric to evaluate what order I should do methods in
- [ ] Add method to deal with negative GTN ratios
- [ ] Revise assumption for assumed GTN value
- [-] Fix issue where net generation > gross generation

## EIA-923 Data Cleaning
- [ ] Allocate 99999 data from EIA 923
- [ ] Check performance of allocation when net generation is negative


# LOWER PRIORITY

## General data cleaning
- [ ] For heat rate validation test, group by fuel and PM
- [ ] Figure out what to do with retired BAs
- [ ] Figure out what to do with energy storage
- [ ] Remove steam-only plants from CEMS (or plants that have no generation)
- [ ] Move all manual update dictionaries to csv tables
- [ ] Instead of creating manual crosswalk, fork the PDSC repo and make changes that can be committed

## Gross to net generation conversion
- [ ] Filter out regression values if slope > 1?
- [ ] if regression intercept is positive, re-run regression forcing intercept through zero
- [ ] Improve method for assumed values in GTN conversion (identigy similar plant types)
- [ ] Look into whether imputation of missing gross generation is needed
   - very small amount of emissions
   - if filling, need to match on unit, rather than plant (57075 - Ivanpah)
- [ ] Look into combined cycle units (GTN > 1) to see if this is working correctly

## Emissions data
- [ ] Finish identifying cems_ids with missing fuel codes
- [ ] Update fuel type assignment in CEMS to use monthly and/or weighed average EFs for filling
- [-] clean "other" fuel codes based on EPA static tables
- [ ] ID plants missing from egrid that have near-zero data

## Validation
- [ ] Filter known issues from egrid validation metrics
- [ ] Ensure proper BA sorting
- [ ] Put validation functions inside functions
- [ ] Add prime mover code to EIA-923 allocation output so that we can validate incorrect PM

## Documentation
- [-] Write up gross to net generation calculation / explanation
- [-] Draft document proposing hierarchy for residual profiles and risks for each fuel type
- [ ] Draft documentation of updated fuel/net gen allocation methodology

## Gross to Net Regression (long-run issues / low priority)
Look into other regressors
- [ ] Investigate annual fixed effects
- [ ] Determine if/when major change in equipment (repower, new environmental controls)
- [ ] Investigate monthly fixed effects
- [ ] Capacity factor
- [ ] Binary variable for months where unit is operating and months when it is not
Other analysis
- [ ] If no net gen data for plant, apply standard efficiency for similar plants
   - do regression on annual values by fuel/PM/Envirocontrol
   - will there even be data for these plants if no net gen data?





