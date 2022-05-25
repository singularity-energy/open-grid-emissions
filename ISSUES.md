## General data cleaning
- [ ] fix issue where subplant id is NA - maybe assign 0?
- [ ] Figure out what to do with retired BAs
- [ ] Figure out what to do with energy storage
- [ ] Remove steam-only plants from CEMS (or plants that have no generation)
- [ ] Move all manual update dictionaries to csv tables
- [ ] Instead of creating manual crosswalk, fork the PDSC repo and make changes that can be committed

## Gross to net generation conversion
- [ ] Filter out regression values if slope > 1?
- [ ] if regression intercept is positive, re-run regression forcing intercept through zero
- [ ] Add method to deal with negative GTN ratios
- [ ] Improve method for assumed values in GTN conversion (identigy similar plant types)
- [ ] Look into whether imputation of missing gross generation is needed
   - very small amount of emissions
   - if filling, need to match on unit, rather than plant (57075 - Ivanpah)
- [ ] Look into combined cycle units (GTN > 1) to see if this is working correctly

## Emissions data
- [ ] Implement biomass adjustments
- [ ] Finish identifying cems_ids with missing fuel codes
- [ ] Update fuel type assignment in CEMS to use monthly and/or weighed average EFs for filling
- [-] clean "other" fuel codes based on EPA static tables
- [ ] ID plants missing from egrid that have near-zero data

## EIA-923 Data Cleaning
- [ ] Allocate 99999 data from EIA 923
- [ ] Check performance of allocation when net generation is negative

## EIA-930 data
- [-] In `gridemissions/src/gridemissions/clean.py`, lines 92-100, chalendar currently assigns a static flat profile to geothermal and biomass emissions. These lines should be removed for now. We can add flat profiles for these fuel types in the main data pipeline.

## Distributing EIA-923 to Hourly
- [ ] Distribute data at plant level instead of BA-fuel level
- [ ] Assign flat profile to Geothermal and biomass generation
- [ ] Identify BAs where "Other" category includes a single fuel type
- [ ] Assign a profile to missing wind and solar data if none provided in 930
For plants with negative net generation
 - [ ] If no fuel consumption, assign a flat profile
 - [ ] If fuel consumption, adjust hourly profile

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





