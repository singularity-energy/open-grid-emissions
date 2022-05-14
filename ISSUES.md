
## General data cleaning
- [ ] Update `identify_emissions_data_source` to match on subplant - ensure that all active units are reporting
- [ ] Convert UTC to local timezone using TZ name rather than GMT offset to account for DST
- [ ] Move all manual update dictionaries to csv tables
- [ ] Instead of creating manual crosswalk, fork the PDSC repo and make changes that can be committed
- [ ] Figure out what to do with retired BAs
- [ ] Figure out what to do with energy storage
- [ ] Remove steam-only plants from CEMS (or plants that have no generation)

## Emissions data
- [ ] ID plants missing from egrid that have near-zero data
- [ ] Finish identifying cems_ids with missing fuel codes
- [ ] Update fuel type assignment in CEMS to use monthly and/or weighed average EFs for filling
- [ ] clean fuel codes based on EPA static tables
- [ ] Implement biomass adjustments

## EIA-923 Data Cleaning
- [ ] Allocate 99999 data from EIA 923
- [ ] Check performance of allocation when net generation is negative

## EIA-930 data
- [ ] Need to check that all timestamps in cleaned chalendar data have been converted to hour starting timestamps (instead of hour ending) to ensure compatibility with hourly CEMS data (which uses hour beginning timestamps)
- [ ] In `gridemissions/src/gridemissions/clean.py`, lines 92-100, chalendar currently assigns a static flat profile to geothermal and biomass emissions. These lines should be removed for now. We can add flat profiles for these fuel types in the main data pipeline.

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

## Documentation
- [ ] Write up gross to net generation calculation / explanation
- [ ] Draft document proposing hierarchy for residual profiles and risks for each fuel type
- [ ] Draft documentation of updated fuel/net gen allocation methodology


## Gross to net generation conversion
### Short-run
Multi-year regression
- [ ] Increase to 5 years of data
- [ ] When dividing by number of units, ensure that all are reporting data
- [ ] Figure out how to deal with units that retire over time
- [ ] Look into combined cycle units (GTN > 1)
Data pipeline
- [ ] Implement hierarchy of methods for converting gross to net
- [ ] Look into whether imputation of missing gross generation is needed
   - very small amount of emissions
   - if filling, need to match on unit, rather than plant (57075 - Ivanpah)
- [ ] Fix imputation of missing hourly gross generation in CEMS

Validation
- [ ] Compare totals to actual values at monthly and annual level
- [ ] Compare methods for calculating net generation

### Long-run
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





