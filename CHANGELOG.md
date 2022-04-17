TODOs:
- [ ] Continue manually fixing prime mover codes
- [ ] Allocate net generation for electricity production, not steam
- [ ] For net generation, allocate difference between gen and genfuel tables
- [ ] Adjust data for CHP and biomass
- [ ] clean fuel codes based on EPA static tables
- [ ] find all of the missing plants




# Done

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