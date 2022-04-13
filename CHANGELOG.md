--------------------------------------------------
April 13, 2022
--------------------------------------------------
 - Assign each plant to a balancing authority
 - Aggregate all calculated generation and emissions data from CEMS and EIA to the annual level for each plant
 - Developed functions to compare calculated totals to the reported eGRID totals
 - Added nuclear generation to EIA data

Next steps:
 - Calculate gross to net generation ratio based on updated net generation calculations
 - update fuel types and locations based on static tables and SQL
 - Investigate missing generators
 - Emission adjustments for CHP and biomass

--------------------------------------------------
April 12, 2022
--------------------------------------------------
### CEMS
 - the `co2_mass_tons` and `heat_content_mmbtu` columns are actually hourly rates, rather than measurements (see https://github.com/catalyst-cooperative/pudl/issues/1581). Thus, in `load_data.load_cems_data`, these columns are corrected by multiplying them by `operating_time_hours`

### EIA
 - Remove non-grid connected plants