---
stoplight-id: input_data
---

## Input Data Sources

**EPA Continuous Emissions Monitoring System (CEMS) data**
* **What is it**: Measured hourly gross generation, fuel consumption, and emissions (CO2, NOx, SO2) data for emitting power generation units > 25MW
* **How we use it**: Primary source for hourly emissions and generation data
* **More information**: [https://ampd.epa.gov/ampd/](https://ampd.epa.gov/ampd/)

**EIA Form 923**
* **What is it**: Reported monthly net generation and fuel consumption data for power generators > 1 MW
* **How we use it**: To convert gross generation data from CEMS to net generation, and to fill in generation and emissions data for generator-months that are not reported in CEMS
* **More information**: [https://www.eia.gov/electricity/data/eia923/](https://www.eia.gov/electricity/data/eia923/)

**EIA Form 860**
* **What is it**: Inventory of all generators and plants and their static characteristics
* **How we use it**: to transform and aggregate the data reported in CEMS and EIA-923 based on plant and generator characteristics; calculating NOx and SO2 emissions based on plant-specific emissions control equipment and boiler types
* **More information**: [https://www.eia.gov/electricity/data/eia860/](https://www.eia.gov/electricity/data/eia860/)

**EIA Form 861**
* **What is it**: Data about electric utility and power marketer sales and operations
* **How we use it**: Calculating grid gross loss based on utility disposition
* **More information**: [https://www.eia.gov/electricity/data/eia861/](https://www.eia.gov/electricity/data/eia861/)

**Public Utility Data Liberation (PUDL) Project**
* **What is it**: an open source data processing pipeline that makes U.S. energy data easier to access and use programmatically
* **How we use it**: To access pre-cleaned and standardized versions of CEMS, EIA-923, and EIA-860 data, as well as the EPA-EIA Power Sector Data Crosswalk.
* **More information**: [https://catalystcoop-pudl.readthedocs.io/en/latest/](https://catalystcoop-pudl.readthedocs.io/en/latest/)

**EPA-EIA Power Sector Data Crosswalk**
* **What is it**: Maps EPA plant IDs and unit IDs to EIA plant IDs and generator IDs
* **How we use it**: To crosswalk data between CEMS and EIA-923
* **More information**: [https://www.epa.gov/airmarkets/power-sector-data-crosswalk](https://www.epa.gov/airmarkets/power-sector-data-crosswalk)

**EIA Form 930 / Hourly Electric Grid Monitor**
* **What is it**: Reported hourly net generation by fuel category, demand, and interchange for each Balancing Area in the U.S.
* **How we use it**: To assign an hourly profile to the monthly generation and fuel data reported in EIA-923 and to calculate hourly consumed emission factors
* **More information**: [https://www.eia.gov/electricity/gridmonitor/about](https://www.eia.gov/electricity/gridmonitor/about) and [https://www.eia.gov/survey/#eia-930](https://www.eia.gov/survey/#eia-930)

**EPA eGRID database**
* **What is it**: Reports annual-level generation and emissions statistics at the plant and BA level
* **How we use it**: to validate our outputs, and as a source for certain static tables (emission factors, conversion factors, etc)
* **More information**: [https://www.epa.gov/egrid](https://www.epa.gov/egrid)

**gridemissions repository**
* **What is it**: Tools for tracking EIA-930 power sector data
* **How we use it**: To clean raw EIA-930 data and to calculate consumption-based emissions totals
* **More information**: [https://gridemissions.jdechalendar.su.domains/#/code](https://gridemissions.jdechalendar.su.domains/#/code)