---
stoplight-id: input_data
---

## Input Data Sources

**EPA Continuous Emissions Monitoring System (CEMS) data**
* **What is it**: Measured hourly gross generation, fuel consumption, and emissions (CO2, NOx, SO2) data for emitting power generation units > 25MW
* **How we use it**: Primary source for hourly emissions and generation data
* **More information**: [https://ampd.epa.gov/ampd/](https://ampd.epa.gov/ampd/)

**EPA Mercury and Air Toxics Standards (MATS) data**
* **What is it**: Measured hourly mercury (Hg) emissions reported to the EPA by coal- and petroleum coke-fired generating units that continuously monitor Hg, from 2015 onward
* **How we use it**: Source of hourly Hg emissions data, loaded from PUDL. Hg is only reported for subplants and plants where the data is complete ([methodology](../Methodology/Emissions%20Calculations/Hg%20Emissions.md))
* **More information**: [https://www.epa.gov/stationary-sources-air-pollution/mercury-and-air-toxics-standards](https://www.epa.gov/stationary-sources-air-pollution/mercury-and-air-toxics-standards)

**EIA Form 923**
* **What is it**: Reported monthly net generation and fuel consumption data for power generators > 1 MW
* **How we use it**: To convert gross generation data from CEMS to net generation, and to fill in generation and emissions data for generator-months that are not reported in CEMS. EIA-923 reports Puerto Rico plants on a separate sheet, which PUDL does not yet include, so we load these data from the raw EIA-923 files.
* **More information**: [https://www.eia.gov/electricity/data/eia923/](https://www.eia.gov/electricity/data/eia923/)

**EIA Form 860**
* **What is it**: Inventory of all generators and plants and their static characteristics
* **How we use it**: to transform and aggregate the data reported in CEMS and EIA-923 based on plant and generator characteristics; calculating NOx and SO2 emissions based on plant-specific emissions control equipment and boiler types
* **More information**: [https://www.eia.gov/electricity/data/eia860/](https://www.eia.gov/electricity/data/eia860/)

**EIA Form 860M**
* **What is it**: Monthly update to the EIA-860 inventory of generators. It is the only EIA source of generator data for Puerto Rico, which is not included in the annual EIA-860.
* **How we use it**: For Puerto Rico generators, from 2018 onward. PUDL currently only includes the most recent EIA-860M data for Puerto Rico, so we use the December EIA-860M file of each year, as eGRID does. Plants that report generation to EIA-923 before they appear in EIA-860M use the generators in PUDL's most recent EIA-860M data instead.
* **More information**: [https://www.eia.gov/electricity/data/eia860m/](https://www.eia.gov/electricity/data/eia860m/)

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
* **How we use it**: To crosswalk data between CEMS and EIA-923, and to identify the CEMS units that the EPA excluded from the crosswalk because they are not connected to the grid
* **More information**: [https://www.epa.gov/airmarkets/power-sector-data-crosswalk](https://www.epa.gov/airmarkets/power-sector-data-crosswalk)

**EPA Clean Air Markets Program Data (CAMPD) facility attributes**
* **What is it**: Annual list of the units that report to the EPA's Clean Air Markets programs, with attributes such as their unit type, fuel types, emissions controls, and associated generators
* **How we use it**: To fill boiler firing types and bottom types that are missing from EIA-860 when calculating NOx and SO2 emissions, as a backstop for the fuel type of CEMS units, and to identify updates to the crosswalk between CEMS units and EIA generators
* **More information**: [https://campd.epa.gov/data/bulk-data-files](https://campd.epa.gov/data/bulk-data-files)

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
* **How we use it**: Our cleaning of raw EIA-930 data (`oge.physics_reconciliation`) and our calculation of consumption-based emissions totals are adapted from this package
* **More information**: [https://gridemissions.jdechalendar.su.domains/#/code](https://gridemissions.jdechalendar.su.domains/#/code)