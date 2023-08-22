---
stoplight-id: pipeline_overview
---

# Overview of the data pipeline

1. **Download data**: Download data, including CEMS (via PUDL), EIA Forms 860, 861, and 923 (raw data and via PUDL), EPA-EIA Power Sector Data Crosswalk, EIA-930 data
2. **Identify subplants**: Identify distinct “subplants” which represent interconnected groupings of units, generators, and boilers within a plant using graph analysis, based on relationships reported in the power sector data crosswalk and the EIA-860 boiler-generator association table. ([methodology](../Methodology/Data%20Aggregation/Subplant%20Aggregation.md))
3. **Clean EIA-923 data**: Clean and allocate monthly generation and fuel data from EIA-923 ([methodology](../Methodology/Data%20Cleaning/EIA-923%20Data.md))
    1. Allocate monthly net generation and fuel consumption data reported at the prime-mover level in the EIA-923 generation and fuel table to individual generators. Data is allocated to each generator-month in proportion to data in the generator table and boiler fuel table, if reported, otherwise it is allocated in proportion to the nameplate capacity of each generator.
    2. Assign an energy source code to any “other” fuel consumed based on the heat content of the “other” fuel and likely fuels based on the plant type ([methodology](../Methodology/Emissions%20Calculations/Assigning%20Energy%20Source%20Codes.md))
    3. Identify an annual primary fuel for each generator and plant based on fuels consumed. ([methodology](../Methodology/Data%20Aggregation/Plant%20Primary%20Fuel.md))
    4. Calculate monthly CO2, CH4, and N2O emissions for each generator based on its fuel consumption of each fuel type ([methodology](../Methodology/Emissions%20Calculations/GHG%20Emissions.md))
    5. Calculate monthly NOx ([methodology](../Methodology/Emissions%20Calculations/NOx%20Emissions.md)) and SO2 ([methodology](../Methodology/Emissions%20Calculations/SO2%20Emissions.md)) emissions from each generator based on its fuel consumption of each fuel type, as well as the boiler firing type and any emissions controls
    6. Calculate biomass-adjusted emissions ([methodology](../Methodology/Emissions%20Calculations/Adjusting%20Emissions%20for%20Biomass.md)), CHP-adjusted emissions ([methodology](../Methodology/Emissions%20Calculations/Adjusting%20Emissions%20for%20CHP.md)), and biomass- and CHP-adjusted emissions
    7. Calculate CO2e emissions using 100-year GWPs from the appropriate IPCC Assessment Report. ([methodology](../Methodology/Emissions%20Calculations/GHG%20Emissions.md))
    8. Aggregate allocated data by generator
    9. Remove non-grid connected plants
4. **Clean CEMS data**: Clean hourly generation, fuel, and emissions data from CEMS ([methodology](../Methodology/Data%20Cleaning/CEMS%20Data.md))
    1. Load data from pudl and crosswalk CAMD plant IDs to EIA plant IDs
    2. Remove data for non grid-connected plants, plants in Puerto Rico, and certain units that only report steam load and do not report to EIA
    3. Assign a monthly "report_date" to each hourly observation based on the date of the local timestamp (this allows us to match the data to EIA-923 monthly report dates)
    4. Remove data for any unit-months where there is incomplete hourly data.
    5. Assign a fuel type to each unit based on the power sector data crosswalk
    6. Fill in missing hourly emissions data using the assigned fuel type and reported hourly fuel consumption data
    7. Remove all observations for each unit-month when no operation is reported for that unit in that month (allows us to fill this data using EIA-923 data if available).
    8. Calculate biomass-adjusted emissions
5. **Assign static attributes to all plants** This includes primary fuel, data sources, state, balancing authority (both commercial and physical), whether the plant is connected to the distribution grid or transmission grid, a fuel category, and local timezone
6. **Crosswalk the EIA-923 and CEMS data**, identifying for each subplant-month whether there is complete CEMS data available, partial CEMS data available (available for some but not all units that consist a subplant), or only EIA data available.
    1. For subplant-months with complete CEMS data, we will use the hourly CEMS data directly
    2. For subplant-months with partial CEMS data, we will use the partial CEMS data to assign an hourly profile to the monthly EIA-923 data for that subplant ([methodology](../Methodology/Assigning%20Hourly%20Profiles%20to%20Monthly%20Data/Scaling%20Partial%20CEMS%20Subplant%20Data.md))
    3. For plant-months where certain subplants have CEMS data, we will use the CEMS data from those subplants to shape the EIA-923 data for the subplants for which no CEMS data is available ([methodology](../Methodology/Assigning%20Hourly%20Profiles%20to%20Monthly%20Data/Shaping%20Partial%20CEMS%20Plant%20Data.md))
    4. For plant-months where no CEMS data is available, we impute the hourly profile based on reported EIA-930 net generation data
7. **Aggregate CEMS data to subplant**: Aggregate hourly CEMS data from units to subplants
8. **Shape partial CEMS data**: For subplants with partial CEMS data, use the available CEMS data to assign an hourly profile to the EIA-923 data reported for that subplant-month.
9. **Convert CEMS hourly gross generation to hourly net generation** ([methodology](../Methodology/Converting%20Gross%20to%20Net%20Generation.md))
    1. Calculate various gross-to-net conversion factors based on gross generation reported in CEMS and net generation reported in EIA, including subplant- and plant-specific, monthly and annual gross to net ratios, shift factors, and linear regressions
    2. Filter out any GTN conversion factors that would result in unreasonable net generation numbers
    3. Apply GTN conversion factors hierarchically, using the best available factor for each subplant
10. **Adjust CEMS emissions**: Calculate hourly CHP-adjusted emissions for the CEMS data ([methodology](../Methodology/Emissions%20Calculations/Adjusting%20Emissions%20for%20CHP.md)), as well as CO2e emissions ([methodology](../Methodology/Emissions%20Calculations/GHG%20Emissions.md))
11. **Export monthly and annual plant data:** Export plant-level results at the monthly and annual aggregations
12. **Clean the EIA-930 data** using a physics-based reconciliation algorithm that ensure that reported data respects conservation of energy laws ([methodology](../Methodology/Data%20Cleaning/EIA-930%20Data.md))
13. **Calculate residual net generation profiles**: Estimate hourly profiles for monthly EIA data that is not reported to CEMS using EIA-930 data ([methodology](../Methodology/Assigning%20Hourly%20Profiles%20to%20Monthly%20Data/Shaping%20Using%20Fleet-Specific%20Profiles.md))
    1. For BAs where wind and solar data is missing, impute profiles based on the wind and solar profiles in neighboring BAs that are in the same time zone. For all other BA-fuels that are missing, assume a flat profile
    2. Calculate the total net generation for each BA-fuel reported in CEMS
    3. Calculate a residual net generation profile by subtracting CEMS net generation from EIA-930 net generation for each BA fuel. This should theoretically represent the hourly net generation profile of all generators that do not report to CEMS.
14. **Shape monthly EIA-923 data**: Assign the best available hourly profile for each BA-fuel to the monthly EIA-923 data ([methodology](../Methodology/Assigning%20Hourly%20Profiles%20to%20Monthly%20Data/Shaping%20Using%20Fleet-Specific%20Profiles.md))
15. **Combine hourly plant-level data** from all sources (CEMS, partial CEMS, and shaped EIA-923) and export hourly plant-level results (currently, shaped EIA-923 data is aggregated to the BA-fuel level instead of individual plants to prevent false precision and keep data size reasonable)
16. **Export power sector (BA-fuel level) data**
17. **Calculate and export consumed emission factors** based on hourly interchange between BAs reported in EIA-930 ([methodology](../Methodology/Emissions%20Calculations/Consumption-based%20Emissions.md))
