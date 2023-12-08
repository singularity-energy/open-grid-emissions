---
stoplight-id: assign_esc
---

There are 41 unique `energy_source_code` for specific fuels that may be consumed by or associated with each boiler, generator, or plant (a complete list of the `energy_source_code` used by this project can be found [here](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/energy_source_groups.csv)).

## Assigning energy source codes to CEMS data
Although CEMS reports the total amount of fuel consumed in each hour, it does not include information about the type of fuel being consumed in each hour. Energy source codes are assigned to each unit in a three step process:
1. Use the energy source code assigned to the unit in the EPA's Power Sector Data Crosswalk
2. If the entire plant only consumed a single fuel in a month (according to the EIA-923 GF table) use that fuel type
3. Otherwise use the `plant_primary_fuel` determined by [this method](../Data%20Aggregation/Plant%20Primary%20Fuel.md).

## Updating "Other" (OTH) Energy Source Codes
Sometimes the `energy_source_code` reported for a specific data point is "Other" (OTH). A fuel type of other may be reported when the fuel consumed does not exactly fit with one of the 40 other fuel types used by the EIA. Other fuel codes are problematic because there are no emissions factors associated with them, so emissions associated with other fuel consumption will be other. Thus, the data pipeline manually updates all "other" fuel types with specific named fuel types.

All of the manual fuel type updates that we make can be found in [this table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/updated_oth_energy_source_codes.csv). We manually identified these updated codes by examining the heat content of the OTH fuel reported for each plant in the IEA-923 GF table and comparing it to the range of heat content of other named fuels. In certain cases, the plant name also provided some clues (e.g. if the generator was located at an oil refinery, or a paper mill, or a phosphate plant).

![fuel_heat_content.png](https://stoplight.io/api/v1/projects/cHJqOjE1MzAxNA/images/xD62XbYfh7Q)


## Energy Source Code for Municipal Solid Waste

Historically, any generators that burned municipal solid waste (MSW) reported this fuel consumption under a single fuel code. In recent years, however, EIA-923 began reporting these data under two separate codes for the biogenic portion (MSB) and non-biogenic portion (MSN). This is important because each portion has different emission rates. When calculating emissions from fuel consumption for EIA-923 data, these two separate fuel codes are used. However, when imputing missing emissions data in CEMS, the fuel type assigned to each unit is sourced from the Power Sector Data Crosswalk, which still only uses the MSW fuel code. This will be fixed in future versions of the database, and progress on this issue can be tracked [here](https://github.com/singularity-energy/open-grid-emissions/issues/51).

## Future Work, Known Issues, and Open Questions
- implement programmatic approach to updating OTH fuel codes ([details](https://github.com/singularity-energy/open-grid-emissions/issues/48))
- Use more specific energy source codes when available ([details](https://github.com/singularity-energy/open-grid-emissions/issues/195))
- Consider using startup-specific fuel codes to fill missing CEMS data during unit startup ([details](https://github.com/singularity-energy/open-grid-emissions/issues/155))
- Do we need to add digester gas (DG) fuel code and emissions factor to the data ([details](https://github.com/singularity-energy/open-grid-emissions/issues/72))