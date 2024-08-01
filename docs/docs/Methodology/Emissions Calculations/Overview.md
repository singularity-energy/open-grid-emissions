---
stoplight-id: emissions_overview
---

## Overview

The Open Grid Emissions Initiative includes power sector emissions data for the three primary greenhouse gases (GHGs), CO2, CH4, and N2O, as well as other air pollutants including NOx and SO2. These data represent only direct power sector emissions from the combustion or consumption of fossil, biogenic, and other fuels. This dataset does not include any lifecycle emissions related to the extraction or transportation of these fuels, or the construction or operation of the power plant. Emissions are reported both on a mass basis, and on an intensity basis (in the form of an emission factor that represents emissions mass normalized by generation of consumption of electricity). In addition to reporting total emissions, the OGE also includes emissions data that has been adjusted for biomass emissions and for combined heat and power (CHP) plants.

The published data reflects a mix of measured and estimated emissions data. The CEMS dataset reports measured emissions mass values for CO2, NOx, and SO2. All CH4 and N2O emissions, as well as CO2, NOx, and SO2 emissions for plants that do not report to CEMS are estimated based on the reported heat input of each fuel consumed and fuel-specific emission factors.

> TODO: insert table with overview of emissions data sources

## Emissions Factors
In addition to providing total emissions mass, OGE also provides emissions factors that describe the emissions intensity of the electricity being generated or consumed. These emissions factors are calcuated by dividing total emissions mass (lb) by total electric generation (MWh).

There are two different types of emissions factors: generated and consumed.

Generated emission factors describe the emissions intensity of generated electricity injected into the grid, and are calculated as

emission (in lb) &frasl; net generation (in MWh)

Consumed emission factors describe the emissions intensity of electricity consumed in a region, which is a mix of the electricity generated in a region and the electricity imported into the region. More details about this methodology can be found [here](../Emissions%20Calculations/Consumption-based%20Emissions.md).

### A note on emissions factors in small balancing authorities during power plant startup
In small or generation-only balancing authorities (with a small number of generators), the operational characteristics of a single plant can have a large impact on the regional emission factor. This is especially relevant during power plant startup, during which time a plant starts consuming fuel but may not be generating a large amount of electricity, thus leading to an emission factor that is orders of magnitude higher than the emission factor during normal operation. This is typically only an issue when observing plant-level hourly data, but this can also affect the regional emission factor for some small BAs.

## Future Work, Known Issues, and Open Questions
- How should we account for energy storage in emission factors ([details](https://github.com/singularity-energy/open-grid-emissions/issues/60))
- How should we account for plant startup in our emissions calcuations ([details](https://github.com/singularity-energy/open-grid-emissions/issues/155)