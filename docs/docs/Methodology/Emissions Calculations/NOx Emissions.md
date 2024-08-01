---
stoplight-id: nox_emissions
---

## Calculating NOx emissions
Nitrogen Oxides (NOx) are air pollutants emitted by electric generators as the result of combusting fuels, and according to the U.S. EPA are precursors to the formation of ozone (smog), fine particulates, acid rain, and other environmental and human health impacts.

Measured NOx data is reported in CEMS, but otherwise must be estimated. Unlike greenhouse gases, whose emissions depend solely on the heat input and type of fuel consumed, NOx emissions additionally depend on the prime mover type, several boiler design parameters, and any NOx control equipment that is installed in the stack. The process of estimating NOx emissions thus consists of three primary steps:
1. Assign a fuel-, prime mover-, boiler firing type-, and boiler bottom type-specific NOx emission factor to each observation.
2. Calculate an uncontrolled emissions mass by multiplying this emission factor by `fuel_consumed_mmbtu`
3. Calculate controlled NOx emissions using reported controlled emissions factors if available.

## Data inputs for NOx emissions calculations
Calculating uncontrolled and controlled NOx emissions relies on boiler design, environmental attribute, and environmental control data that is reported in EIA-860 and EIA-923, but which are contained in tables that are not currently included in the PUDL data pipeline. Until these tables are integrated into PUDL (track progress [here](https://github.com/singularity-energy/open-grid-emissions/issues/154)), the OGE data pipeline downloads the raw EIA-860 and EIA-923 webfiles and loads these tables directly from those files.

## Detailed steps for calculating NOx emissions

1.  Load boiler design-, prime mover-, and fuel-specific [uncontrolled NOx factors](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/oge/reference_tables/emission_factors_for_nox.csv). These factors are derived from [Table A-2](https://www.eia.gov/electricity/annual/html/epa_a_02.html) of the U.S. EIA's Electric Power Annual.
2. Assign these uncontrolled NOx emission factors to each boiler based on the relevant parameters. If an emission factor is not available for a specific boiler design, we default to using the value for a dry-bottom, "other" firing type boiler, which is consistent with the assumption made by the EIA's Electric Power Annual (https://www.eia.gov/electricity/annual/pdf/tech_notes.pdf)
3. Because the fuel consumption data is reported at the generator-fuel level, we crosswalk each boiler-fuel to each generator-fuel and then calculate calculate a simple average NOx emission factor for all boilers associated with each generator.
4. Because emission factors are reported in lb NOx per unit mass or unit volume of fuel consumed (instead of per mmbtu), we then convert these factors to lb per mmbtu. In the EIA-923 boiler fuel table, each boiler reports the month-specific heat content (mmbtu) per mass of fuel consumed for each boiler. If there is not a boiler-specific heat content value available, we use the month-specific national average heat content of that fuel. If a month-specific national average heat content value is not available, we use the national and annual average heat content of that fuel.
5.  After calculating uncontrolled NOx emissions for each generator, we then calculated controlled NOx emissions. Actual generator-specific and control equipment specific controlled emission rates are reported by certain generators, for both the ozone season (May-Sept) and for the annual average, in EIA-923 Schedule 8. We calculate a non-ozone season average rate using the annual average and the ozone season average.
6.  Certain generators use multiple types of control equipment, but the total number of hours per year that each control equipment was operating is reported. Thus, when aggregating the controlled emission factor up to the generator level, we are able to calculate a weighted average emission factor using the number of operating hours as the weighting factor.
7. We then calculate controlled NOx emissions for each month using the appropriate factor (for ozone or non-ozone season). If there are not season-specific factors available, we use the annual average factor.
8. If there is a controlled emission total available, it is used as the final NOx value. Otherwise, we assume emissions are uncontrolled, and we use the calculated uncontrolled NOx value.


## NOx emissions from Geothermal Generators
According to the EPA's eGRID documentation, certain types of geothermal plants can have NOx emissions (see [this issue](https://github.com/singularity-energy/open-grid-emissions/issues/69)):
> While CO2 is a gas in the geothermal reservoir, SO2 and NOx result from hydrogen sulfide combustion.

These factors are loaded separately from our [manual geothermal emissions factor table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/geothermal_emission_factors.csv). These factors are added before calculating uncontrolled NOx emissions.

## Future Work, Known Issues, and Open Questions
- Missing NOx emissions in CEMS are not currently being imputed ([details](https://github.com/singularity-energy/open-grid-emissions/issues/153))
- Need to determine standard fallback behavior when boiler design unknown ([details](https://github.com/singularity-energy/open-grid-emissions/issues/150))
- We are uncertain whether Nox data reported in CEMS already accounts for emissions controls ([details](https://github.com/singularity-energy/open-grid-emissions/issues/151))
- Data for some plants that only report NOx data to CEMS is currently being dropped in favor of EIA-923 data ([details](https://github.com/singularity-energy/open-grid-emissions/issues/134))