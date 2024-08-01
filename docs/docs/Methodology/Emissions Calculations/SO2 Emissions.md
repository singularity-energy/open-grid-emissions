---
stoplight-id: so2_emissions
---

## Calculating SO2 emissions
Sulfur Dioxide (SO2) is an air pollutant emitted by electric generators as the result of combusting fuels (especially coal), and according to the U.S. EPA is a precursor to the formation fine particulates, acid rain, and other environmental and human health impacts.

Measured SO2 data is reported in CEMS, but otherwise must be estimated. Unlike greenhouse gases, whose emissions depend solely on the heat input and type of fuel consumed, NOx emissions additionally depend on the sulfur content of the fuel, the prime mover type, the boiler firing type, and any SO2 control equipment that is installed in the stack. The process of estimating SO2 emissions thus consists of three primary steps:
1. Assign month- and generator-specific uncontrolled SO2 emission factors based on boiler design parameters of all associated boilers and sulfur content of consumed fuels.
2. Calculate an uncontrolled emissions mass by multiplying this emission factor by `fuel_consumed_mmbtu`
2. Adjust SO2 emissions using boiler-specific SO2 removal efficiencies if available. Otherwise, assume SO2 emissions are uncontrolled.

### Detailed Approach
1.  Load boiler design-, prime mover-, and fuel-specific [uncontrolled SO2 factors](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/oge/reference_tables/emission_factors_for_so2.csv). These factors are derived from [Table A-1](https://www.eia.gov/electricity/annual/html/epa_a_01.html) of the U.S. EIA's Electric Power Annual.
2. Assign these uncontrolled SO2 emission factors to each boiler based on the relevant parameters. If an emission factor is not available for a specific boiler design, we default to using the value for an "other" firing type boiler, which is consistent with the assumption made by the EIA's Electric Power Annual (https://www.eia.gov/electricity/annual/pdf/tech_notes.pdf)
3. Some SO2 emission factors depend on the specific sulfur content of the fuel burned by each boiler. The sulfur content of the fuel burned in each boiler-month is reported in EIA-923, so we use that to calculate the boiler-month-fuel specific SO2 emission factor. If boiler-month specific sulfur content data is reported in EIA-923, we use that value. If those values are missing, we utilitize backstop values in the following order, if available: month- and fuel-specific state-average values, year- and fuel-specific national-average values, fuel-specific national-average values from the previous year. Note that fuel sulfur content is not available for years 2005-2008 and 2008-2012 data is used to derive a plant, state and national fuel sulfur content average.
4. Because the fuel consumption data is reported at the generator-fuel level, we crosswalk each boiler-fuel to each generator-fuel and then calculate calculate a simple average SO2 emission factor for all boilers associated with each generator.
5. Because emission factors are reported in lb SO2 per unit mass or unit volume of fuel consumed (instead of per mmbtu), we then convert these factors to lb per mmbtu. In the EIA-923 boiler fuel table, each boiler reports the month-specific heat content (mmbtu) per mass of fuel consumed for each boiler. If there is not a boiler-specific heat content value available, we use the month-specific national average heat content of that fuel. If a month-specific national average heat content value is not available, we use the national and annual average heat content of that fuel.
6.  After calculating uncontrolled SO2 emissions for each generator, we then calculate controlled SO2 emissions. Certain boilers report the SO2 removal efficiency of all of the the SO2 control equipment that they have installed in EIA-923. Certain generators use multiple types of control equipment, but the total number of hours per year that each control equipment was operating is reported. Thus, when aggregating the SO2 removal efficiency up to the generator level, we are able to calculate a weighted average emission factor using the number of operating hours as the weighting factor. The exception is boilers with a fluidized bed boiler type, whose emission factors already represent controlled emissions in EIA (see: https://www.eia.gov/electricity/annual/pdf/tech_notes.pdf)
7. If there is an SO2 removal efficiency value available, final so2 emissions values are calculated by multiplying the total uncontrolled so2 emissions by (1 - so2 removal efficiency) for each generators. Otherwise, we assume SO2 emissions are uncontrolled, and we use the calculated uncontrolled SO2 value.


## SO2 emissions from geothermal generators
According to the EPA's eGRID documentation, certain types of geothermal plants can have SO2 emissions (see [this issue](https://github.com/singularity-energy/open-grid-emissions/issues/69)):
> While CO2 is a gas in the geothermal reservoir, SO2 and NOx result from hydrogen sulfide combustion.

These factors are loaded separately from our [manual geothermal emissions factor table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/geothermal_emission_factors.csv). These factors are added before calculating uncontrolled SO2 emissions.

## Future Work, Known Issues, and Open Questions
- Issues with SO2 emission factors for coal plants with fluidized bed boilers ([details](https://github.com/singularity-energy/open-grid-emissions/issues/248))
- Issues with SO2 emission factors for landfill gas generators ([details](https://github.com/singularity-energy/open-grid-emissions/issues/218))
- Missing SO2 emissions in CEMS are not currently being imputed ([details](https://github.com/singularity-energy/open-grid-emissions/issues/153))