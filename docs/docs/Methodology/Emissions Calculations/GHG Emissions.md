---
stoplight-id: ghg_emissions
---

## Calculating GHG Emissions

The Open Grid Emissions Initiative includes emission mass data for the three primary greenhouse gases (GHGs) emitted by electric generators that contibute to human-driven climate change: carbon dioxide (CO2), methane (CH4), and nitrous oxide (N2O).

Measured CO2 data is reported in CEMS, but otherwise these emissions values must be estimated by multiplying the total heat input (mmBTU) of each fuel consumed by a [fuel-specific emission factor](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/emission_factors_for_co2_ch4_n2o.csv) (lb per mmBTU). Generally, these emissions factors for GHGs only depend on the type of fuel being consumed (i.e. they are matched with fuel consumption solely based on `energy_source_code`). The exception to this rule is geothermal emissions.

## Geothermal emissions
Geothermal plants do not combust fuel, although geothermal reservoirs do contain dissolved CO2. Binary geothermal power plants (prime mover is `BT`) operate as a closed loop system, so have no direct emissions to the atmosphere. Flash steam and dry steam power plants operate as an open loop, and thus some of this CO2 is emitted to the atmosphere. Dry steam geothermal plants have a higher CO2 emission factor than flash steam plants. Thus, the appropriate CO2 emission factor to use for geothermal plants depends not only on the fuel type (`GEO`), but also the prime mover code.

Although EIA data includes a prime mover code for binary turbines (`BT`), there is only a single code used for steam turbines (`ST`) so this alone cannot be used to assign the correct emission factor. However, based on information from the EIA and data from NREL, most steam geothermal plants in the U.S. are flash steam, and only plants located at the Geysers Complex in California operate as dry steam plants.

Some geothermal plants may use multiple prime mover types at a single plant. If we have generator-specific fuel consumption data, emissions are calculated on the basis of the prime mover used by each generator. If only plant-level fuel consumption data is available, then we calculate a weighted average emission factor based on the nameplate capacity of each generator at a plant, based on the prime mover associated with each generator.

All geothermal power plants are assumed to have no CH4 or N2O emissions.

## Calculating CO2-equivalent emissions
Although CO2, CH4, and N2O are all GHGs, they each contribute to global warming differently, as described by its Global Warming Potential (GWP). These GWPs can be used to calculate a CO2-equivalent (CO2e) value. This formula is:

Mass<sub>CO2e</sub> = (Mass<sub>co2</sub> * GWP<sub>co2</sub>) + (Mass<sub>ch4</sub> * GWP<sub>ch4</sub>) + (Mass<sub>n2o</sub> * GWP<sub>n2o</sub>)


The Intergovernmental Panel on Climate Change (IPCC) regularly updates these values over time in published Assessment Reports (AR). The values change over time both due to improvements in our scientific understanding of global warming, and because as the atmospheric concentration of GHGs changes over time, the GWP of each gas also changes. Thus, for each year, we use the most recently-published GWP values that were available in the data year. For example, AR5 was published in 2014, and AR6 was published in 2021. Thus, emissions occuring in years 2014-2020 will be converted using the AR5 values, and emissions occuring in 2021 and later will use the AR6 values.

In addition, the IPCC publishes GWPs on both a 20-year time horizon and a 100-year time horizon. The 100-year values are used as the default for all published CO2e values in the OGEI.

Finally, in AR5, the IPCC published GWPs that both included and excluded the impact of what they called the "climate-carbon feedback." As a default, all CO2e values published as part of the OGEI that use AR5 values use the GWPs with climate-carbon feedback. Starting with AR6, GWPs by default include climate carbon feedback.

## Future Work, Known Issues, and Open Questions
- Consider using fuel-weighted emissions factors for filling missing emissions data ([details](https://github.com/singularity-energy/open-grid-emissions/issues/163))