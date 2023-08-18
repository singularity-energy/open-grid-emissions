---
stoplight-id: cleaning_eia923
---


## Loading data from PUDL
All of the generation and fuel consumption data that the OGEI uses from EIA-923 comes from pre-cleaned tables from the PUDL project, rather than from raw EIA-923 webfiles. The PUDL transformation process cleans the input data so that it is adjusted for uniformity, corrected for errors, and ready for bulk programmatic use. More information about the transformations that PUDL applies to the raw data can be found [here](https://catalystcoop-pudl.readthedocs.io/en/latest/data_sources/eia923.html#pudl-data-transformations).

Because annual environmental information from EIA-923 has not yet been integrated into PUDL, these tables are loaded from the raw files, and column names are changed so that they are consistent with the column names used in PUDL.

## Allocating Generation and Fuel Data

The EIA-923 dataset reports net generation and fuel consumption data in three different tables, as summarized in the below table.

<table>
  <tr>
   <td>Table Name / abbreviation
   </td>
   <td>Aggregation level
   </td>
   <td>Required reporters
   </td>
   <td>Data reported
   </td>
  </tr>
  <tr>
   <td>Generation and fuel (GF)
   </td>
   <td>Prime-mover / fuel
   </td>
   <td>Any grid-connected plant >1MW nameplate capacity
   </td>
   <td>Monthly net generation and fuel consumption (separated by total and consumed for electricity generation)
   </td>
  </tr>
  <tr>
   <td>Boiler Fuel (BF)
   </td>
   <td>Boiler ID / fuel
   </td>
   <td rowspan="2" >Only units with steam turbines (ST) driven by steam from a boiler fueled by combustible fuels. Mostly consists of units > 10MW.
   </td>
   <td>Monthly fuel consumption totals
   </td>
  </tr>
  <tr>
   <td>Generator (G)
   </td>
   <td>Generator ID
   </td>
   <td>Monthly net generation totals
   </td>
  </tr>
</table>




Because different plants (and parts of plants) report to different tables at different aggregations, the generation and fuel data reported in one table does not always match the generation and fuel data reported in another. Although the GF table is more complete, the BF and G tables are more granular.

Thus, the data pipeline distributes the complete data from the GF table (at the prime mover level) to the generator level using data reported in the G and BF tables, as well as using data from other EIA tables. The allocation process assumes that:
* The GF table is the authoritative source of information about how much generation and fuel consumption is attributable to an entire plant.
* The G table is the authoritative source of information about how much generation is attributable to an individual generator, if it reports in that table.
* The BF table is the authoritative source of information about how much fuel consumption is attributable to an individual boiler, if it reports in that table.

Net generation data from the GF table is allocated to individual generators proportionally to the the reported generation data for each generator in the G table, if available, and otherwise is allocated proportionally to the nameplate capacity of each generator with the same fuel type and prime mover. Likewise, fuel consumption data from the GF table is allocated to individual generators proportionally to the the reported fuel consumption data for each boiler in the BF table, if available, and otherwise is allocated proportionally to the nameplate capacity of each generator with the same fuel type and prime mover.

Because fuel consumption in the boiler_fuel_eia923 table is reported per `boiler_id`, we must first map this data to generators using the Boiler-Generator Association table from EIA-860. For boilers that have a 1:m or m:m relationship with generators, we allocate the reported fuel to each associated generator based on the nameplate capacity of each generator. So if boiler "1" was associated with generator A (25 MW) and generator B (75 MW), 25% of the fuel consumption would be allocated to generator A and 75% would be allocated to generator B.

At the end of this process, each generator is associated with multiple records for each month, which represent the total amount of fuel consumption and net generation associated with each `energy_source_code` consumed by a generator. Before aggregating the data to the generator level, we perform several intermediate calculations that rely on information about the volume of each specific fuel consumed in a generator.

## Updating energy source codes

For each generator, the EIA-860 generator file lists up to 6 energy_source_code, which represents all of the energy sources consumed to produce power, as well as up to 6 startup energy source codes, which represent all of the energy sources consumed during plant startup. The unique set of `energy_source_code` and `startup_source_code` form part of the key on which generation and fuel data from the GF table is allocated to each generator.

Any "Other" (OTH) energy source codes are manually changed to specific fuel types using the methodology described [here](../Emissions%20Calculations/Assigning%20Energy%20Source%20Codes.md).

## Assigning a primary fuel value

While fuel-level data is available for each generator, we assign an annual primary fuel for each generator and plant using the methodology described [here](../Data%20Aggregation/Plant%20Primary%20Fuel.md).

## Calculating emissions
While fuel-level data is available for each generator, we calculate total emissions totals for each fuel consumed. We calculate:
 - GHG emissions ([methodology](../Emissions%20Calculations/GHG%20Emissions.md))
 - NOx emissions ([methodology](../Emissions%20Calculations/NOx%20Emissions.md))
 - SO2 emissions ([methodology](../Emissions%20Calculations/SO2%20Emissions.md))
 - Biomass-adjusted emissions ([methodology](../Emissions%20Calculations/Adjusting%20Emissions%20for%20Biomass.md))
 - CHP-adjusted emissions ([methodology](../Emissions%20Calculations/Adjusting%20Emissions%20for%20CHP.md))

## Validation Checks

While cleaning the EIA-923 data, we run several validation checks to ensure that the cleaned data does not contain any unexpected or anomalous values. These checks include:
 - Check that there are no missing energy source codes associated with records with non-zero fuel consumption.
 - Check that there are no negative fuel consumption or emissions values (net generation may be negative)

## Final data cleaning steps
Once all of the preceding steps are complete, we aggregate the EIA-923 data to the generator level. A `subplant_id` (see [methodology](../Data%20Aggregation/Subplant%20Aggregation.md)), `plant_primary_fuel`, and `prime_mover_code` are then assigned to each generator record.

