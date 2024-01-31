---
stoplight-id: egrid_comparison
---
Although the OGE methodology is based on the EPA's eGRID methodology, there are some important differences. Thus, if comparing OGEI data to eGRID data, it is important to keep the following differences in mind:


<table>
  <tr>
   <td>
Method/approach
   </td>
   <td>eGRID2020
   </td>
   <td>Open Grid Emissions
   </td>
   <td>Impact
   </td>
  </tr>
  <tr>
   <td>CAMD data used
   </td>
   <td>Uses annual total CEMS data
   </td>
   <td>Uses hourly CEMS data
   </td>
   <td>Allows OGE to identify missing and anomalous values in CEMS
   </td>
  </tr>
  <tr>
   <td>Allocation of EIA-923 fuel and generation data to generators
   </td>
   <td>Based on generator nameplate capacity
   </td>
   <td>Based on reported fuel and generation data if available, then based on nameplate capacity
   </td>
   <td>May result in different calculated values
   </td>
  </tr>
  <tr>
   <td>Assigning geothermal geotype (flash, steam, or binary) to geothermal plants
   </td>
   <td>Assigns a single geotype to each plant
   </td>
   <td>Assigns a generator-specific geotype (there are sometimes multiple types at a single plant), and uses updated geotype classifications based on most recent EIA data
   </td>
   <td>Geothermal CO2 emissions about 3% lower than values calculated in eGRID2020
   </td>
  </tr>
  <tr>
   <td>Fuel cell emissions
   </td>
   <td>Assumes zero co2 emissions for generators with fuel cell prime movers
   </td>
   <td>Applies fuel-specific emission factor for fuel cells that consume fossil fuels
   </td>
   <td>CO2 emissions for plants with a fuel cell prime mover are more than 800% higher than eGRID values.
   </td>
  </tr>
  <tr>
   <td>NOx emission factor for flared landfill gas
   </td>
   <td>0.02 lb NOx per MMbtu
   </td>
   <td>0.078 lb NOx per MMbtu
   </td>
   <td>Adjusted emissions from LFG will be lower than the eGRID values.
   </td>
  </tr>
  <tr>
   <td>CHP electric allocation factor
   </td>
   <td>A single annual average electric allocation factor is calculated
   </td>
   <td>Month-specific electric allocation factors are calculated
   </td>
   <td>May result in different calculated values
   </td>
  </tr>
  <tr>
   <td>Plant primary fuel determination
   </td>
   <td>Based on fuel consumption reported in CAMD data, otherwise on nameplate capacity.
<p>
Uses MSW rather than MSN or MSB.
   </td>
   <td>Same method, but based on fuel consumption reported in EIA-923. Also uses net generation as a tiebreaker when nameplate capacity for two fuels is equal.
   </td>
   <td>Certain plants are assigned a different primary fuel code
   </td>
  </tr>
  <tr>
   <td>OTH fuel types
   </td>
   <td>Updates OG fuel types, but not OTH fuel types
   </td>
   <td>Manually assigns an energy source code to generation with OTH fuel type
   </td>
   <td>Improves coverage of emissions data for these plants (there is no emission factor for OTH fuel, so these emissions would otherwise be zero)
   </td>
  </tr>
  <tr>
   <td>Global Warming Potential (GWP)
   </td>
   <td>Has used AR4 GWPs since eGRID2018 (still using AR4 as of eGRID2022)
   </td>
   <td>AR5 GWPs are used starting in 2019 (the earliest year of OGE data available, although they apply as far back as 2014). AR6 GWPs have been used since data year 2021. 
   </td>
   <td>CO2-eq factors from eGRID will underestimate the GWP of CH4, and overestimate the GWP of N2O relative to the currently-recognized GWPs of these gases.
   </td>
  </tr>
</table>