---
stoplight-id: gross_to_net
---

## Background of Gross Generation and Net Generation
Gross generation describes the total amount of electricity, in MWh, that is generated as measured at the terminal of a generator. However, the amount of generation that a power plant injects into the grid is referred to as “net generation,” which is gross generation net any losses between the generator and the point at which the power plant is interconnected to the grid. These losses result from electricity consumed by auxiliary equipment used to generate electricity (pumps, compressors, feed systems), control emissions, and operate the building (lighting, air conditioning, SCADA systems), as well as electrical losses in the high voltage power system ([EPRI 2011](https://www.epri.com/research/products/1024651)). These losses are known by various names including “house loads,” “station use,” “parasitic loads,” and “auxiliary loads.” Certain losses, such as those related to drivepower systems, scale with gross generation, and can be thought of as an efficiency penalty for gross generation. Other losses, such as those related to station use for controls and building operation, can be thought of as “always on” loads that are independent of the amount of electricity being generated at any moment.

Generally in power systems emissions analysis, we primarily care about net generation. The monthly generation reported in EIA-923 represents net generation, but the hourly generation reported in CEMS represents gross generation. In order to calculate hourly production-based emission factors, hourly gross generation must be converted to hourly net generation.

## General approaches to converting gross to net generation

The general approach for this gross-to-net generation conversion requires matching gross generation reported in CEMS with net generation reported in EIA-923 for each subplant, and calculating a conversion factor between the values. The goals of the gross-to-net generation conversion are twofold:
1. Ensure that the calculated hourly net generation values sum at the annual level to total net generation reported in EIA-923 for each subplant or plant
2. Ensure that the hourly shape of the calculated net generation is realistic

There are three primary methods that we have identified for calculating this gross-to-net conversion factor, each of which have different advantages/disadvantages in achieving the above two goals:

**Method A: Gross-to-net ratio**.
* **How the factor is calculated**: This approach involves dividing net generation by gross generation.
* **How it is used to calculate net generation**: This ratio would be multiplied by hourly gross generation to estimate hourly net generation.
* **Implicit assumptions**: This approach assumes that all losses scale with generation, such that if gross generation in an hour were zero, net generation would also be zero.
* **Advantages**: Calculated net generation values will always match reported net generation values in total. While this method does not reflect any static losses, it is likely that most power plant losses vary with generation.
* **Disadvantages**: This method cannot be used for converting positive gross generation to negative net generation, since it will result in a negative ratio that will invert the hourly shape of the gross generation profile. This method also has the disadvantage of potentially distorting the generation shape in some cases if the ratio is large (>>1, which leads to peaky profiles) or small (close to zero, which leads to flattened profiles).
* **Previous uses in the academic literature**: [Graff-Zivin et al. (2014)](http://www.sciencedirect.com/science/article/pii/S0167268114000808), [Bushnell and Novan (2021)](https://doi.org/10.1086/713249), [Cicala (2022)](https://pubs.aeaweb.org/doi/10.1257/aer.20172034)

**Method B: Gross-to-net shift factor**.
* **How the factor is calculated**:This approach involves subtracting gross generation from net generation, and then dividing by the number of hours in the period to which it will be applied (i.e. if annual gross generation were subtracted from annual net generation, the resulting value would be divided by 8,760 since that is the number of hours in a non-leap year).
* **How it is used to calculate net generation**: This shift factor would then be added to hourly gross generation to estimate hourly net generation.
* **Implicit assumptions**: This approach assumes that all losses are independent of generation, and are static throughout the period. Thus, if gross generation in an hour were zero, net generation would be negative.
* **Advantages**: Calculated net generation values will always match reported net generation values in total. In addition, this method always preserves the shape of the gross generation profile, since it only shifts the profile up or down along the y axis. This makes it useful for preserving the shape of gross generation even when total net generation is negative.
* **Disadvantages**: It is unlikely that this method’s implicit assumption is realistic in most situations. For plants that tend to cycle on and off, it also performs poorly because during startup periods of low gross generation, it tends to calculate negative net generation, which can lead to negative emission factors in those hours.
* **Previous uses in the academic literature**: None yet identified

**Method C: Gross-to-net regression**.
* **How the factor is calculated**: This approach involves calculating a linear regression of the form _N<sub>m</sub>_ = _G<sub>m</sub>_ _r_ + _s_, where _N<sub>m</sub>_ represents the average hourly net generation in month _m_ (in units of MW), _G<sub>m</sub>_ represents the average hourly gross generation in month _m_ (in MW), _r_ represents the gross-to-net ratio and _s_ represents the gross-to-net shift factor. Monthly gross and net generation totals (in MWh) are divided by the number of hours in each month before regressing so that the resulting shift factor _s_ represents the hourly shift value rather than the monthly shift value.
* **How it is used to calculate net generation**: Net generation would be calculated by multiplying hourly gross generation by _r_ and then adding _s_.
* **Implicit assumptions**: This approach assumes that there are both losses that scale linearly with generation and static losses that are independent of generation. However, there may be multiple other predictors of this relationship which may or may not be related in a linear fashion.
* **Advantages**: Of these three methods, the regression approach probably best represents the physical reality of plant operations (that losses are some combination of variable and static losses). We find for many subplants that the adjusted r2 values are quite high (> 0.9). Because this conversion includes a static loss term, it means that net generation can be negative in certain hours.
* **Disadvantages**: Although the linear model seems to fit the data well, the calculated net generation, when summed to the annual level, will not match the total volume of net generation reported in EIA-923, thus failing to achieve the first goal of this conversion. In addition, for some (sub)plants the static loss term s is positive, which does not have an intuitive physical explanation (it would mean that a plant could have positive net generation even when gross generation is zero). Due to these disadvantages, these regression values are not currently used in the data pipeline (although the values are calculated and exported as intermediate outputs).
* **Previous uses in the academic literature**: None yet identified

## Specific approach used to convert gross to net generation
Following the approach primarily used in previous academic literature, the OGE data pipeline primarily relies on using a gross-to-net ratio to convert gross to net generation.

**Resolution for conversion factors**

Ratios are calculated both at the subplant and plant level and both at the monthly and annual resolution. Monthly-resolution conversion factors are not currently used in the data pipeline because we are not always confident that the monthly values reported in EIA-923 are accurate. For example, in some cases, generators report annual totals and EIA distributes the data to each month. The use of annual-resolution conversion factors means that the annual total calculated net generation will match the annual total reported net generation, but monthly total calculated net generation will not always match monthly total reported net generation.

**Hierarchy for applying gross to net conversion factors**

When calculating net generation, the pipeline uses the best available conversion factor for each subplant. Certain conversion factors may not be available for a subplant for multiple reasons, including if matching net generation data in EIA-923 cannot be identified, or if the factor is anomalous and filtered out (see next section). In order to ensure that the total calculated net generation matches the total reported net generation, the same conversion method must be used for all subplants and months within a plant. Thus, for each plant, the pipeline uses the best method that is available for all subplants and all months at that plant.

The conversion factors are applied in the following hierarchical order:

1. Subplant-specific, annual gross-to-net ratio
2. Plant-specific, annual gross-to-net ratio
3. Subplant-specific, annual gross-to-net shift factor
4. Plant-specific, annual gross-to-net shift factor
5. Fleet-specific, annual gross-to-net ratio
6. Prime mover-specific, default gross-to-net ratios from EIA 
7. Apply an assumed gross to net ratio of 0.97

If no plant-specific factors are available (methods 1-4), the pipeline uses a fleet-specific ratio that represents the average gross-to-net ratio for all subplant (nationally) that consume the same fuel category (natural gas, coal, etc). Otherwise, we use default gross to net ratios published in the EIA Electric Power Monthly Technical Notes (Appendix C). As a final backstop, we use an assumed gross to net ratio of 0.97

Shift factors are used explicitly instead of ratios in cases where the reported gross generation data in CEMS is zero, but EIA reports non-zero net generation data. In these cases, multiplying by a ratio would result in 0 calculated net generation. Applying a shift factor (which is added to the gross generation) effectively applies a flat generation profile to the data, but ensures that the annual total net mwh will match the data in EIA. 

The following table shows what percent of gross generation reported in CEMS was converted to net generation using each method.


<table>
  <tr>
   <td>Method
   </td>
   <td>2021
   </td>
   <td>2022
   </td>
  </tr>
  <tr>
   <td>1. Subplant ratio
   </td>
   <td>87.2%
   </td>
   <td>87.4%
   </td>
  </tr>
  <tr>
   <td>2. Plant ratio
   </td>
   <td>9.6%
   </td>
   <td>9.0%
   </td>
  </tr>
  <tr>
   <td>3. Fleet-specific ratio
   </td>
   <td>2.6%
   </td>
   <td>2.9%
   </td>
  </tr>
  <tr>
   <td>4. Default EIA ratio
   </td>
   <td>0.6%
   </td>
   <td>0.6%
   </td>
  </tr>
  <tr>
   <td>5. Assumed 97% ratio
   </td>
   <td>0.1%
   </td>
   <td>0.1%
   </td>
  </tr>
</table>


**Filtering of GTN Ratios**

Before calculating net generation, the conversion factors are filtered to remove any factors that would lead to anomalous net generation values being calculated.

Based on analysis of historical data, about 90% of the gross to net generation values fall between approximately 0.5 and 1.25.

The first filter removes factors that would lead to net generation values that significantly exceed the subplant’s nameplate capacity. In this case, the factor is filtered out if the 98th percentile of calculated hourly net generation values in a month exceed 125% of the subplant’s nameplate capacity. We use the 98th percentile instead of the maximum value to allow for a small handful of hours (approximately 14 hours in a 30-day month) to exceed the value. The use of the 125% exceedance threshold allows for the fact that a plant’s nameplate capacity can vary throughout the year and sometimes be exceeded, and the use of the 98th percentile prevents a small number of anomalous hours from causing the factor to be filtered out.

The second filter removes any ratios that are less than 0.5 (including negative ratios), since multiplying a negative ratio by gross generation would invert the shape of the hourly gross generation data.

After these filters have been applied for individual factors for each subplant-month, the pipeline removes all factors of a certain type for an entire plant if there are any missing factors for any subplant-month. For example, if a subplant ratio for a single month at a single subplant is filtered out, all subplant ratios for all other months and subplants at the same plant will also be removed.