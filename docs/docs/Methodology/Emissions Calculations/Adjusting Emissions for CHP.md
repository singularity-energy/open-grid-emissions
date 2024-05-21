---
stoplight-id: chp_adjustments
---

## Background on Combined Heat and Power generation
Combined Heat and Power (CHP) plants produce both electricity and useful thermal energy in the form of heat or steam that may be used for industrial, commercial, or space conditioning purposes. Thus, when considering the emission intensity of generated electricity, the portion of the fuel that was consumed to produce useful thermal output should be excluded from this calculation. At its most basic, this adjustment involves calculating an electric allocation factor (which represents the percentage of fuel consumption was used to generate electricity) and then multiplying it by fuel and emissions totals for each plant. All totals that have been adjusted in this way use the `_for_electricity` string in the column name.

The formula and assumptions used for calculating the electric allocation factor are copied directly from the [eGRID2020 methodology](https://www.epa.gov/system/files/documents/2022-01/egrid2020_technical_guide.pdf) and is summarized below (track progress on updates to this methodology [here](https://github.com/singularity-energy/open-grid-emissions/issues/23)).

Although both total `fuel_consumed_mmbtu` and `fuel_consumed_mmbtu_for_electricity` are reported in EIA-923, the eGRID methodology does not directly use these values to calculate an electric allocation factor. Although this is not explicitly explained in the eGRID methodology, this is ostensibly because they are trying to allocate fuel based on the proportion of energy outputs, not energy inputs. This means that the `fuel_consumed_mmbtu_for_electricity` values calculated using this method will not match the reported `fuel_consumed_mmbtu_for_electricity` values from EIA-923.

When performing the below calculations, all calculations are performed at the subplant level, rather than at the generator or unit level. This is due to the fact that in some CHP systems, all of the reported fuel consumption at the subplant may be associated with one of the subplant's generators, while all of the reported generation data may be associated with another of the subplant's generators. 

### Calculating useful thermal output
In form EIA-923, generators report both the total quantity of each fuel consumed, and the quantity consumed for electric generation. By subtracting these two values, we can calculate the total fuel consumed for heating purposes:

`fuel_consumed_for_heating_mmbtu = fuel_consumed_mmbtu - fuel_consumed_for_electricity_mmbtu`

To calculate the gross thermal output, we use an assumed efficiency factor of 80%:

`gross_thermal_output_for_heating_mmbtu = 0.8 * fuel_consumed_for_heating_mmbtu`

However, not all gross thermal output can be captured for useful purposes such as space heating or industrial processes. Thus, the gross thermal output is further adjusted using an assumed efficiency factor of 75%:

`useful_thermal_output_mmbtu = 0.75 *gross_thermal_output_for_heating_mmbtu`

This assumed efficiency factor assumes that the CHP units generate electricity first and use the waste heat for other purposes, also known as “topping.” While there are some units that generate and use heat first and then use the waste heat to generate electricity, also known as “bottoming, data from the EIA shows that over 90% of CHP facilities are topping facilities (track progress on updating these assumptions [here](https://github.com/singularity-energy/open-grid-emissions/issues/23)).

### Calculating the electric allocation factor
The electric allocation factor is the ratio of the electricity output to the sum of the combined electricity and thermal outputs. Because electricity net generation is in units of MWh, it must be converted to mmbtu using a conversion factor of 3.412142 mmbtu per mwh.

`Electric_Allocation_Factor =  (3.412142 * net_generation_mwh) / (useful_thermal_output_mmbtu + (3.412142 * net_generation_mwh))`

An electric allocation factor is calculated and applied to all generators, regardless of reported CHP status. In the case that total fuel consumed equals total fuel consumed for electricity, the useful thermal output will be zero, and thus the electric allocation factor will be 100%.

## Adjusting CEMS data for CHP
Because the calculation of the electric allocation factor requires `net_generation_mwh` and `fuel_consumed_for_electricity_mmbtu`, neither of which are reported in CEMS, several additional steps are required to adjust CEMS data for CHP.

First, we estimate `fuel_consumed_for_electricity_mmbtu` using data reported in EIA-923. For each subplant-month, we calculate the ratio between `fuel_consumed_for_electricity_mmbtu` and `fuel_consumed_mmbtu`, as reported in EIA-923. This subplant-month specific ratio is then multiplied by the reported hourly `fuel_consumed_mmbtu` in CEMS to estimate hourly `fuel_consumed for_electricity` ([Deetjen and Azevedo 2019](https://doi.org/10.1021/acs.est.9b02500)). If a subplant-specific ratio cannot be calculated, we use a plant-month specific ratio. Once this value is calculated, we can calculate `useful_thermal_output_for_heating_mmbtu` as described above.

Since electric allocation factor also relies on `net_generation_mwh`, we must first convert `gross_generation_mwh` to `net_generation_mwh`. That methdology is described in detail [here](../Converting%20Gross%20to%20Net%20Generation.md).

For each CEMS subplant, a unique electric allocation factor is calculated for for each hour.

Once the electric allocation factor is available, a new value for `fuel_consumed_for_electricity_mmbtu` is calculated by multiplying the factor by the original `fuel_consumed_mmbtu` value.

## Future Work, Known Issues, and Open Questions
- Update assumed efficiency factors for calculating the electric allocation factor ([details](https://github.com/singularity-energy/open-grid-emissions/issues/23))
- Update electric allocation factor calculation to consider topping vs bottoming cycles ([details](https://github.com/singularity-energy/open-grid-emissions/issues/23))
- Confirm how CEMS data reports fuel consumed for electricity ([details](https://github.com/singularity-energy/open-grid-emissions/issues/23)) and how to interpret reported steam load ([details](https://github.com/singularity-energy/open-grid-emissions/issues/103))