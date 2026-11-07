---
stoplight-id: hg_emissions
---

## Mercury (Hg) emissions
Mercury (Hg) is a toxic heavy metal that is released when fuels containing trace amounts of it, especially coal, are burned. Starting in 2015, coal- and oil-fired electric generating units larger than 25 MW were subject to the EPA's [Mercury and Air Toxics Standards (MATS)](https://www.epa.gov/stationary-sources-air-pollution/mercury-and-air-toxics-standards), and many of these units report hourly Hg emissions to the EPA.

OGE currently only includes the Hg emissions that are measured and reported under MATS. Unlike CO2, NOx, and SO2, Hg emissions are not estimated for units that don't report them, so Hg data is only available for some subplants and plants. This is consistent with eGRID, which reports Hg emissions for individual units but leaves its plant-level and regional Hg fields blank, because Hg data is not available for every unit at a plant.

### Which units report hourly Hg data
Hourly Hg data is reported by coal- and petroleum coke-fired units that monitor Hg with a continuous emissions monitor or sorbent traps. Several types of units do not report hourly Hg data:
- units that qualify as low emitting EGUs (LEE), which test their emissions periodically instead of monitoring them continuously
- oil-fired units, which demonstrate compliance through periodic stack tests
- natural gas units and combustion turbines, which are not subject to MATS
- units of 25 MW or smaller, and units that burn municipal solid waste or biomass

In addition, no Hg data is reported for hours when a monitor isn't working, and MATS does not require units to substitute data for these hours. In 2024, 344 units reported hourly Hg data, with total Hg emissions of 4,431 lb, but only 6 of these units reported Hg in every hour that they operated.

### Loading Hg data
We load hourly MATS data from PUDL's `core_epamats__hourly_emissions` table, which starts in 2015, and match it to the hourly CEMS data for each unit. Like other CEMS data, negative Hg values are treated as missing, and Hg data for plants that are removed from the CEMS data (such as plants that are not connected to the grid, or plants in Puerto Rico) is also removed. Hourly unit-level Hg data is included in the `cems_cleaned` output file.

### Identifying complete Hg data
Because Hg is only reported for some units, summing Hg data to subplants and plants could produce totals that silently leave out the emissions of units that don't report Hg. For example, at a plant with both a coal unit that reports Hg and a natural gas unit that doesn't, the sum of the plant's reported Hg would only represent the coal unit. So Hg data is only reported where it is complete.

Hg data is most complete for individual unit-hours, which are the least aggregated, and least complete for plant-years, which require complete data for every unit at the plant in every month of the year. Each level of aggregation adds its own requirement:

| Data | Hg data is reported when... |
|---|---|
| Unit-hour | the unit reported Hg in that hour |
| Subplant-month | the subplant reports Hg to MATS, every unit in the subplant that consumed fuel during the month reported Hg in at least one hour of that month, and all of the subplant's fuel consumption that month is reported in CEMS |
| Subplant-hour | the subplant-month is complete, and every unit in the subplant that consumed fuel in that hour reported Hg |
| Plant-hour | every subplant at the plant that consumed fuel in that hour has Hg data, and no subplant at the plant consumed fuel that month that is only reported in EIA-923 |
| Plant-month | every subplant at the plant that consumed fuel during the month has complete Hg data for that month |
| Subplant-year or plant-year | the Hg data is complete in every month of the year |

A few notes on these requirements:
- **Missing hours within a month**: Monthly and annual Hg values sum the hours in which Hg was reported. If a unit reports Hg in some, but not all, of the hours that it operated during a month, the month is still considered complete. This is consistent with how eGRID calculates annual Hg emissions for each unit. If every hour were required, almost no months would be complete.
- **Zero vs. blank values**: At subplants that report Hg, Hg is zero in hours when the subplant didn't consume any fuel. Hg is blank when a subplant or plant consumed fuel and doesn't have complete Hg data, and for subplants and plants that don't report Hg at all, such as natural gas, solar, or wind plants.
- **Non-CEMS data**: The hourly values of data that comes from EIA-923 are shaped using an hourly profile, so if a subplant consumed fuel that is only reported in EIA-923 in a month, the Hg data for its plant is blank in every hour of that month.
- **Plants with units that don't report Hg**: Because the requirements apply to each hour and month, a plant with units that don't report Hg can still have Hg data in the hours and months when those units weren't operating. For example, a coal plant with a natural gas peaking unit that doesn't report Hg has hourly Hg data in hours when the peaking unit is off, and monthly Hg data in months when it didn't run at all. However, it won't have annual Hg data if the peaking unit operated at any point during the year.

Because a plant-hour requires complete data from every subplant at the same time, but a plant-month only requires each unit to report Hg at some point during the month, the hourly plant data for a year may add up to less than the annual plant data. In 2024, the hourly data of the plants with annual Hg data added up to 3,704 lb, 2% less than their annual data (3,782 lb).

In 2024, Hg data was complete for 311 subplants (4,332 lb) and 103 plants (3,782 lb) for the entire year, for 1,358 plant-months at 130 plants, and for at least some hours at 147 plants.

### Hg in the published data
- Hg is reported as `hg_mass_lb` (total Hg) and `hg_mass_lb_for_electricity` (Hg adjusted for CHP using the electric allocation factor, see [Adjusting Emissions for CHP](Adjusting%20Emissions%20for%20CHP.md)). There is no biomass adjustment for Hg, so there are no `_adjusted` Hg columns. eGRID does not adjust Hg for CHP.
- Hg data is included in the hourly `cems_subplant` output file, and in the hourly, monthly, and annual plant and subplant data.
- Because Hg data is not complete for every plant, Hg is not included in the power sector data or the carbon accounting data, and no Hg emission rates are calculated.
- The share of Hg mass that was measured or imputed under MATS is summarized in the CEMS pollutant measurement quality metrics (see [Data Quality Metrics](../../Data%20Validation/Data%20Quality%20Metrics.md)).

## Future Work, Known Issues, and Open Questions
- Hg is not yet estimated for units that don't report hourly Hg data. In the future, these values could be estimated from each unit's monthly fuel consumption using Hg emission factors from the EPA's AP-42 or other sources and the Hg control data reported in EIA-860 and EIA-923, similar to how NOx and SO2 are estimated for units that don't report to CEMS. This would allow Hg to be reported for every plant and included in the power sector data.
- Missing Hg values in hours when a unit operated are not filled.
- Several large coal units, such as those at Intermountain (6481), Hunter (6165), and Bonanza (7790), report hourly CEMS data but no MATS Hg data, likely because they qualify as low emitting EGUs. eGRID's unit file has the same gap.
