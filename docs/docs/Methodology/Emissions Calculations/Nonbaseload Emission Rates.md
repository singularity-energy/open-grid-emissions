---
stoplight-id: nonbaseload_emission_rates
---

## Nonbaseload emission rates

Nonbaseload emission rates are a type of marginal emission rate which estimate the emission rate of the generation that is most likely to change when demand changes, such as when energy efficiency or renewable energy displaces other generation. They are sometimes used to estimate the emissions avoided by such projects. They should not be used to assign emissions to electricity use in carbon footprints or GHG inventories (scope 2 emissions); use OGE's [consumed emission rates](../Emissions%20Calculations/Consumption-based%20Emissions.md) or total generated emission rates for that.

OGE calculates nonbaseload rates using the method used in EPA's eGRID. OGE's results include each plant's nonbaseload factor, and the [eGRID-format workbook](../../Overview/Using%20the%20Data.md#files) includes nonbaseload generation, emission rates, generation by fuel type, and resource mix for each state, balancing authority, NERC region, and the U.S.

### Capacity factor
Each plant's annual capacity factor (`capacity_factor` in the annual plant data) is its annual net generation divided by what it could have generated at nameplate capacity while its generators were operating, leaving out energy storage:

capacity factor = net generation of the plant's subplants other than storage (MWh) &frasl; &sum; (nameplate capacity (MW) of each generator other than storage &times; hours in the months the generator operated)

A generator is counted as operating in a month if either:
- its EIA-860 operating and retirement dates cover the month, i.e. it is existing or retired, it began operating by the end of the month, and it had not retired before the month began. A generator that starts or retires during a month counts for the whole month, since EIA-860 reports these dates by month.
- it produced electricity in the month (positive net generation in the EIA-923 data allocated to generators), even if it has not reached its commercial operation date. Generators that are still being tested before commercial operation (EIA-860 status `TS`) are reported as proposed, but can generate for months before they are considered operating.

Energy storage, including pumped storage, is left out of both the net generation and the capacity: storage subplants are those whose primary prime mover is a storage technology, and storage generators are those with a storage prime mover. Storage mostly shifts electricity from one time to another, and its net generation is mostly the energy lost in charging and discharging, so counting its capacity would lower the capacity factor of the plant's other generators and make them look more nonbaseload than they are. Plants that only have storage have a blank capacity factor. eGRID counts storage in both. In 2024, leaving out storage gives 366 storage-only plants a blank capacity factor instead of zero, and changes the nonbaseload factor of 14 plants that also have batteries, reducing U.S. nonbaseload generation by another 1.0 TWh on top of the change from counting only operating months (below). For example, Moss Landing (260) has 750 MW of batteries next to 1,398 MW of gas generators: its capacity factor is 0.21 with the batteries and 0.33 without them, and its nonbaseload factor is 0.99 with them and 0.79 without them. A battery that is grouped into the same subplant as other generators can't be separated from them, which was the case for one subplant in 2024.

The number of hours in a month accounts for leap years, so a full year has 8,784 hours in a leap year and 8,760 otherwise. As in eGRID, the capacity factor is zero for plants with negative net generation. Plants that have no generator data with a capacity are assumed to have operated their whole nameplate capacity (`capacity_mw`, which includes any storage) for the whole year. The capacity factor is blank for plants without a nameplate capacity.

eGRID divides by the plant's nameplate capacity times 8,760 hours, so plants that started operating or retired during the year look mostly idle, and plants with generators that are still testing can have capacity factors far above 1. For example, Roseland Solar (65028) generated 950 GWh in 2024 from a 500 MW solar generator that was still testing all year, while its only operating generator was a 51.8 MW battery. Its capacity factor is 2.09 when only the battery counts, and 0.22 when only the testing solar generator counts. In 2024, counting only operating months (with storage still included) changes the capacity factor of 591 plants and the nonbaseload factor of 25, reducing U.S. nonbaseload generation by 4.4 TWh (0.4%) compared with eGRID's full-year capacity factor.

Capacity factors above 1 point to problems in the input data, such as missing capacity or net generation that is too high. They are not corrected, but the pipeline lists them in a warning.

### Nonbaseload factor
The capacity factor is used as a surrogate for where a plant falls in the dispatch order: plants that run most of the time are considered baseload, and plants that run rarely are considered to follow load. Each plant's nonbaseload factor (`nonbaseload_factor` in the annual plant data) is the share of its generation and emissions that is considered nonbaseload:

| Capacity factor | Nonbaseload factor |
|---|---|
| Less than 0.2 | 1 (all generation is nonbaseload) |
| 0.2 to 0.8 | -5/3 &times; capacity factor + 4/3, which falls linearly from 1 to 0 |
| Greater than 0.8 | 0 (all generation is baseload) |

The nonbaseload factor is blank for:
- plants whose primary fuel is geothermal (`GEO`), energy storage (`MWH`), nuclear (`NUC`), purchased steam (`PUR`), solar (`SUN`), hydro including pumped storage (`WAT`), or wind (`WND`). These resources are not dispatched to follow load, so their capacity factor does not indicate whether they are on the margin.
- plants whose capacity factor is zero or blank, including plants with negative net generation, which have no nonbaseload generation.

The factor uses each plant's [primary fuel](../Data%20Aggregation/Plant%20Primary%20Fuel.md), which OGE assigns differently from eGRID, so a few plants are excluded in one dataset and not the other.

The nonbaseload factor applies to the whole plant. Subplants and months do not have their own factors.

### Regional nonbaseload values
Each plant's nonbaseload generation and emissions are its annual values times its nonbaseload factor. Plants without a factor have no nonbaseload generation or emissions. In the eGRID-format workbook, the regional values are summed from the plant sheet:

- **Nonbaseload generation** is the sum of the plants' nonbaseload generation (`PLNGENNB`).
- **Nonbaseload output emission rates** divide the sum of the plants' nonbaseload emissions by the region's nonbaseload generation, using the same emissions as OGE's other eGRID rates (adjusted for biomass and CHP). The ozone season NOx rate applies each plant's annual factor to its ozone season (May to September) NOx and generation. The Hg rate is blank, as in eGRID.
- **Nonbaseload generation by fuel type** assigns each plant's whole nonbaseload generation to the fuel type of its plant fuel category, as eGRID does. Since nuclear, hydro, wind, solar, and geothermal plants have no nonbaseload factor, those fuel types always have zero nonbaseload generation.
- **Nonbaseload resource mix** divides each fuel type's nonbaseload generation by the region's nonbaseload generation.

Regions without any plant with a nonbaseload factor have blank nonbaseload values.

Because the rates and the generation come from the same plant values, a region's nonbaseload generation always equals the sum of its plants' nonbaseload generation, and its nonbaseload emissions divided by its rate. eGRID's 2024 and 2025 data files are not consistent in this way, so their published regional nonbaseload rates are higher than the rates calculated from their own plant data (see [Comparison to eGRID](../../Overview/Comparison%20to%20eGRID.md#2b-egrid2023-bugs)).

### Future work
- A single annual factor treats a plant that runs all summer and is offline the rest of the year the same as a plant that runs at part load all year. Monthly or subplant-level factors could better separate baseload from nonbaseload generation.
