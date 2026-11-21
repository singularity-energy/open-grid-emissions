---
stoplight-id: nonbaseload_emission_rates
---

## Nonbaseload emission rates

Nonbaseload emission rates are a type of marginal emission rate which estimate the emission rate of the generation that is most likely to change when demand changes, such as when energy efficiency or renewable energy displaces other generation. They are sometimes used to estimate the emissions avoided by such projects. They should not be used to assign emissions to electricity use in carbon footprints or GHG inventories (scope 2 emissions); use OGE's [consumed emission rates](../Emissions%20Calculations/Consumption-based%20Emissions.md) or total generated emission rates for that.

OGE calculates nonbaseload rates using the method used in EPA's eGRID. OGE's results include each plant's nonbaseload factor, and the [eGRID-format workbook](../../Overview/Using%20the%20Data.md#files) includes nonbaseload generation, emission rates, generation by fuel type, and resource mix for each state, balancing authority, NERC region, and the U.S.

### Capacity factor
Each plant's annual capacity factor (`capacity_factor` in the annual plant data) is its annual net generation divided by its nameplate capacity times the number of hours in the year:

capacity factor = net generation (MWh) &frasl; (nameplate capacity (MW) &times; hours in the year)

The number of hours in the year is 8,784 in a leap year and 8,760 otherwise. As in eGRID, the capacity factor is zero for plants with negative net generation, and it is blank for plants without a nameplate capacity.

The nameplate capacity is the plant's total capacity in the year (`capacity_mw`), which includes generators that started operating or retired during the year. Plants that operated for only part of the year therefore have a lower capacity factor than their operation while online would suggest.

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
- Plants that started operating or retired during the year are mostly counted as nonbaseload because of their low annual capacity factor. Counting only the hours in which each generator was operating would give these plants a more accurate capacity factor.
- A single annual factor treats a plant that runs all summer and is offline the rest of the year the same as a plant that runs at part load all year. Monthly or subplant-level factors could better separate baseload from nonbaseload generation.
