---
stoplight-id: primary_fuel
---

## Assigning Plant Primary Fuel

The primary fuel of a plant is assigned based on the `energy_source_code` that has the highest annual volume of `fuel_consumed_mmbtu`, as reported in the EIA-923 generation and fuel table.

If no fuel consumption data is reported or if there are multiple fuels with the same volume of annual fuel consumption, the primary fuel is based on the `energy_source_code_1` associated with the highest combined generator nameplate capacity, as reported by the EIA-860 Generators file.

If multiple fuels are associated with the same amount of nameplate capacity at a plant, then the primary fuel is determined based on the energy source code that is associated with the highest annual volume of net generation.

In rare cases where all of the preceding methods of determining plant primary fuel fail, the plant primary fuel is determined based on the primary fuel associated with the most number of generators at a plant.

The same methodology is used to assign a primary fuel to each subplant (see [subplant aggregation](Subplant%20Aggregation.md)), using only the generators that belong to that subplant.

## Assigning Primary Prime Mover

While most energy storage resources report an energy source code of `MWH`, not all do. For example, pumped storage hydroelectric plants typically report an energy source code of `WAT` (water). The primary fuel alone therefore cannot reliably identify energy storage, so we also assign a primary prime mover to each subplant and plant to assist with this identification.

The primary prime mover is the `prime_mover_code` associated with the highest combined generator nameplate capacity at each subplant or plant, based on data from the EIA-860 Generators file. If multiple prime movers are associated with the same amount of nameplate capacity, the prime mover associated with the greatest number of generators is used. If a tie still remains, the prime mover that comes first alphabetically is used so that the result is consistent every time the data is generated.

The primary prime mover is assigned independently of the primary fuel, so the two may describe different generators at the same plant. For example, a natural gas plant whose battery storage capacity is greater than its combustion turbine capacity will have a primary fuel of `NG` but a primary prime mover of `BA` (battery), and will be categorized as a storage plant (see below).

## Assigning Fuel Categories

Each plant and subplant is assigned a `fuel_category` based on its primary fuel, using [this mapping table](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/oge/reference_tables/energy_source_groups.csv). Any plant or subplant whose primary prime mover is one of the following energy storage prime movers is then assigned a `fuel_category` of `storage`, regardless of its primary fuel:

Prime mover code | Description
---------|---------
BA | Battery energy storage
CE | Compressed air energy storage
ES | Energy storage, other
FW | Flywheel energy storage
PS | Pumped storage hydroelectric

Plants and subplants with a primary fuel of `MWH` are assigned a `fuel_category` of `storage`, regardless of prime mover. 
