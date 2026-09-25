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

The primary prime mover is assigned after the primary fuel, so that both describe the same generators. The primary prime mover is the `prime_mover_code` associated with the highest combined generator nameplate capacity at each subplant or plant, based on data from the EIA-860 Generators file, considering only the generators whose primary energy source code matches the primary fuel. For example, a hybrid solar plant whose battery storage capacity is greater than its PV capacity will have a primary fuel of `SUN` and a primary prime mover of `PV`.

If multiple prime movers are associated with the same amount of nameplate capacity, the prime mover associated with the greatest number of generators is used. If a tie still remains, the prime mover that comes first alphabetically is used so that the result is consistent every time the data is generated.

In some cases, no generator's energy source code matches the primary fuel. This can happen when the energy source code reported in EIA-923 differs from the one reported in EIA-860 (for example, a coal unit that reports burning subbituminous coal in EIA-923 but lists bituminous coal as its primary energy source in EIA-860). In these cases, the primary prime mover is the prime mover associated with the highest combined nameplate capacity across all generators at the subplant or plant.

## Assigning Fuel Categories

Each plant and subplant is assigned a `fuel_category` based on its primary fuel, using [this mapping table](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/oge/reference_tables/energy_source_groups.csv). Any plant or subplant whose primary prime mover is one of the following energy storage prime movers is then assigned a `fuel_category` of `storage`, regardless of its primary fuel:

Prime mover code | Description
---------|---------
BA | Battery energy storage
CE | Compressed air energy storage
ES | Energy storage, other
FW | Flywheel energy storage
PS | Pumped storage hydroelectric

Compressed air energy storage (`CE`) is the exception. It is only assigned a `fuel_category` of `storage` if its primary fuel is `MWH`. Older compressed air storage plants burn natural gas to supplement the stored energy when discharging, so these keep the fuel category of their primary fuel (for example, `natural_gas`), and their emissions are reported with that fuel category. They are identified as `hybrid` storage instead (see [Energy Storage](Energy%20Storage.md)).

Plants and subplants with a primary fuel of `MWH` are assigned a `fuel_category` of `storage`, regardless of prime mover.

Plants and subplants that contain energy storage are also assigned a storage type (standalone, co-located, or hybrid), which is described in [Energy Storage](Energy%20Storage.md).
