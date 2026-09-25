---
stoplight-id: energy_storage
---

## Identifying Energy Storage Resources

Energy storage resources do not generate electricity from a fuel. Instead, they charge from the grid or from a co-located generator and later discharge that energy. How a storage resource relates to other generators also affects how its data is reported. For example, a battery that is metered together with a solar array may have its discharge reported as part of the solar generation. To help users interpret and aggregate data for these resources, OGE identifies energy storage generators and categorizes each one by how it is connected to other generators.

Storage generators are identified based on the `prime_mover_code` reported for each generator in the EIA-860 Generators file for the data year:

Prime mover code | Description
---------|---------
BA | Battery energy storage
CE | Compressed air energy storage
ES | Energy storage, other
FW | Flywheel energy storage
PS | Pumped storage hydroelectric

Storage resources are defined specifically as resources that use electrical energy as an input (either from the grid or a specific generator). This means that (concentrated) solar thermal generators with integrated storage (prime mover `CP`) are not classified as energy storage since they directly store solar thermal energy prior to any electricity generation occuring. 

## Storage Categories

Each storage resource is assigned one of three storage categories:

Storage category | Description | Example
---------|---------|---------
`standalone` | The storage is an independent power plant, with no non-storage generators at the same plant or location. | A battery that is its own plant and charges from the grid
`co_located` | The storage is at the same plant or physical location as a non-storage generator, but is metered separately from it. | A battery at a solar or wind plant, or at a natural gas plant
`hybrid` | The storage is an integral part of a generator and is metered together with it. | A compressed air energy storage plant that burns natural gas when discharging, or a battery that is tightly DC-coupled to a solar array

## Assigning a Storage Category to Each Storage Generator

Each storage generator is assigned the storage category of the first of the following rules that applies. The rule that was used for each subplant is reported in the `subplant_storage_category_method` column of the subplant attributes table.

\# | `storage_category_method` | Storage category | Rule
---------|---------|---------|---------
1 | `hybrid_prime_mover` | `hybrid` | The generator's prime mover is `CE` (compressed air), and it reports an energy source other than `MWH`. This indicates that the storage is integrated into a generator that uses another energy source, such as compressed air storage that burns natural gas when discharging.
2 | `dc_coupled_tightly` | `hybrid` | The storage is reported as tightly DC-coupled in the EIA-860 Energy Storage file.
3 | `same_plant` | `co_located` | The same plant has at least one non-storage generator that is operating in the data year.
4 | `direct_support_other_plant` | `co_located` | The storage is reported as directly supporting a generator, and the supported generator is at a different plant.
5 | `same_location` | `co_located` | A different plant with at least one operating non-storage generator has exactly the same latitude and longitude.
6 | `eia860_flag` | `co_located` | The storage is reported as directly supporting another generator, or as being used to firm co-located renewable generation.
7 | `is_independent` | `standalone` | The storage is reported as independent.
8 | `no_evidence` | `standalone` | None of the preceding rules apply.

Rule 1 is based on the reported energy source rather than the prime mover alone because compressed air storage technologies differ. Older compressed air storage, such as the McIntosh plant in Alabama (currently the only operating compressed air storage plant in the EIA-860 data), uses grid electricity to compress air that is then used to supplement a natural gas combustion turbine. Newer compressed air technologies do not burn fuel and operate more like other standalone or co-located storage. If a `CE` generator reports an energy source of `MWH`, it is categorized using the remaining rules instead, and the pipeline logs a warning so that the assumptions for that technology can be reviewed.

The storage flags used in rules 2, 4, 6, and 7 are reported in the EIA-860 Energy Storage file (accessed through the PUDL `core_eia860__scd_generators_energy_storage` table). A flag is only treated as true if it is explicitly reported as true, since many flags are left blank rather than reported as false. The direct support rule (4) is only applied if the storage is reported as providing direct support, since some storage generators list a supported generator without being reported as providing direct support.

The rules are ordered so that evidence about how storage is metered (rules 1-2) and where it is physically located (rules 3-5) takes precedence over the operational flags reported in EIA-860 (rules 6-7). These can disagree. For example, a battery at a natural gas plant may be reported as independent because it is dispatched independently of the gas turbines, but it is still physically co-located with them, so it is categorized as `co_located`.

A generator is considered to be operating if it is reported as existing in the EIA-860 Generators file for the data year, or if it reported generation or fuel data to EIA-923 or CEMS in the data year. This includes generators that are producing energy while testing before they begin commercial operation, such as a solar array that is still being tested after its co-located battery has entered service. Other retired and proposed generators are not considered when determining whether storage is co-located.

## Assigning Storage Categories to Subplants and Plants

Each subplant that contains at least one storage generator is assigned a `subplant_storage_category`. If a subplant contains multiple storage generators, its `subplant_storage_category_method` is the method of the generator that was categorized by the earliest rule in the list above. Subplants that do not contain a storage generator do not have a storage category.

Each plant that contains at least one storage generator is assigned a `plant_storage_category`, regardless of whether storage is the primary fuel of the plant. For example, a natural gas plant with a co-located battery will have a `plant_storage_category` of `co_located`, even though its `fuel_category` is `natural_gas`.

We expect all of the storage generators at a plant to have the same storage category. The pipeline checks this for every plant, and raises an error if storage generators at the same plant are assigned different storage categories so that the discrepancy can be investigated.

When storage is co-located with a generator at a different plant (identified by rule 4 or rule 5), the `plant_id_eia` of the other plant is recorded for the storage subplant in `co_located_plant_ids`. For example, a battery that is reported as its own plant but is located at the same coordinates as a wind farm will list the `plant_id_eia` of the wind farm. If the storage is co-located with more than one other plant, the IDs are listed as a comma-separated string. This column is blank for all other subplants, including storage that is co-located with a generator at its own plant.

The `plant_storage_category` is included in the plant static attributes table, and the `subplant_storage_category`, `subplant_storage_category_method`, and `co_located_plant_ids` are included in the subplant attributes table.

## Energy storage dispatch data
Because of roundtrip efficiency losses, the cumulative energy discharged by a storage resource over time will always be less than the amount of energy charged. This means that over the course of a month or year, the "net generation" of a storage resource will generally be negative. 

## Known Limitations

The current method has several known limitations:

- **Storage flags are only available in recent years.** The tightly DC-coupled, direct support, and independent flags are only reported in EIA-860 data from 2023 onwards, and the co-located renewable firming flag is only reported from 2016 onwards. For earlier years, storage that is tightly DC-coupled to a solar array cannot be identified as `hybrid`, and will instead usually be categorized as `co_located`. As a result, the storage category of some resources may change between data years even if the resource itself did not change.
- **Pumped storage is not reported in the EIA-860 Energy Storage file.** Pumped storage plants can only be categorized using rules 3, 5, and 8. Many pumped storage plants also include conventional hydroelectric units, so these plants are categorized as `co_located`.
- **Locations are matched exactly.** EIA-860 reports coordinates for each plant rather than each generator, and rule 5 requires an exact match of latitude and longitude. This may miss storage that is located next to another plant but reported with slightly different coordinates, and could incorrectly match plants that share placeholder coordinates.
- **Storage categories do not yet affect other calculations.** Storage categories are currently informational, and are not yet used when assigning fuel categories or hourly profiles.
