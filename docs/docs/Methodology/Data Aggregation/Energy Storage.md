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

Storage category | Description
---------|---------
`standalone` | The storage resource is an independent plant and charges from the grid.
`co_located` | A storage resource located at the same plant as one or more generators. The storage may charge from the co-located generator and/or the grid.
`hybrid` | The storage is an integral part of a generator is treated like a single resource.

## Assigning a Storage Category to Each Storage Generator

Each storage generator is assigned the storage category of the first of the following rules that applies. The rule that was used for each subplant is reported in the `subplant_storage_category_method` column of the subplant attributes table.

\# | `storage_category_method` | Storage category | Rule
---------|---------|---------|---------
1 | `hybrid_prime_mover` | `hybrid` | The generator's prime mover is `CE` (compressed air), and it reports an energy source other than `MWH`. This indicates that the storage is integrated into a generator that uses another energy source, such as compressed air storage that burns natural gas when discharging. This rule was designed to reflect the McIntosh CAES generator in Alabama.
2 | `pumped_storage_with_inflow` | `hybrid` | The generator is pumped storage (`PS`), and the plant reports discharging at least as much energy as it charges over the year (and some discharge) in the EIA-923 Energy Storage file. This indicates that the plant also generates electricity from natural inflow to its reservoir.
3 | `dc_coupled_tightly` | `hybrid` | The storage is reported as tightly DC-coupled in the EIA-860 Energy Storage file, indicating that it cannot charge from the grid.
4 | `same_plant` | `co_located` | The storage resource shares a `plant_id_eia` with at least one non-storage generator.
5 | `direct_support_other_plant` | `co_located` | The storage is reported as directly supporting a generator in the EIA-860 Energy Storage file, and the supported generator has a different `plant_id_eia`.
6 | `same_location` | `co_located` | The storage resource has exactly the same latitude and longitude as a non-storage generator.
7 | `eia860_flag` | `co_located` | The EIA-860 Energy storage file reports that the resource directly supports another generator, or is being used to firm co-located renewable generation.
8 | `is_independent` | `standalone` | The storage is reported as independent in the EIA-860 Energy Storage file.
9 | `no_evidence` | `standalone` | None of the preceding rules apply.

Rule 1 is based on the reported energy source rather than the prime mover alone because compressed air storage technologies differ. Older compressed air storage, such as the McIntosh plant in Alabama (currently the only operating compressed air storage plant in the EIA-860 data), uses grid electricity to compress air that is then used to supplement a natural gas combustion turbine. Newer compressed air technologies do not burn fuel and operate more like other standalone or co-located storage. If a `CE` generator reports an energy source of `MWH`, it is categorized using the remaining rules instead, and the pipeline logs a warning so that the assumptions for that technology can be reviewed.

Rule 2 identifies pumped storage plants that both pump water to store energy and generate electricity from water that flows naturally into their upper reservoir. Pumped storage that only discharges energy that it previously pumped will always discharge less energy than it charges, due to round-trip efficiency losses, so discharging at least as much energy as is charged indicates that some of the generation comes from natural inflow. These plants keep a `fuel_category` of `storage`, but their discharge data is left blank (see below).

The storage flags used in rules 3, 5, 7, and 8 are reported in the EIA-860 Energy Storage file (accessed through the PUDL `core_eia860__scd_generators_energy_storage` table). A flag is only treated as true if it is explicitly reported as true, since many flags are left blank rather than reported as false. The direct support rule (5) is only applied if the storage is reported as providing direct support, since some storage generators list a supported generator without being reported as providing direct support.

The rules are ordered so that evidence about how storage is metered or operates (rules 1-3) and where it is physically located (rules 4-6) takes precedence over the operational flags reported in EIA-860 (rules 7-8). These can disagree. For example, a battery at a natural gas plant may be reported as independent because it is dispatched independently of the gas turbines, but it is still physically co-located with them, so it is categorized as `co_located`.

A generator is considered to be operating if it is reported as existing in the EIA-860 Generators file for the data year, or if it reported generation or fuel data to EIA-923 or CEMS in the data year. This includes generators that are producing energy while testing before they begin commercial operation, such as a solar array that is still being tested after its co-located battery has entered service. Other retired and proposed generators are not considered when determining whether storage is co-located.

## Assigning Storage Categories to Subplants and Plants

Each subplant that contains at least one storage generator is assigned a `subplant_storage_category`. Subplants that do not contain a storage generator do not have a storage category. Each plant that contains at least one storage generator is assigned a `plant_storage_category`, regardless of whether storage is the primary fuel of the plant. We expect all of the storage generators at a plant to have the same storage category. The pipeline checks this for every plant, and raises an error if storage generators at the same plant are assigned different storage categories so that the discrepancy can be investigated.

When storage is co-located with a generator at a different plant (identified by rule 5 or rule 6), the `plant_id_eia` of the other plant is recorded for the storage subplant in `co_located_plant_ids`. If the storage is co-located with more than one other plant, the IDs are listed as a comma-separated string. This column is blank for all other subplants, including storage that is co-located with a generator at its own plant.

## Energy Storage Charging and Discharging Data

Because of roundtrip efficiency losses, the cumulative energy discharged by a storage resource over time will always be less than the amount of energy charged. This means that over the course of a month or year, the "net generation" of a storage resource will generally be negative. The `net_generation_mwh` of a storage resource therefore does not show how much energy it charged and discharged, so the monthly and annual plant and power sector results also include:

Column | Description
---------|---------
`storage_charge_mwh` | The electricity used to charge energy storage resources, in MWh
`storage_discharge_mwh` | The electricity discharged by energy storage resources, in MWh

These columns are only reported for energy storage resources. They are blank for all other resources, and are zero for storage resources that did not charge or discharge in a given month. The `net_generation_mwh` column is unchanged, and for most storage resources is equal to `storage_discharge_mwh` minus `storage_charge_mwh`. These columns are not yet included in the hourly results.

### Source data

Monthly charging and discharging data is reported for each plant, prime mover, and energy source in the EIA-923 Energy Storage file (accessed through the PUDL `core_eia923__monthly_energy_storage` table). 
- For storage that reports its energy input in MWh, `storage_charge_mwh` is the reported fuel consumed in physical units (MWh in this case) for electric generation, and `storage_discharge_mwh` is the reported gross generation.
- For compressed air storage that burns natural gas when discharging (such as McIntosh), the fuel consumed value represents the natural gas combusted during discharge rather than the electricity used to compress air. For these resources, `storage_discharge_mwh` is the reported gross generation of the combustion turbine, and `storage_charge_mwh` is the gross generation minus the net generation. Because these resources keep the fuel category of their fuel (see [Plant Primary Fuel](Plant%20Primary%20Fuel.md)), their charging and discharging data is included in the `natural_gas` fleet in the power sector results.
- For pumped storage with natural inflow (rule 2 above), `storage_charge_mwh` is the reported energy used for pumping, but `storage_discharge_mwh` is left blank, since the discharge of previously pumped water cannot be separated from generation from natural inflow.

Energy storage data is reported for each plant and prime mover, but not for individual generators. When a plant has multiple storage generators with the same prime mover, the charging and discharging data is allocated to each generator based on its share of the nameplate capacity. This matches how net generation is allocated to these generators, since EIA-923 does not report generator-level data for them.

Charging and discharging data is taken directly from EIA-923 for all storage resources, including storage located at plants that report to CEMS. It is not shaped using the hourly profiles of other resources at the plant, since the operation of storage is unlikely to follow the operation of combustion generators. As a result, the charging and discharging data for a storage resource at a plant that reports to CEMS may not exactly reconcile with the net generation of that resource, which may be calculated from CEMS data.

### Data quality

The pipeline checks whether `net_generation_mwh` equals `storage_discharge_mwh` minus `storage_charge_mwh` for each storage resource and month. If a month reports zero discharge, but filling the discharge with the net generation plus the charge would make that month consistent, and all other months for that resource are consistent, the discharge is filled with this value. This addresses what appear to be missing values, such as a pumped storage plant that reports charging and net generation in a month, but no discharge.

Other inconsistencies are reported in the pipeline log, but the reported values are not changed. These include:

- **Consistent differences in every month.** Some plants report a small, consistent difference between net generation and discharge minus charge every month. This likely reflects gross generation being reported before station service or auxiliary loads are subtracted.
- **Annual discharge greater than annual charge.** A storage resource cannot discharge more energy than it charges. This is reported for some solar-plus-storage plants, where the battery may charge partly from the on-site solar array in a way that is not metered, and for plants in their first months of operation.

Because the discharge of pumped storage with natural inflow is left blank, and because compressed air storage is included in the `natural_gas` fleet, `storage_discharge_mwh` minus `storage_charge_mwh` will not equal `net_generation_mwh` when aggregated to the fleet or balancing authority level.

As with other EIA-923 data, monthly charging and discharging data for plants that only report annually to EIA is estimated by EIA. The share of energy storage data that comes from annual reporters is included in the `annually_reported_eia_data` data quality metrics table. The EIA-923 Energy Storage file is not available for the earliest years of data, so these columns are blank for those years.

## Known Limitations

The current method has several known limitations:

- **Storage flags are only available in recent years.** The tightly DC-coupled, direct support, and independent flags are only reported in EIA-860 data from 2023 onwards, and the co-located renewable firming flag is only reported from 2016 onwards. For earlier years, storage that is tightly DC-coupled to a solar array cannot be identified as `hybrid`, and will instead usually be categorized as `co_located`. As a result, the storage category of some resources may change between data years even if the resource itself did not change.
- **Pumped storage is not reported in the EIA-860 Energy Storage file.** Pumped storage plants can only be categorized using rules 2, 4, 6, and 9. Pumped storage plants that also include conventional hydroelectric units are categorized as `co_located` unless rule 2 applies.
- **Pumped storage with natural inflow is identified from annual totals.** A pumped storage plant that reports exactly equal charging and discharging (for example, due to how the data was reported) will be categorized as `hybrid` and have its discharge left blank.
- **Locations are matched exactly.** EIA-860 reports coordinates for each plant rather than each generator, and rule 5 requires an exact match of latitude and longitude. This may miss storage that is located next to another plant but reported with slightly different coordinates, and could incorrectly match plants that share placeholder coordinates.
