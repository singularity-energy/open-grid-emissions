---
stoplight-id: shaping_fleet_data
---

## Shaping EIA-only without plant-specific profiles

One of the primary innovations of the Open Grid Emissions Initiative is its approach to assigning an hourly profile to monthly generation and fuel data reported in EIA-923. This is accomplished using hourly regional generation fleet data reported in EIA-930 (also known as the “Hourly Electric Grid Monitor”). For each regional balancing authority, EIA-930 reports the hourly net generation from all generators of each fuel category (e.g. coal, natural gas, hydro, solar, etc) in that region (which we will refer to as a “fleet”). Since we know the total net generation profile of each fleet, as well as the net generation profile reported to CEMS, we can calculate a residual profile that should theoretically reflect the aggregate profile of all generators in a fleet that do not report to CEMS. This fleet-specific profile can then be used to shape the monthly total data reported to EIA-923. Although this method still has several issues (which will be discussed later), we believe this to be the best currently available method for imputing these shapes, since it is based on observed data.

> Definition: a group of all plants in a balancing authority that have the same fuel category (e.g. coal plants in MISO) are referred to as a "fleet"


## Calculating a residual profile

At its most basic, calculating a residual hourly profile that reflects the generation profile of the part of each fleet that does not report to CEMS involves subtracting the hourly CEMS net generation profile for that fleet from the hourly total EIA-930 profile for that fleet.

To prepare the EIA-930 data for this calculation, it is cleaned and reconciled using a process described elsewhere in this documentation. To prepare the CEMS data for this calculation, first all previously-shaped data (CEMS, partial CEMS subplant, and partial CEMS plant) is added together, and aggregated by balancing authority and fuel category. Each plant may be assigned one of 41 unique primary fuels, but EIA-930 only reports generation totals for 8 broader fuel categories (solar, wind, hydro, nuclear, natural gas, coal, petroleum, and other), each specific energy source type is mapped to one of these broader categories for aggregation.

In some cases, the CEMS profile for a fleet will be larger than the EIA-930 fleet total. This may be a result of several inconsistencies in the way that generation data is reported to EIA-930 versus EIA-860 and EIA-923. These issues are discussed in more depth below, but include instances when generation is reported to a different balancing authority, is categorized under a different fuel category, or not reported to EIA-930 because the plant it is associated with is connected to the distribution grid, rather than the transmission grid.

In these cases when CEMS generation is greater than EIA-930 generation, the residual profile would be negative. Thus, we also calculate a shifted residual profile by shifting the CEMS profile down by a constant amount such that the minimum residual value would be zero. This approach preserves the shape of the residual profile while ensuring no negative values. This approach works well when only small shifts are needed (relative to the magnitude of the EIA-930 profile) to correct the data, but in other cases, this approach would lead to residual profiles that are much larger than the original EIA-930 profile. In such cases the shifted residual profile is not used.

For non-emitting resources like wind and solar, which report data to EIA-930 but not CEMS, the residual profile will always equal the reported EIA-930 value (since the CEMS profile is zero).


## Imputing missing wind and solar profiles

For BAs where wind or solar generation is reported in EIA-923, but not in EIA-930, we use wind or solar data from neighboring BAs to infer the wind or solar profile. Specifically, we use data from all directly-interconnected balancing authorities (other BAs with which the BA exchanges electricity) that are located in the same EIA-930-identified region, and which are located in the same local timezone. This approach minimizes the impact of spatial variation on wind and solar resource profiles, and minimizes any impacts of sunrise/sunset times on variation.

We then average together the wind or solar profiles from these neighboring regions to use as a proxy profile for the wind or solar profile in the region in question. Averaging these profiles without first normalizing them more heavily weights the shape of the profile for regions where there is a greater volume of generation from these sources.

If no data is available from neighboring BAs, we instead use a national average wind or solar profile. This approach is less accurate but only needs to be applied in a few situations. Data quality metrics about the number of BA-months where each method is applied, and validation metrics evaluating how well each of these imputation methods perform when cross-validated with known data are published alongside the dataset.


## Use of other profiles

If neither a good-quality residual profile or shifted residual profile are available for a given fleet, we first fall back to using the total EIA-930 fleet profile to estimate the profile. If a complete EIA-930 fleet profile is not available, we fall back to using the fleet profile for those plants that report to CEMS. Both of these approaches assume that plants of a certain fuel type in a given region generally operate in similar patterns (which may not always be the case). The EIA-930 profile is preferred to the CEMS profile because we assume that the average profile of the entire fleet better represents any single member of the fleet than the profile of a non-random sample of the fleet (CEMS data only represents relatively large generators > 25MW).

If all other attempts at imputing a reasonable hourly profile fail, we assign a flat hourly profile to the data (which is functionally the same as using a monthly average value. For certain fuels, like geothermal, nuclear, biomass, and waste, which tend to run as baseload resources, this may be a reasonable assumption.


## Known issues with using EIA-930 data for hourly profiles

One challenge of comparing data from EIA-860/923 to data from EIA-930 is the uncertainty about whether specific plants are members of the same fleet. EIA-860 and EIA-923 report information about which BA a plant belongs to and what fuel types it consumes, but these might not always be consistent with how balancing authorities report data to EIA-930.

We have identified three primary areas why net generation calculated from EIA-923 and CEMS might not match the net generation reported in EIA-930:



* The generation might be reported under a different balancing authority
* The generation might be reported under a different fuel code
* The generation might not be reported for generators connected to the distribution grid

**Balancing authority assignment**

According to the [EIA-930 instructions](https://www.eia.gov/survey/form/eia_930/instructions.pdf), there are two primary ways that an individual generator can be assigned to a balancing authority:



* The “physical” definition would assign a generator to a balancing authority based on whether that generator is “physically embedded within the tie line boundary of [a] balancing authority.”
* The “commercial” definition would assign a plant to the balancing authority that owns, operates, and/or dispatches that generator.

Our understanding is that the BA reported in EIA-860 represents the plant’s commercial BA, while the plant would be reported in EIA-930 as part of its physical BA.

**Fuel code assignment**

The primary fuel codes assigned to each plant in the OGEI pipeline may not match the generator primary fuel code used by the balancing authority. Thus, it is possible that the set of plants classified as “natural gas”, for example, in EIA-930 may be different from the set of plants that we classify as natural gas.

**Inclusion of distribution-connected plants**

Our understanding is that the data published in EIA-930 only reflects plants have metered telemetry that communicates with the grid operator. In many cases, this may exclude plants that are connected to distribution grids, rather than directly interconnected at the transmission level. Thus, it is possible that the set of plants reported in EIA-930 may not reflect the full set of plants that report to EIA-923 (which includes plants that are connected to distribution grids as well).

## Shaping the Monthly data
Once an hourly profile for each fleet-month has been calculated, it is converted to a percentage of the monthly total value. These percentages are then multiplied by each monthly total value to get the hourly profile for each fleet. This approach ensures that the shaped hourly values, when aggregated back to the monthly level, will equal the total reported monthly value that was shaped.

Currently, we only report hourly data for these EIA-only plants at the fleet level, rather than the individual plant level. We do this for several reasons:
1. Although the hourly CEMS data represents approximately 90% of all electricity-related emissions, the EIA data that is being shaped accounts for approximately 80% of all subplant operating hours in the dataset, meaning that there is a large amount of small plants that need to be shaped. Thus, creating an hourly value for each of these plants would result in a huge dataset that many computers do not have the ability to store in their memory (RAM). Thus, shaping fleet-level data helps keep the size of this dataset managable.
2. Because the method used to shape the data relies on fleet-level observations (rather than plant-specific imputation), we feel that providing a plant-level hourly value may create a sense of false precision in the end result.

Because aggregating the plant-level data to the fleet level removes the `plant_id_eia` identifier, we create a synthetic `shaped_plant_id` to represent these aggregated plants. These synthetic ids are 6-digit identifiers that follow the format `9BBBFF` where `BBB` is the three-digit ba number identified in [this table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/ba_reference.csv), and `FF` represents a two-digit number unique to each fuel category, and defined [here](https://github.com/singularity-energy/open-grid-emissions/blob/afb3ddec0dc93003c21f655b90300c17344107f8/src/impute_hourly_profiles.py#L11). A mapping of `plant_id_eia` to `shaped_plant_id` for each plant can be found in the `plant_static_attributes` table that is included in the dataset.

## Future Work, Known Issues, and Open Questions
- Infer missing hourly profiles for hydro generation ([details](https://github.com/singularity-energy/open-grid-emissions/issues/37)
- Infer hourly profiles for energy storage charge and discharge ([details](https://github.com/singularity-energy/open-grid-emissions/issues/59)
- Should we model hourly shapes for missing peaker or load following genration? ([details](https://github.com/singularity-energy/open-grid-emissions/issues/96)
- Improve imputation of missing wind and solar generation profiles ([details](https://github.com/singularity-energy/open-grid-emissions/issues/171)