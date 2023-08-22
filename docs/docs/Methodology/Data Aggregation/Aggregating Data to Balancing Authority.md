---
stoplight-id: ba_aggregation
---

## Aggregating plant data to balancing authority

In addition to plant-level data, OGEI includes data aggregated to the balancing authority level.

As explained in the [eGRID technical support documentation](https://www.epa.gov/system/files/documents/2022-01/egrid2020_technical_guide.pdf):


> A balancing authority is a portion of an integrated power grid for which a single dispatcher has operational control of all electric generators. A balancing authority is the responsible entity that integrates resource plans ahead of time, maintains demand and resource balance within a BA area, and supports interconnection frequency in real time. The balancing authority dispatches generators in order to meet an area’s needs and can also control load to maintain the load-generation balance.

Balancing authorities are assigned to each plant based on plant-level data reported in EIA-860. If a plant is missing a reported BA code, we attempt to infer the ba based on the reported balancing authority name, the reported utility name, and the reported transmission or distribution system owner name, in that order. This mapping is based on [this table](https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/utility_name_ba_code_map.csv) which associates the full balancing authority name with the ba_code.

For plants that do not belong to a specific BA, we assign a miscellaneous BA code based on the two-letter state code of the state where the plant is located. For example, many plants in Alaska will be assigned the code “AKMS” and many plants in Hawaii will be assigned a code of “HIMS.” For any plants that are associated with a balancing authority that is no longer active, we also replace the retired BA code with a state-based miscellaneous code.

## Commercial vs physical balancing authorities

According to the [EIA-930 instructions](https://www.eia.gov/survey/form/eia_930/instructions.pdf), there are two primary ways that an individual generator can be assigned to a balancing authority:
- The “physical” definition would assign a generator to a balancing authority based on whether that generator is “physically embedded within the tie line boundary of [a] balancing authority.”
- The “commercial” definition would assign a plant to the balancing authority that owns, operates, and/or dispatches that generator.

In general, it seems that the “commercial” definition of balancing authority is used, whereas EIA-930 “is attempting to represent electric system operations in as purely a physical way as possible. Ownership and dispatch are irrelevant to the determination of what is associated with a balancing authority.”

More specifically, the EIA-930 instructions state the following requirements for reporting data:

> Generators physically embedded within the tie line boundary of your balancing authority, but owned, operated, or dispatched by another balancing authority, are to be included, for the purposes of the EIA-930, in your reporting of net generation. The transmission connection of that plant to your system is not considered to be a tie line boundary.

In EIA Form 860, information is collected both about a plant’s balancing authority, and about who owns the transmission lines that a plant is connected to. The [instructions](https://www.eia.gov/survey/form/eia_860/instructions.pdf) for Form 860 note that “A balancing authority manages supply, demand, and interchanges within an electrically defined area. It may or may not be the same as the Owner of Transmission/Distribution Facilities.”

It would thus seem that the “balancing authority” reported in Form 860 would represent the plant’s “commercial” balancing authority, while the “transmission owner” would represent the plant’s “physical” balancing authority.

Thus, in the OGE dataset, all references to `ba_code` represent commercial balancing authorities. However, we also assign each plant a `ba_code_physical` based on the reported “Transmission or Distribution System Owner” in Form 860, Schedule 2. However, `ba_code_physical` is not yet used in the data pipeline until we can better understand how it relates to the EIA-930 data.
