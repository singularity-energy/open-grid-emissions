---
stoplight-id: ba_aggregation
---

## Aggregating plant data to balancing authority

In addition to plant-level data, OGE includes data aggregated to the balancing authority level.

As explained in the [eGRID technical support documentation](https://www.epa.gov/system/files/documents/2022-01/egrid2020_technical_guide.pdf):


> A balancing authority is a portion of an integrated power grid for which a single dispatcher has operational control of all electric generators. A balancing authority is the responsible entity that integrates resource plans ahead of time, maintains demand and resource balance within a BA area, and supports interconnection frequency in real time. The balancing authority dispatches generators in order to meet an area’s needs and can also control load to maintain the load-generation balance.

Balancing authorities are assigned to each plant based on plant-level data reported in EIA-860. If a plant is missing a reported BA code, we attempt to infer the ba based on the reported balancing authority name, the reported utility name, and the reported transmission or distribution system owner name, in that order. This mapping is based on [this table](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/oge/reference_tables/utility_name_ba_code_map.csv) which associates the full balancing authority name with the ba_code.

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

## EIA-930A generator inventories

Since 2023, each BA that reports to EIA-930 also files [Form EIA-930A](https://www.eia.gov/survey/#eia-930), an annual inventory of the generators reflected in its EIA-930 data, matched by EIA to EIA-860 plant and generator IDs. It shows which BA's EIA-930 net generation includes each generator, which is the boundary that OGE's BA-level data needs to match to be consistent with EIA-930. OGE downloads the published files and cleans three of their schedules:
- **Schedule 2:** the generators included in the BA's EIA-930 net generation at the end of the year.
- **Schedule 3:** the generators that the BA plans to include in the following year.
- **Schedule 4:** the generators with which the BA has a pseudo-tie or dynamic scheduling arrangement, the partner BAs, and whether the BA includes the generator in its Schedule 2.

In each schedule, plant and generator IDs are the EIA-860 IDs where EIA matched the generator to EIA-860 (`validated_eia860`), and the IDs submitted by the BA otherwise. Generator IDs are upper case. Schedule 4 partners are converted to EIA-930 BA codes (e.g. "CAISO" to "CISO", "BPA" to "BPAT"), and values that aren't a BA are kept in separate columns (`arrangement_type` for "Dynamic import" and "Pseudo tie", `partner_note` for e.g. "WAPA", which operates several BAs).

The cleaned schedules are written to `outputs/{year}/eia930/eia930a_schedule_{2,3,4}.parquet`. EIA-930A is published about a year after the end of the year, so it is often not available when a year is first run. In that case, the latest published year of EIA-930A is used (and 2023, the first year, for earlier years), and the year used is in the `report_year` column. The published years are set in `constants.py`. EIA-930A is not yet used to assign plants to BAs.

The notebook `notebooks/validation/validate_ba_assignment.ipynb` uses EIA-930A to check OGE's assignment of plants to BAs. It compares OGE net generation by BA and fuel with EIA-930, for all generation and for the generators listed in Schedule 2, and flags plants that are far from every other plant in their BA and much closer to another BA's plant, on a map.
