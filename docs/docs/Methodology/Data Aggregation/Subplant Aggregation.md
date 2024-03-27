---
stoplight-id: subplant_aggregation
---

## Background on power plant configuration and data crosswalking

> NOTE: There are many different terms used to describe the different parts of a power plant: boiler, generator, turbine (prime mover), and unit. We recommend reading [this resource](https://github.com/USEPA/camd-eia-crosswalk#what-is-the-lowest-level-of-spatial-aggregation-for-each-data-set) from the U.S. EPA about how these terms are defined and how they are related.

CEMS reports data at the unit level (which generally represents a combination of boilers and smokestacks), while EIA-923 data is reported at the generator level (which represents the source of electrical generation itself), or at the boiler level (which for steam plants represents where the fuel is actually combusted). Sometimes units and generators are related in simple one-to-one relationships (i.e. a single boiler powers a single generator, and emits pollution via a single stack), but in other cases these units and generators can be configured in complex arrangements (such as one-to-many, many-to-one, or many-to-many). These relationships are described in the [EPA’s Power Sector Data Crosswalk](https://www.epa.gov/airmarkets/power-sector-data-crosswalk) and in the EIA-860 Boiler-Generator Association table. Whenever comparing data between CEMS and EIA, such as when identifying which data is missing from CEMS, or when calculating gross-to-net generation conversions, it is important to ensure that we are comparing the same data.

## Subplant aggregation methodology

One way to compare data from the two sources would be to aggregate the data to the plant level, which includes all units and generators. However, at the plant level we lose important details about the operation patterns and emissions from each generator at a plant.

Thus, using a methodology developed as part of the PUDL project, we aggregate the data into “subplants,” which represent distinct, independent groupings of units and generators at a plant. In some cases, a subplant might consist of a single unit and single generator, but in other cases, it might consist of multiple generators and units that are interconnected (such as for combined-cycle generators). For example, in figure 1 of the [above referenced]((https://github.com/USEPA/camd-eia-crosswalk#what-is-the-lowest-level-of-spatial-aggregation-for-each-data-set)) EPA explainer, each of the four examples (A, B, C, D) would be considered its own subplant.

The PUDL module for identifying subplants uses graph analysis to identify these unique groupings. As explained in [their documentation](https://catalystcoop-pudl.readthedocs.io/en/latest/autoapi/pudl/analysis/epa_crosswalk/index.html):

> In graph analysis terminology, the [Power Sector Data] crosswalk is a list of edges between nodes (combustors and generators) in a bipartite graph. The `networkx` python package provides functions to analyze this edge list and extract disjoint subgraphs (groups of combustors and generators that are connected to each other). These are the distinct power plants. To avoid a name collision with plant_id, we term these collections ‘subplants’, and identify them with a subplant_id that is unique within each plant_id. Subplants are thus identified with the composite key (plant_id, subplant_id).

Because the pudl methodology only relies on the power sector data crosswalk to identify these groupings and only identifies subplants based on unit to generator relationships (ignoring how generators may be connected to boilers), a subplant_id is only generated for those units/generators that exist in the crosswalk. Thus, the `subplant_id`s created by pudl may be incomplete or incorrect. Because the `subplant_id` is used as one of the primary keys for crosswalking CEMS and EIA data, we want to ensure that there is a subplant_id assigned for every generator in the EIA data and every unit in the CEMS data. In the future, we plan to update the pudl epa_crosswalk code to enable the creation of a complete set of `subplant_id` (see [this issue](https://github.com/singularity-energy/open-grid-emissions/issues/49) for updates). However, in the meantime, we have implemented a methodology to expand and correct the subplant_id mapping in our own data pipeline.

Because the current subplant_id code does not take boiler-generator associations into account, there may be instances where the code assigns generators to different subplants when in fact, according to the boiler-generator association table, these generators are grouped into a single unit based on their boiler associations. PUDL uses the same network analysis method that it uses to identify subplants to assign a unique `unit_id_pudl` to each boiler-generator grouping. The `unit_id_pudl` is similar to a `subplant_id` except that `unit_id_pudl` is based on boiler-generator groups, while `subplant_id` is based on unit-generator groups. Thus, our goal is to harmonize these `unit_id_pudl` with `subplant_id` such that a single `subplant_id` is assigned to any generators with overlapping subplant_id or unit_id_pudl.

At a high level, we update and expand the `subplant_id` mapping based on the following steps:
1. Use the PUDL `subplant_id` if available. In the case where a `unit_id_pudl` groups several subplants, we overwrite these multiple existing `subplant_id` with a single `subplant_id`.
2. Where there is no PUDL `subplant_id`, we use the `unit_id_pudl` to assign a unique `subplant_id`
3. Where there is neither a pudl `subplant_id` nor `unit_id_pudl`, we use the generator ID to assign a unique `subplant_id`, assuming that this generator represents its own subplant.
4. All of the new unique ids are renumbered in consecutive ascending order


> Note: `subplant_id` should be stable across all data years within a single version of OGE. However, `subplant_id` may change from one version of OGE to the next as we get better data about subplant mappings. 