---
stoplight-id: egrid_comparison
---

The Open Grid Emissions (OGE) methodology is based on the U.S. EPA's [eGRID](https://www.epa.gov/egrid) methodology, and both datasets are built mostly from the same public data (EPA CAMD/CEMS, EIA-860, and EIA-923). This page compares the two in three parts:

1. **[What OGE provides that eGRID does not](#1-new-features-and-data-added-in-oge)**
2. **[Where both datasets report the same outputs, how and why the values differ](#2-how-does-oge-differ-from-egrid)**
3. **[What eGRID provides that OGE does not (yet)](#3-what-egrid-provides-that-oge-does-not-yet)**

**Versions compared.** All results are for data year 2023.

- **eGRID:** eGRID2023 (published June 2025), including the [eGRID2023 Technical Guide](https://www.epa.gov/egrid) and the [eGRID R production model](https://github.com/USEPA/egrid) at commit [`08eddc7`](https://github.com/USEPA/egrid/commit/08eddc7) (June 12, 2025).
- **OGE:** data release v0.8.0, with the methodology as of OGE commit [`b1ffdb0`](https://github.com/singularity-energy/open-grid-emissions/commit/b1ffdb0) (v0.8.0 plus later changes).

## 1. New features and data added in OGE

- **Hourly and monthly data.** eGRID publishes annual values only. OGE publishes hourly (2019 onward), monthly, and annual data for plants, subplants, and balancing areas, with timestamps in both UTC and local time.
- **A longer, consistent time series.** OGE covers 2005 to the present. Every release reprocesses all years with the current methodology, so years can be compared directly.
- **Consumption-based emission rates.** OGE calculates hourly, monthly, and annual *consumed* emission rates for each balancing area, using EIA-930 interchange data and a multi-region input-output model ([methodology](../Methodology/Emissions%20Calculations/Consumption-based%20Emissions.md)). eGRID only reports rates for generation within each region.
- **Subplant-level data.** OGE links CEMS units, EIA boilers, and EIA generators into "subplants" ([methodology](../Methodology/Data%20Aggregation/Subplant%20Aggregation.md)) and publishes data at that level. This sits between eGRID's unit/generator files and its plant file.
- **Every adjustment variant at every level.** For each pollutant OGE publishes unadjusted, CHP-adjusted (`_for_electricity`), biomass-adjusted (`_adjusted`), and both-adjusted values, for plants, subplants, and balancing areas. eGRID publishes unadjusted values only at the unit and plant level; its aggregated files are adjusted only.
- **Energy storage detail.** OGE has a separate storage fuel category, storage charging and discharging (MWh), and a classification of each storage resource as standalone, co-located, or hybrid ([methodology](../Methodology/Data%20Aggregation/Energy%20Storage.md)).
- **Data quality metadata.** For each subplant-month, OGE records the input data source, hourly profile method, and gross-to-net method. It also publishes summary metrics such as the share of data from CEMS vs EIA, the share of EIA data from annual reporters, and CEMS measurement quality ([details](../Data%20Validation/Data%20Quality%20Metrics.md)).

## 2. How does OGE differ from eGRID?

### Summary of the 2023 comparison

Plants were matched on the eGRID plant ID (12,540 matched plants). Values are unadjusted unless noted.

| Metric | eGRID2023 | OGE 2023 | OGE vs eGRID | Plants within 1% |
|---|---|---|---|---|
| Net generation | 4,173 TWh | 4,170 TWh | −0.1% | 98.8% |
| Heat input | 35.24 B MMBtu | 35.24 B MMBtu | 0.0% | 98.2% |
| CO2 | 1,829 Mt | 1,857 Mt | +1.5% | 94.1% |
| NOx | 1.162 Mt | 1.252 Mt | +7.7% | 90.6% |
| SO2 | 0.908 Mt | 0.934 Mt | +2.9% | 81.2% |
| CH4 / N2O | — | — | +1.2% / +1.2% | 85.7% / 82.0% |
| CO2, biomass + CHP adjusted | 1,593 Mt | 1,601 Mt | +0.5% | 92.2% |
| NOx, biomass + CHP adjusted | 0.913 Mt | 1.038 Mt | +13.6% | 89.8% |

eGRID also includes 62 plants in Puerto Rico (17.5 TWh, 13.5 Mt CO2), which OGE does not yet cover (see section 3).

In most cases, the OGE and eGRID data are identical. However, where the outputs differ, they are often driven by bug fixes and methodological enhancements introduced in the OGE dataset. Specific differences are described in the tables below, which fall into four groups:

- **[2a. Convention differences](#2b-convention-differences):** choices where neither approach is more accurate, such as IDs, fuel categories, and BA labels.
- **[2b. eGRID2023 bugs](#2c-egrid2023-bugs):** places where the eGRID code doesn't do what the Technical Guide describes or intends.
- **[2c. OGE methodological enhancements](#2a-oge-methodological-enhancements):** places where eGRID's method is deliberate and OGE's is more accurate or complete.
- **[2d. Areas where eGRID is currently more accurate, and OGE fixes](#2d-areas-where-egrid-is-currently-more-accurate-and-oge-fixes).**

In the tables below, **Effect** describes the implication of the methodological difference. **Observed impact** describes the overall impact of this difference on outputs between eGRID2023 and OGE 2023 data. **Example** lists a specific example of this difference.

### 2a. Convention differences

These differ between the datasets, but neither approach is more accurate.

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| Plant IDs | Some EIA plant IDs are changed to EPA IDs, and some plants are combined (eGRID Table C-5) | EIA plant IDs (`plant_id_eia`) | No effect on totals. | 23 plants matched through the crosswalk; 23 OGE plants (14.0 TWh, 9.3 Mt CO2) have no eGRID match, mostly because of 2023 ID changes or different exclusion lists | Indiana Harbor West + E 5 AC Station (10397 + 54995 in OGE): reported as 10474 in eGRID (4.4 Mt CO2) |
| Gross to net generation conversion limits | Uses raw EIA-923 net generation, even if gross-to-net ratio is anomalous | Uses EIA-923 net generation unless ratio is anomalous; backstops to fleet-average, EIA default. Implicitly trusts gross generation data from CEMS to ensure consistency with reported CEMS emissions. | Net generation may differ for some plants | 14 plants: 13.4 vs 12.7 TWh (OGE −4.8%); median plant difference 37%. Other CEMS plants: 96.9% within 0.5%. | Dearborn Industrial Generation (55088): 5.26 vs 3.30 TWh (OGE −37%) |
| Regional fuel totals | Each plant's generation by fuel comes from the net generation EIA-923 reports for each fuel (Generation and Fuel data), and is summed to regions. A plant that burns several fuels contributes to several categories. | BA fuel totals assign each subplant's entire generation to its primary fuel category | OGE reflects fuel totals based on generator fleet aggregation, which better aligns with how grid operators report this data; eGRID reports based on electricity generation from combustion of each fuel | Relative to EIA-923 fuel-level generation (eGRID's method) for the same plants: natural gas −6.7 TWh, petroleum −1.5 TWh, biomass +5.4 TWh, coal +3.4 TWh; other categories within 0.5 TWh. About 8.8 TWh (0.2% of US generation) moves between categories. | Mansfield Mill (54091): biomass 0.37 vs 0.60 TWh, natural gas 0.49 vs 0.28 TWh.<br>Edwardsport (1004), a coal gasification combined cycle plant: coal 2.09 vs 3.27 TWh, natural gas 1.18 vs 0 TWh. |
| Fuel categories | **Categories:**<br>gas<br>oil<br>other fossil<br>(none)<br>(none)<br><br>**Category mapping:** see table below | **Categories:**<br>natural_gas<br>petroleum<br>(none)<br>waste<br>storage<br><br>**Category mapping:** see table below | Convention difference; totals are unaffected. OGE's [`energy_source_groups.csv`](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/oge/reference_tables/energy_source_groups.csv) includes an eGRID category column for comparisons. | Generation classified differently: 11.8 TWh (BFG/OG/PRG/SG), 7.4 TWh (TDF/MSN), 5.7 TWh (MSW/MSB) | US Steel Gary Works (50733): 0.64 TWh "other fossil" in eGRID, counted as natural_gas in OGE |

**Fuel category mapping.** Energy source codes that are categorized differently; all others map to equivalent categories.

| Energy source | eGRID category | OGE category |
|---|---|---|
| Propane Gas (PG) | gas | petroleum |
| Butane Gas (BU) | gas | petroleum |
| Blast Furnace Gas (BFG) | other fossil | natural_gas |
| Other Gas (OG) | other fossil | natural_gas |
| Process Gas (PRG) | other | natural_gas |
| Other Syngas (SG) | not categorized | natural_gas |
| Tire-derived Fuels (TDF) | other fossil | waste |
| Municipal Solid Waste, nonbiogenic (MSN) | other fossil | waste |
| Municipal Solid Waste (MSW) | biomass | waste |
| Municipal Solid Waste, biogenic (MSB) | biomass | waste |
| Energy Storage (MWH) | other | storage |
| Refinery Gas (RG) | not categorized | petroleum |
| Coal Synfuel (SC) | not categorized | coal |

### 2b. eGRID2023 bugs

These bugs represent differences between the methodology described in the eGRID Technical Guide and the implemented [eGRID R production model](https://github.com/USEPA/egrid). 

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| Landfill gas NOx flaring adjustment: Mismatched units | When adjusting NOx for the flaring baseline, eGRID subtracts biomass NOx in **lb** from unadjusted NOx in **short tons**. This results in an adjusted NOx value of zero at every plant that burns any landfill gas. | Calculation uses consistent units | OGE is more accurate | 291 plants. Adjusted NOx 0 vs 88.2 kt, about 70% of the national gap in adjusted NOx. Unadjusted landfill gas NOx is also higher in OGE (55.9 vs 91.9 kt), largely because of different emission factors. | Curtis H. Stanton Energy Center (564), a coal plant that co-fires landfill gas: adjusted NOx 0 vs 3,522 t |
| Landfill gas non-CO2 adjustments: documentation vs code | The Technical Guide (section 2.2) says the eGRID2023 landfill gas adjustment "was updated to remove adjustments made to NOx, SO2, CH4, and N2O," leaving CO2 only. The code still adjusts all four (`plant_file_create.R`): NOx and SO2 by subtracting flaring baselines (0.08 and 0.0115 lb/MMBtu), and CH4 and N2O by subtracting landfill gas emissions. | NOx reduced by a flaring baseline; SO2, CH4, and N2O from landfill gas set to zero | eGRID's adjusted NOx, SO2, CH4, and N2O for landfill gas plants don't follow its documented method, so users relying on the guide will misread them. OGE's approach is closer to what the eGRID code does than to the guide. | 291 plants. Adjusted SO2 3.0 kt (eGRID) vs 1.6 kt (OGE). For NOx, see the unit-conversion bug above. | — |
| Missing CO2 and SO2 data for ozone-season CEMS reporters | For units that only report to CEMS during the ozone season (May–Sep), eGRID fills in their Oct–Apr fuel consumption and NOx from EIA-923 (see 2c), but never fills in CO2 and SO2 values for these months, meaning that these units are missing 7 months of CO2 and SO2 data. | Full-year CO2 and SO2: Oct–Apr from EIA-923 fuel × emission factors, May–Sep from CEMS | OGE is more accurate: eGRID understates annual CO2 and SO2 for these units | At the 80 plants with ozone-season reporters: CO2 61.2 vs 75.5 Mt (OGE +23%). Nearly all of the gap is at plants that also have zero-reported CO2 (2c); the other 22 plants differ by only 0.1 Mt. 47 of the 80 differ by more than 1%. | Eastman Chemical Company (50481): 9 ozone-season units; CO2 0 vs 2.43 Mt |
| Ozone-season CEMS reporters at plants with year-round units: ozone-season heat input subtracted twice | At plants that also have year-round CEMS units, eGRID fills in Oct–Apr heat input using the EIA-923 heat input left over after subtracting CEMS data. The code (`epa_q_oz_reporters_dist` in `unit_file_create.R`) subtracts the CEMS total, which already includes the ozone-season units' May–Sep heat input, and then also subtracts the EIA-923 May–Sep heat input for that plant and prime mover. Ozone-season heat input is removed twice, so the leftover is understated. When it's zero or negative, the unit isn't gap-filled. | No equivalent step (Oct–Apr months use each subplant's own EIA-923 data) | OGE is more complete, but the effect is small | 16 plant/prime-mover groups have ozone-season units alongside year-round units; eGRID gap-fills 3 of them. Removing the double subtraction would gap-fill 1–4 more small gas turbine groups (about 1.6k–23k MMBtu, depending on how the leftover is defined). The six coal steam groups would still have nothing to distribute, because their CEMS heat input already exceeds EIA-923 fuel use. | Kendall Green Energy (1595): leftover of −2,418 MMBtu, so the ozone-season turbine (S6) isn't gap-filled; without the double subtraction the leftover is +1,613 MMBtu |
| Geothermal type table | The NREL table eGRID uses codes every non-binary plant as dry steam ("ST") and has no flash ("F") entries, although the emission factor table has flash factors. Flash plants therefore get dry-steam factors: CO2 26.0 vs 17.6 lb/MMBtu, SO2 0.00006 vs 0.103 lb/MMBtu (`geothermal_emission_factors.csv`). | Flash, dry steam (Geysers only), and binary types | OGE is more accurate: eGRID overstates CO2 and nearly eliminates SO2 at flash plants | 67 plants: CO2 499 vs 402 kt (OGE −19%); SO2 about 1 vs 924 t | Dixie Valley (52015), a flash plant coded as steam: CO2 20.9 vs 13.2 kt (OGE −37%); SO2 0.05 vs 77 t |


### 2c. OGE methodological enhancements

These are deliberate methodological choices in eGRID that OGE improves on.

#### CEMS data processing

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| CEMS data resolution | Annual totals from CAMD daily/quarterly data | Hourly CEMS data | OGE is more accurate: hourly data allows anomaly screening and hour-level gap filling | Unable to determine isolated impact | — |
| Zero-reported CO2 in CAMD data | eGRID's CO2 estimation step fills only missing (NA) values (`rows_patch` in `unit_file_create.R`), so units that report heat input with zero CO2 emissions keep zero CO2. This mostly affects units in NOx-only reporting programs (SIP NOx, CSAPR), which often report zero (rather than missing) CO2. | CO2 is treated as missing when it's zero but fuel use is positive, and is estimated from fuel-specific combustion emission factors | OGE is more accurate: these units' CO2 is missing from eGRID | 363 units at 140 plants; 379 M MMBtu of heat input; about 25 Mt of CO2 not counted in eGRID. At these plants OGE's unadjusted CO2 is 31.4 Mt higher and adjusted CO2 13.7 Mt higher. | US Steel Corp – Gary Works (50733): CO2 0 vs 5.07 Mt |
| Other missing CEMS values | Missing CO2 filled with the primary-fuel emission factor | Missing CO2 filled with fuel-specific emission factors; months where CEMS shows zero generation and fuel are filled from EIA-923 | Similar; OGE is slightly more complete | Unable to determine isolated impact | — |
| Units that report to CEMS only during the ozone season (May–Sep): gap-fill method | Oct–Apr heat input comes from EIA-923 plant/prime-mover heat input, split among the ozone-season units by their share of ozone-season heat input; | Each month without CEMS data uses that subplant's allocated EIA-923 fuel, generation, and emissions for all pollutants; May–Sep uses CEMS | OGE is more accurate and complete: all 12 months and all pollutants, allocated by generator rather than by ozone-season share | 80 plants: 86 gap-filled units and 70 never gap-filled. Heat input (−2%) and NOx (−1%) are similar in total. | Covanta Niagara (50472): eGRID assigns a natural gas boiler (BLR05) with 22k MMBtu of May–Sep CEMS heat input another 1.97M MMBtu for Oct–Apr, from a plant/prime-mover pool that includes fuel from the plant's non-CEMS municipal solid waste boilers. OGE doesn't include this boiler. |
| Steam-only CEMS units (units reporting heat input and steam output but no electricity) | Included in plant totals; only scaled down if the plant is flagged as CHP | Excluded if the unit isn't linked to an EIA generator | OGE is more accurate at industrial and district-steam facilities, where these boilers don't make electricity. eGRID may be more complete for startup/auxiliary boilers at power plants. | 19 plants, 46 units. eGRID includes 54.3 M MMBtu and 2,405 t NOx from 36 of them; at these plants OGE heat input is 3% lower and NOx 4% lower. | Purdue University–Wade Utility (50240): heat input 3.38 M vs 2.63 M MMBtu (OGE −22%) |
| Other CEMS unit filters | Excludes units listed as retired, future, or in long-term cold storage | Keeps any unit that reports data in the year | OGE is slightly more complete | Small: 8 CEMS units (7 M MMBtu) at eGRID plants are not in the eGRID unit file | — |

#### EIA-923 data allocation

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| Allocating EIA-923 generation and fuel to generators | **Generation:** from EIA-923 generator-level data where reported, otherwise plant/prime-mover generation distributed by nameplate capacity.<br><br>**Fuel:** not allocated to generators by fuel type. Each unit gets a single total heat input, summed across all fuels (EIA-923 boiler fuel for boilers; plant/prime-mover heat input distributed by nameplate capacity for other generators), and a single primary fuel. | **Generation and fuel:** allocated to each generator separately for each fuel it burns (one record per generator per energy source code per month), in proportion to reported generator- and boiler-level data, otherwise by nameplate capacity. | Plant totals are the same by design. OGE is more accurate at the generator level, and keeping each fuel separate is what allows fuel-specific emissions (see Fuel resolution below). | EIA-only plants: 92.4% match generation within 1%; the US total differs by 12 GWh | — |


#### Combining CEMS and EIA data

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| Subplant crosswalk | Not used | Links CEMS units, EIA boilers, and generators into subplants | OGE is more accurate: prevents double counting and gaps when combining CEMS and EIA data | Enables the next section; no separate impact | — |
| Plants where only some combustion units report to CEMS | Generators at CEMS-reporting plants are included only if they're wind, solar, hydro, nuclear, geothermal, or on a list of biomass units. Fuel and emissions from other non-CEMS units at these plants are not counted, although their generation is. | Fuel and emissions for non-CEMS subplants come from EIA-923 | OGE is more accurate: eGRID understates heat input, emissions, and emission rates at these plants | 176 plants. OGE adds 262 M MMBtu, 16.9 Mt CO2, and 22.3 TWh from non-CEMS units. Plant heat input +4.9% and CO2 +11% in OGE (higher at 81 plants, lower at 28). | APS West Phoenix (117): heat input 16.3 M vs 26.4 M MMBtu (OGE +62%); CO2 0.97 vs 1.56 Mt (OGE +61%) |
| Subplants where only some units report to CEMS | Only the reporting units are counted, using their CEMS values. As in the row above, non-reporting units at a CEMS plant are left out of the unit file, so their fuel and emissions are missing, but the generation they share with the reporting units is counted in full from EIA-923. | EIA-923 totals for the whole subplant, shaped by the partial CEMS data | OGE is more complete: eGRID understates heat input and emissions relative to generation, so emission rates are too low. eGRID's values for the measured units themselves are more precise. | Included in the row above | Polk (7242): heat input 39.8 M vs 41.7 M MMBtu (OGE +4.7%); CO2 2.37 vs 2.44 Mt (OGE +3.3%) |

#### Emissions estimates for units that don't report to CEMS

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| Fuel resolution | A unit's total heat input × the emission factor of its single primary fuel | Emissions calculated separately for each fuel burned by each generator | OGE is more accurate for units that burn more than one fuel | 197 multi-fuel EIA-only plants: totals nearly equal (8.01 vs 7.96 Mt), but 63% differ by more than 1% because the errors offset. At gas-primary plants, eGRID's CO2 per MMBtu equals the gas factor (0.0584 t) vs 0.0614 t in OGE, which includes oil burned. | West Riverside Energy Center (64020): CO2 1.65 vs 1.70 Mt (OGE +3.0%) |
| CH4 and N2O | Emissions for all plants calculated from EIA-923 fuel use × emission factors | For CEMS-reporting units, CH4 and N2O calculated from CEMS heat input × emission factors; All other units use EIA-923 fuel for calculation. | OGE is slightly more accurate for CEMS plants (consistent with measured heat input) | CEMS plants: OGE +1.1%. EIA-only single-fuel plants: within 0.4%. | — |
| SO2 emission estimates | Calculated for each unit's primary fuel only. Fuel sulfur content is the boiler's annual average from EIA-923; if missing, a state average for coal or a national average for other fuels. Removal efficiency is the highest value reported for the boiler in EIA-923. | Calculated for each fuel burned, by month, using the same emission factor sources. Fuel sulfur content is the boiler's month-specific reported value from EIA-923; if missing, a state-month average, then a national annual average, then the previous year's national average. Removal efficiency is weighted by each control device's operating hours, and fluidized bed boilers use factors that already include controls. | OGE is more precise: monthly, fuel-specific sulfur content and emissions from all fuels, not just the primary fuel | EIA-only single-fuel plants: SO2 +4.7% in total, mostly from different sulfur-content values | Bethel (6566): SO2 59 vs 107 t (OGE +81%) |
| Geothermal type assignment | One geothermal type per plant | Geothermal type assigned to each generator: binary (`BT` prime mover), dry steam only at the Geysers, flash otherwise | OGE is more accurate for plants with more than one type | Not isolated from eGRID's geothermal table error (2c) | — |
| Fuel cells | CO2 set to zero. Fuel cells are also left out of CH4 and N2O emissions and out of combustion heat input (`plant_file_create.R`). | Fuel-specific emission factors for CO2, CH4, and N2O (most fuel cells run on natural gas) | OGE is more accurate | 151 plants: CO2 44 kt vs 1,102 kt. 130 plants differ by more than 1% for this reason alone. | Red Lion Energy Center (58433): no CO2 vs 71 kt |
| Municipal solid waste fuel codes | Biogenic and non-biogenic portions (MSB, MSN) are combined into MSW | MSB and MSN kept separate, each with its own biomass treatment | OGE is more precise for both total CO2 and for biogenic CO2 adjustments. | 59 plants have primary fuel MSW in eGRID and MSN in OGE | — |
| "Other" (OTH) fuel | Natural gas emission factor applied to OTH | OTH is reassigned to specific fuels based on heat content and plant type | OGE is more accurate for the reassigned plants | 13 reassigned plants: CO2 6.2 vs 8.9 Mt (OGE +43%) | Swift Creek Chemical Complex (50474): CO2 201 vs 353 kt (OGE +76%) |

#### Adjustments

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| Biomass CO2 adjustment at plants that don't report to CEMS | Plant level: CO2 from EIA-923 biomass fuel use is subtracted from total plant CO2, with a floor of zero | CO2 from each biomass fuel is set to zero for each generator prior to aggregation | OGE is more accurate. eGRID can over-subtract where its total CO2 was estimated from a fossil primary fuel (see fuel resolution). See 2d for CEMS plants. | 621 plants with biomass: CO2 removed 132.4 vs 136.4 Mt. 288 plants differ by more than 1%. | WestRock Covington (50900): unadjusted CO2 1.23 vs 2.66 Mt; adjusted 0 vs 0.78 Mt |
| CHP electric allocation | Applied to plants flagged as CHP, using one annual factor based on total plant generation | Applied to every subplant using monthly factors (hourly for CEMS data) based on each subplant's own generation | OGE is more accurate: the factor reflects the combustion units that produce the heat, and it doesn't depend on a CHP flag | 919 eGRID CHP plants: CO2 removed by the CHP adjustment 103.8 vs 119.8 Mt. Median allocation factor 0.44 vs 0.43; 43 plants differ by more than 0.05. 150 plants differ by more than 1% for this reason alone. | East River (2493): allocation factor 0.63 vs 0.36 |
| Global warming potentials for CO2e | IPCC AR5 without climate-carbon feedback. Earlier data years used other assessments; see table below. | IPCC AR6. Earlier data years used other assessments; see table below. | OGE is more current: it uses the most recent IPCC assessment available in each data year, while eGRID has lagged behind (eGRID follows the EPA GHG Inventory convention). The two datasets don't use the same GWPs in any data year. | US CO2e: +0.008%. Biomass plants: median +0.1%, up to +2%. | Kentucky Mills (55429): OGE CO2e is 2,133 t with AR5 vs 2,174 t with AR6 (+1.9%) |

**IPCC assessment report used for CO2e GWPs, by data year.** Year in parentheses indicates the year that the assessment report was published.

| Data years | eGRID | OGE |
|---|---|---|
| 2005–2006 | SAR (1995) | TAR (2001) |
| 2007–2013 | SAR (1995) | AR4 (2007) |
| 2014–2017 | SAR (1995) | AR5 (2014) with climate-carbon feedback |
| 2018–2020 | AR4 (2007) | AR5 (2014) with climate-carbon feedback |
| 2021–2022 | AR4 (2007) | AR6 (2021) |
| 2023 | AR5 (2014) without climate-carbon feedback | AR6 (2021) |

#### Classification and aggregation

| Topic | eGRID2023 | OGE | Effect | Observed impact (2023) | Example |
|---|---|---|---|---|---|
| Plant primary fuel | Each unit is assigned one primary fuel, and its entire heat input counts toward that fuel regardless of what it actually burned. For CEMS units this is CEMS-measured heat input and the generator's designated EIA-860 fuel. | EIA-923 fuel consumed for electricity, summed for each fuel actually burned | OGE is more accurate for multi-fuel units and plants. The two are equivalent for single-fuel units. | 173 plants differ (98.9% of EIA-only and 96.6% of CEMS plants agree): 59 MSW vs MSN labels, 39 natural gas vs distillate oil at dual-fuel plants, 9 PRG vs OG, 8 MWH vs SUN | West Lorain (2869): DFO vs NG |
| Balancing area assignment | Plants without a BA are assigned "NA" | Missing BAs are inferred from the BA name, utility, or transmission owner; plants outside a BA get state-based codes (e.g. AKMS, HIMS) | Mostly labeling; OGE is slightly more complete | 178 plants differ (7.1 TWh), nearly all Alaska and Hawaii plants labeled "NA" in eGRID; 18 are assigned to a named BA (HECO, CEA) in OGE | Hamakua Energy Plant (55369): "NA - HI" vs HECO |

## 3. What eGRID provides that OGE does not (yet)

A column-by-column mapping of eGRID2023's 1,114 data fields to OGE outputs finds that 568 fields have a direct or derivable OGE equivalent and 546 don't. The main gaps, grouped by theme:

- **Regional aggregations.**
  - **State (ST, 169 fields):** OGE doesn't publish state totals, but 127 of the 169 fields can be derived by summing OGE plant data by state. The other 42 are nonbaseload metrics (31), mercury (8), the "other fossil" fuel category (2), and the state FIPS code.
  - **eGRID subregion and NERC region (SRL/NRL, 338 fields):** these need a subregion and NERC assignment for each plant. NERC region is reported in EIA-860 (it matches eGRID for 98% of plants). eGRID subregions are assigned by EPA and would need to be adopted or replicated.
- **Puerto Rico:** 62 plants, 17.5 TWh (6.8 TWh oil), and 13.5 Mt CO2 in eGRID2023, e.g. AES Puerto Rico (61082, 2.65 TWh, 3.0 Mt). OGE removes Puerto Rico from its EIA and CEMS inputs; adding it means removing that filter and handling the lack of EIA-930 data for hourly shaping and consumption-based rates ([#79](https://github.com/singularity-energy/open-grid-emissions/issues/79)).
- **Nonbaseload metrics:** 31 fields in each aggregated file (nonbaseload generation, emission rates, and generation and resource mix by fuel), plus plant nonbaseload generation. These are computable from existing plant data using eGRID's capacity-factor method, but OGE doesn't calculate them.
- **Other derived metrics OGE doesn't publish:**
  - ozone-season (May–Sep) values
  - input emission rates (lb/MMBtu)
  - combustion-only heat input, generation, and output rates
  - plant heat rate
  - resource mix percentages, including renewable/nonrenewable and combustion/non-combustion groupings
  
  All of these can be derived from existing OGE monthly or annual outputs or subplant data.
- **Unit- and generator-level files:** OGE's finest published level is the subplant. Unit-level CEMS data and generator-level EIA-923 data exist within the OGE pipeline but aren't published as annual tables. Unit and generator attributes (status, firing type, controls, online and retirement years, capacity) are available in PUDL's EIA-860 tables. A few smaller fields would also need to be added: the EPA/CAPD programs each unit is subject to (e.g. Acid Rain Program, Cross-State Air Pollution Rule), unit stack height (from EIA-860), and generator capacity factor (annual net generation ÷ (nameplate capacity × 8,760)).
- **Plant attributes:** utility, transmission owner, sector, NERC region, ISO/RTO, FIPS codes, CHP and biomass flags, capacity factor. Most are available directly from EIA-860 (via PUDL) or derivable from OGE data.
- **Mercury:** unit-level mercury emissions from EPA CAMD, reported for 357 units (4,757 lb) in eGRID2023. The plant and regional files include mercury fields, but they're blank in eGRID2023.
- **Demographics:** EJScreen demographic indicators within 3 miles of each plant.
- **Grid gross loss:** transmission and distribution loss factors by interconnection, from EIA State Electricity Profiles.
- **Formats and supporting resources:** the Excel workbook, summary tables, a PDF technical guide, the eGRID Explorer, and GIS files for subregions.
