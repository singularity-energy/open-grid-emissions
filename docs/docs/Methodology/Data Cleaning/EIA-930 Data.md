---
stoplight-id: cleaning_eia930
---

## EIA-930 data
All of the EIA-930 data used by OGE is loaded from tables created by the PUDL project. We use three of the many data types provided by EIA-930:

* Hourly fuel-specific generation for each Balancing Authority (BA)
* Hourly interchange between connected BAs
* Hourly demand for each BA

For the EIA-930 output tables published in PUDL, in addition to the original reported electricity demand, PUDL provide a column with anomalous and missing values replaced by imputed values. Check their [website](https://catalystcoop-pudl.readthedocs.io/en/latest/data_sources/eia930.html#notable-irregularities) for a detailed description of the imputation procedure.

EIA has been collecting fuel-specific generation since July 2018, making this the newest of the data sources and the limiting factor in why our provided dataset begins in 2019. The data is provided to EIA from balancing authorities daily. There have been several reporting improvements since the beginning of collection, and some issues remain with the data reported to EIA from BAs. These issues include systematic issues with reported data and transient inconsistencies between data types and regions. We fix systematic issues using a correlation-based manual inspection process, then fix transient inconsistencies using the physics-based reconciliation process developed by [Chalendar and Benson (2021)](https://doi.org/10.1016/j.apenergy.2021.117761)

### Energy sources
EIA-930 originally reported generation for eight energy sources: coal, natural gas, nuclear, petroleum, hydro, solar, wind, and other, plus unknown. Starting in mid-2024, it split some of these into more specific energy sources. BAs started to report them at different times, and some only report part of them (for example, CISO reports geothermal separately but still reports its batteries in "other").

| Original energy source | New energy sources |
| --- | --- |
| Hydro (including pumped storage) | Hydro (excluding pumped storage), pumped storage |
| Solar | Solar (without integrated battery storage), solar with integrated battery storage |
| Wind | Wind (without integrated battery storage), wind with integrated battery storage |
| Other | Other, battery storage, other energy storage, geothermal |

Each original energy source equals the sum of the new energy sources it was split into, so a BA's generation means the same thing before and after it starts reporting them. OGE keeps the new energy sources as their own series in the [physics-based reconciliation](Physics-based%20Reconciliation.md). Storage is reported as negative net generation when it charges more than it discharges. Some BAs (for example, ERCO) report battery charging as "unknown energy storage", which isn't in the EIA-930 instructions, and OGE keeps this as its own series too.

PUDL renamed hydro, solar, and wind to "hydro excluding pumped storage", "solar without integrated battery storage", and "wind without integrated battery storage" on July 1, 2024 for every BA, so the new names don't mean that a BA has started reporting the new energy sources. OGE combines each pair into one series under the original name.

For hourly profiles, the new energy sources are kept separate in the months after each BA started reporting them (see [Cutover dates](#cutover-dates) below and [Shaping Using Fleet-Specific Profiles](../Assigning%20Hourly%20Profiles%20to%20Monthly%20Data/Shaping%20Using%20Fleet-Specific%20Profiles.md#new-eia-930-energy-sources)), and added back to the energy source that included them before until then.

#### Cutover dates
Until a BA reports a new energy source, its generators of that type are still included in the energy source that included them before. For each BA and new energy source, OGE identifies a cutover date: the first day on which the BA reports a non-zero value, but no earlier than January 1, 2025. Zeros don't count, since some BAs (for example, TVA) report zeros for energy sources they don't have. Two shares of the hours after the cutover date describe how consistently it is reported: the share of hours with a value (including zero), which is low for BAs that leave the energy source blank in some hours (for example, FPC's batteries), and the share of hours with a non-zero value, which is low for energy sources that are rarely used or whose generation is mostly reported elsewhere (for example, CPLE's batteries before October 2025). Battery charging that a BA reports as "unknown energy storage" counts as reporting batteries. Data years before 2025 have no cutover dates.

Detected dates can be replaced, or cleared for a BA that hasn't really switched, in [`eia930_energy_source_cutover_overrides.csv`](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/oge/reference_tables/eia930_energy_source_cutover_overrides.csv).

Each BA and new energy source is also compared with OGE's subplants in the BA with the matching `fuel_category_eia930_expanded` (see [Plant Primary Fuel](../Data%20Aggregation/Plant%20Primary%20Fuel.md#eia-930-fuel-categories)), and flagged as:

Flag | Meaning
---------|---------
`reporting_and_resources` | The BA reports the energy source and has subplants of that type
`resources_not_reporting` | The BA has subplants of that type but doesn't report the energy source, so they are still included in the energy source that included them before
`reporting_without_resources` | The BA reports the energy source but OGE has no subplants of that type in the BA, which points to a difference in BA assignment or in how generators are categorized

The cutover dates and flags are written to `outputs/{year}/eia930/eia930_energy_source_cutover.csv`.

### Systematic data issues
We identified two types of issues with the data: 1) data reported with the wrong timestamp and 2) data reported with the wrong sign (this issue was specific to interchange data). To identify the specific time periods, BAs, and data types with each of these issues, we performed iterative lagged and windowed correlation between data sets that we expected to be highly correlated, fixing issues as we found them.

For generation data, we used CEMS data as a reference to identify BAs where EIA-930 data was reported using the wrong timestamp or in the wrong hour. CEMS data aggregated to the BA level should correlate with fossil fuel generation reported in EIA-930. Even if there are data issues with individual CEMS plants, these should not be widespread enough to cause the systemic offsets between CEMS and 930 data seen in balancing authorities with timestamp reporting issues after aggregating the many CEMS plants in large BAs.

For interchange, we compared interchange reported by neighboring BAs. These interchanges should be negatively correlated without a lag.

To identify lags in correlations, we perform lagged correlations between timeseries we expect to be correlated at lags of -12 hours to +12 hours. When the best correlation is not at lag=0, we inspect the dataset for potential issues.

While 930 data should be reported with end-of-hour timestamps ([https://www.eia.gov/survey/form/eia_930/instructions.pdf](https://www.eia.gov/survey/form/eia_930/instructions.pdf)), three BAs (CISO, PJM, and TEPC) instead report data with start-of-hour timestamps. One BA, SC, seems to report data in local time, instead of UTC, resulting in data that is 4 hours ahead during daylight savings hours and 5 hours ahead during the rest of the year (shown below).

![sc_lag_example.png](https://stoplight.io/api/v1/projects/cHJqOjE1MzAxNA/images/hai1utOhmF0)


The EIA specifies that interchange should be positive when electricity is flowing out of a balancing authority and negative when electricity is flowing in. We found that PJM had inverted the signs of its interchange, an issue that was fixed on October 31, 2019, 4:00 UTC for all trading partners except OVEC, which continues to be inverted until at least the end of 2020. In addition, PJM appears to report erroneous timestamps for its interchange data, which is lagged by 3 hours during daylight savings time and 4 hours during the rest of the year. TEPC’s interchange data is lags by 7 hours.

To ensure that we were not performing corrections based on spurious correlations, we corrected issues only when they were consistent across years (2019 and 2020) and between DST and non-DST times (April-October and December-February). In cases where a data reporting issue was fixed at a specific date, for example, the sign change in PJM interchange data, we manually inspected the data to confirm that the issue was not spurious even though it changed across years. In cases where a timestamp issue was DST-related, we ensured that the lags identified by correlation in DST and non-DST windows was one hour different, consistent with the balancing authority erroneously reporting data with a DST shift.

After performing these corrections, we saw some non-zero lags in small BAs that were not consistent across windows and not explicable by a data issue resolved at some date. In these cases, since we could not identify a specific issue, we did not make corrections to the 930 data.

We have reported these timestamp issues with the EIA, which works with BAs on reporting. We will continue to monitor the data stream and remove our corrections for future data as the incoming data is fixed.

The manual adjustments performed can be found [here](https://github.com/singularity-energy/open-grid-emissions/blob/908bfcea07d00b38bfd7ac2503f947dea09b9595/src/oge/eia930.py#L577).


### Transient issues
The procedure described above works only for major, systematic issues. However, EIA-930 data contains a multitude of smaller issues, including hours with missing data, hours with anomalously large or small data, or hours where data is physically impossible; for example, when the demand in a region is smaller than the sum of the region’s generation and import, or when the import from a region disagrees with the export from its neighbor.

To fix these issues, we use a physics-based reconciliation (in `oge.physics_reconciliation`) ported from the [gridemissions package](https://github.com/jdechalendar/gridemissions) developed by Jacques de Chalendar ([Chalendar and Benson, 2021](https://doi.org/10.1016/j.apenergy.2021.117761)). It drops data outside of reasonable bounds, drops anomalous values using a rolling window filter and fills gaps, and then finds the smallest set of adjustments to the data that make it physically consistent in each hour. See [Physics-based Reconciliation](Physics-based%20Reconciliation.md) for details.

