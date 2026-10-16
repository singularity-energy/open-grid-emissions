---
stoplight-id: physics_reconciliation_eia930
---
<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

## Physics-based Reconciliation of EIA-930 Data

EIA-930 data are not always physically consistent. In a given hour, a balancing authority's (BA's) demand plus its net exports may not equal its net generation, its total interchange may not equal the sum of its interchange with each neighbor, two neighboring BAs may disagree about the electricity that flowed between them, and a BA's generation by energy source may not add up to its total net generation. EIA-930 data also contain missing hours and anomalous values.

OGE fixes these transient issues with a physics-based reconciliation, after correcting the systematic timestamp and sign issues described in [EIA-930 Data](EIA-930%20Data.md). The method was developed by Jacques de Chalendar, first to track consumption-based emissions in the US electricity system ([de Chalendar, Taggart and Benson, 2019](https://doi.org/10.1073/pnas.1912950116)) and then as a general data reconciliation framework for real-time electricity and emissions tracking ([de Chalendar and Benson, 2021](https://doi.org/10.1016/j.apenergy.2021.117761)). OGE originally used a fork of de Chalendar's [gridemissions](https://github.com/jdechalendar/gridemissions) package, and now uses a port of its cleaning code in `oge.physics_reconciliation`. The cleaning steps give identical results, and the reconciled values are within about 1 MWh of gridemissions' (see [Differences from the gridemissions package](#differences-from-the-gridemissions-package)).

The reconciliation runs on EIA-930 data from October 1 of the previous year through January 2 of the following year, so that the rolling window has data before and after the year of interest.

PUDL's operations table (demand, net generation, and total interchange) can be missing BAs that report other EIA-930 data. For example, SIKE, a generation-only BA that started reporting in June 2025, has generation and interchange data but no operations data. For such BAs, we derive net generation as the sum of their generation by energy source, total interchange as the sum of their interchange with each neighbor, and demand as net generation minus total interchange. Otherwise, the reconciliation would treat them as BAs outside of the US (see below) and discard the data they report.

The reconciliation has three steps.

### 1. Basic cleaning

The first step drops values that are outside of realistic bounds and adds data for BAs outside of the US.

**Realistic bounds.** Values outside of these ranges are treated as missing:

| Series | Range (MWh) |
| --- | --- |
| Demand, net generation, and generation by energy source (default) | -1 to 200,000 |
| Total interchange and interchange between BAs | -100,000 to 100,000 |
| Solar generation | -1,000 to 200,000, then negative values are set to 0 |
| AZPS demand | 1 to 30,000 |
| BANC demand | 1 to 6,500 |
| BANC total interchange | -5,000 to 5,000 |
| CISO net generation | 5,000 to 60,000 |

Solar generation is often reported at around -50 MWh at night, from station use, so negative solar values are set to zero rather than dropped.

**BAs outside of the US.** Some US BAs exchange electricity with BAs in Canada and Mexico (AESO, BCHA, CEN, CFE, HQT, IESO, MHEB, NBSO, and SPC), which do not report to EIA-930. Their total interchange is set to the opposite of the interchange reported by their US neighbors. Their net exports are treated as their net generation, and their net imports as their demand. Their generation is assumed to be hydro for BCHA, HQT, and MHEB, and "other" for the rest. See [Consumption-based Emissions](../Emissions%20Calculations/Consumption-based%20Emissions.md) for how this generation is used.

**Other adjustments.**
* Every interchange series must have a matching series in the opposite direction. If one is missing, it is added as the opposite of the existing series.
* AVRN only exports to BPAT. If AVRN has no interchange data, its exports to BPAT are set to its net generation.
* If a BA reports no generation by energy source, all of its net generation is assigned to "other".

### 2. Rolling window cleaning

The second step drops anomalous values and fills gaps in each time series.

**Filtering.** For each series, we compute a rolling median and the rolling range between the 5th and 95th percentiles, using a centered 20-day (480-hour) window with at least 100 hours of data. Values more than 3 times this range away from the median are dropped. The range is at least 150 MWh, so that small series are not over-filtered. This filter is applied twice, since very large outliers can mask smaller ones.

The original gridemissions package used a mean and standard deviation filter, dropping values more than 4 standard deviations from the mean over a 10-day window. We found that this missed some large spikes, including in BANC in August 2020, because the spikes inflated the standard deviation of the window (see the figure below). We switched to a median and percentile filter based on [Ruggles et al. (2020)](https://doi.org/10.1038/s41597-020-0483-x), and adjusted its parameters until it filtered the BANC spikes without excluding normal data, even in small BAs where fossil fuel generation can be highly variable.

![spikes in hydro time series](https://stoplight.io/api/v1/projects/cHJqOjE1MzAxNA/images/84uL15gr69I)

**Filling gaps.** Electricity data have strong daily patterns, so each missing value is first filled with the average of the same hour on the previous and next days. This is repeated up to four times, one day at a time, which gives smoother transitions across gaps of several days, and then without a limit. Series with no data at all are set to zero. Any remaining gaps are filled by linear interpolation for gaps of up to 3 hours, and with the rolling median otherwise.

### 3. Reconciliation

The last step adjusts the data in each hour so that it is physically consistent, while changing it as little as possible. For each hour, we solve a quadratic program that finds the smallest weighted adjustment to each value. Using the notation below, for each BA _i_:

* _D<sub>i</sub>_ is demand, _G<sub>i</sub>_ is net generation, and _T<sub>i</sub>_ is total interchange (positive for net exports)
* _I<sub>ij</sub>_ is interchange from BA _i_ to BA _j_ (positive for exports from _i_)
* _G<sub>is</sub>_ is net generation from energy source _s_
* &delta; is the adjustment to each value, and _w_ is its weight

The problem is:

$$
\min_{\delta} \sum_k w_k \delta_k^2
$$

subject to, for every BA _i_, neighbor _j_, and energy source _s_:

$$
\begin{aligned}
(D_i + \delta_{D_i}) + (T_i + \delta_{T_i}) &= G_i + \delta_{G_i} \\
T_i + \delta_{T_i} &= \sum_j \left( I_{ij} + \delta_{I_{ij}} \right) \\
I_{ij} + \delta_{I_{ij}} &= -\left( I_{ji} + \delta_{I_{ji}} \right) \\
G_i + \delta_{G_i} &= \sum_s \left( G_{is} + \delta_{G_{is}} \right) \\
D_i + \delta_{D_i} \geq 1, \quad G_i + \delta_{G_i} &\geq 1, \quad G_{is} + \delta_{G_{is}} \geq 1
\end{aligned}
$$

where the sum _k_ is over every value in the hour. The constraints say that demand plus net exports equals net generation, that total interchange equals the sum of interchange with each neighbor, that neighboring BAs agree on the interchange between them, and that generation by energy source adds up to net generation. EIA-930 did not collect generation by energy source before July 1, 2018, so the constraints and adjustments for generation by energy source are only used from that date onwards.

**Weights.** The weight of each value is inversely proportional to the typical size of its series:

$$
w_k = \frac{10{,}000}{\max\left(10, \left|\bar{x}_k\right|\right)}
$$

where $\bar{x}_k$ is the rolling mean of the series over the same 20-day window used for filtering. As a result, the optimization prefers to adjust large series (for example, demand in a large BA) rather than small ones, and the size of an adjustment is roughly proportional to the size of the series.

The problem is solved independently for each hour using [cvxpy](https://www.cvxpy.org/) and the [Clarabel](https://clarabel.org/) solver. The problem only differs between hours in its data and weights, so it is built and compiled once, and then solved for each hour with that hour's values. Clarabel is run with tight tolerances (10<sup>-10</sup>), so the reconciled values are within about 0.01 MWh of the exact solution of the problem.

### Differences from the gridemissions package

* The rolling filter uses a median and percentile range instead of a mean and standard deviation (see above).
* We do not add the constant geothermal (900 MWh) and biomass (600 MWh) generation that gridemissions adds to CISO.
* gridemissions solves the reconciliation with the OSQP solver at its default tolerances, which leaves errors of up to about 1 MWh compared to the exact solution, and occasionally violates the constraints (for example, negative demand). We use Clarabel with tight tolerances instead, so our reconciled values differ from gridemissions' by up to about 1 MWh, and are closer to the exact solution.
* The optimization requires demand, net generation, and generation by energy source to be at least 1 MWh, so a value that should be zero comes out of the reconciliation as about 1 MWh. After reconciliation, we set values with an absolute value below 1.5 MWh to zero, both in the hourly profiles and in the consumption-based emissions calculation.

### Outputs

Each step writes its results to `outputs/{year}/eia930/` as parquet files, in three tables that follow PUDL's EIA-930 tables: `operations` (demand, net generation, and total interchange for each BA), `net_generation_by_energy_source`, and `interchange` (between each pair of BAs). The files are named `eia930_{step}_{table}.parquet`, where the steps are:

| Step | Description |
| --- | --- |
| `unadjusted` | EIA-930 data as loaded from PUDL |
| `raw` | After the manual timestamp and sign adjustments |
| `basic` | After basic cleaning |
| `rolling` | After rolling window cleaning |
| `weights` | The weight of each value in the reconciliation |
| `reconciled_without_energy_sources` | Reconciled data before July 1, 2018 (only when the data include those hours) |
| `reconciled` | Reconciled data, used to build hourly profiles and to calculate consumption-based emissions |

The value of the optimization objective for each hour is written to `eia930_objective.parquet` (and `eia930_objective_without_energy_sources.parquet` before July 1, 2018).

## Future Work, Known Issues, and Open Questions
- Since late 2024, BAs have started reporting energy storage separately, as negative net generation when charging. Negative values below -1 MWh are dropped by the basic cleaning, and the new energy source categories are currently combined into the old ones before reconciliation.
- The reconciliation does not conserve the total energy reported to EIA-930: values that are dropped are re-filled, and the lower bounds of 1 MWh are met by raising other values.
- Reconciliation can add hour-to-hour variation to series that should be flat, such as nuclear generation ([details](https://github.com/singularity-energy/open-grid-emissions/issues/240)).
- Consider reconciling EIA-930 data against OGE's own generation estimates ([details](https://github.com/singularity-energy/open-grid-emissions/issues/220)).
