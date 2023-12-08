---
stoplight-id: cleaning_eia930
---

## EIA-930 data

We use EIA’s form 930 to calcualte hourly profiles for plants which do not report to CEMS and to calculate consumed emission rates. We use two of the many data types provided by EIA-930:

* Hourly fuel-specific generation for each Balancing Authority (BA)
* Hourly interchange between connected BAs

EIA has been collecting fuel-specific generation since July 2018, making this the newest of the data sources and the limiting factor in why our provided dataset begins in 2019. The data is provided to EIA from balancing authorities daily. There have been several reporting improvements since the beginning of collection, and some issues remain with the data reported to EIA from BAs. These issues include systematic issues with reported data and transient inconsistencies between data types and regions. We fix systematic issues using a correlation-based manual inspection process, then fix transient inconsistencies using the physics-based reconciliation process developed by [Chalendar and Benson (2021)](https://doi.org/10.1016/j.apenergy.2021.117761)

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

**Summary of manual data cleaning steps**

* Adjustments to generation timestamps
    * PJM: + 1 hour
    * CISO: + 1 hour
    * TEPC: + 1 hour
    * SC: - 4 hours during daylight savings hours; -5 hours during standard hours
* Adjustments to interchange timestamps
    * PJM: + 4 hours during standard time, + 3 hours during daylight time
    * TEPC: + 7 hours
    * CFE has zero interchange, so the correlation looks strange/misleading
* Flipped interchange sign (x -1)
    * PJM-{CPLE, CPLW, DUK, LGEE, MISO, NYIS, TVA} before Oct 31, 2019, 4:00 UTC (this is all interchange partners except OVEC)
    * PJM-OVEC, all time. Based on OVEC demand - generation, OVEC should be a net exporter to PJM
* Other interchange issues
    * AZPS-SEC interchange before June 1, 2020 does not agree with SEC-AZPS. We assume the SEC-AZPS interchange is correct, and use its inverse for AZPS-SEC prior to 7/1/2020, 3:00 UTC
* Adjustments to all timestamps (after the above adjustments)
        * - 1 hour to move from end-of-hour to start-of-hour

### Transient issues

The procedure described above works only for major, systematic issues. However, EIA-930 data contains a multitude of smaller issues, including hours with missing data, hours with anomalously large or small data, or hours where data is physically impossible; for example, when the demand in a region is smaller than the sum of the region’s generation and import, or when the import from a region disagrees with the export from its neighbor.

To fix these issues, we use the [gridemissions package](https://github.com/jdechalendar/gridemissions) developed by Jaques de Chalendar ([Chalendar and Benson, 2021](https://doi.org/10.1016/j.apenergy.2021.117761)). This package first performs a basic cleaning step where data outside of reasonable bounds is dropped, followed by a rolling cleaning step which drops anomalous values. Finally, the package performs an optimization algorithm to make the data physically consistent in each hour. Specifically, it finds the smallest set of adjustments to the data that ensure that demand in a region is the sum of generation and interchange into that region, and that interchange between each pair of neighboring regions agrees.

The original gridemissions package used a mean and standard deviation rolling filter, dropping values that were outside of 4 standard deviations from the mean over a 10 day window. We found that there were some cases, including in BANC in August, 2020, where large spikes were not caught using this filter because they altered the standard deviation in the rolling window (see the figure below for an example of these spikes). We switched to a median/IQR filter based on the method developed by [Ruggles et al. (2020)](https://doi.org/10.1038/s41597-020-0483-x) and adjusted the filter parameters until we we able to filter the BANC spikes without excluding normal data, even in small BAs where fossil fuel generation can be highly variable. Our filter excluded values further than 3 times the 5th to 95th percentile range from the median over a 20 day window.

![spikes in hydro time series](https://stoplight.io/api/v1/projects/cHJqOjE1MzAxNA/images/84uL15gr69I)


We make two changes to the gridemissions package to exclude default values introduced by the cleaning process that do not reflect the underlying data. We remove the addition of geothermal and biomass fuel categories in CISO, added by the package with default generation estimates. The package also adds a small amount of generation (1.0 MWh) for fuel types in BAs without those fuel types to simplify the physics-based reconciliation. We drop these small non-zero values after cleaning, since they come from BAs which do not actually have any generation of those fuel types.

