---
stoplight-id: quality_metrics
---
## Data Quality Metrics
To help data users understand the quality of the published data, the dataset includes several tables that summarize different types of metadata.

### Input Data Source
Indicates what percentage of the output data values came from either EIA-923 data or CEMS data. This percentage is calculated for `subplant-hours `(essentially the number of observations), `net_generation_mwh`, `co2_mass_lb`, and `co2_mass_lb_for_electricity`.

Useful for understanding how much of the hourly data represents measured hourly data (CEMS) versus shaped hourly data (EIA).

### Annually reported EIA Data
Summarize the percent of input EIA data that was reported annually, the percent of output data that uses the annually-reported EIA data, and the percent of output data that uses mixed CEMS and annually-reported EIA data.

### CEMS Gross-to-Net Methods
Indicates what percentage of `gross_generation_mwh` in CEMS was converted to `net_generation_mwh` using each method. Methods are listed in order of preference in the hierarchy.

### Hourly Profile Method
Indicates the percentage of hourly data from each shaping method, including reported cems, partial cems, and each of the fleet-specific shaping methods. This percentage is reported for both `net_generation_mwh` and `co2_mass_lb`.

### DIBA Imputation Performance / National Imputation Performance
These two tables indicate report the results of the cross-validation performed for imputing missing wind and solar profiles. Specifically, the tables include the correlation coefficients between the imputed wind and solar profiles for each BA-month and the actual reported profiles for those regions where profiles are available.

### Plant Metadata
Indicates for each subplant-month:
- the input data source
- the hourly profile method
- the gross-to-net generation method