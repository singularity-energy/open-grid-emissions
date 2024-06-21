---
stoplight-id: known_issues
---
## Known Issues
While to our knowledge, the published OGE data represents the highest-quality data of its kind that is publicly available, there are nevertheless a number of known data quality issues with the published OGE data.

We maintain a [list of these known issues on GitHub](https://github.com/singularity-energy/open-grid-emissions/issues), and seek to prioritize and fix these issues over time. 

When running the OGE data pipeline, we run a growing number of [validation checks](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/validation.py) on the data at all stages of the pipeline to ensure that no unexpected transformations are occuring and the data is as complete and accurate as possible. The full results of these validation checks can be viewed in the `data_pipeline.txt` file included with the `data_quality_metrics` results, but we have also summarized warnings raised by these validation checks for the 2022 data below.

### Some subplants only contain a single component of a combined-cycle generator
When identifying subplants that include a combined cycle generator, the subplant should include both parts (i.e. the steam turbine with prime mover CA and the combustion turbine with prime mover CT). In 2022, we identified 120 subplants (out of over 33,000) that only include a single combined cycle part. This could affect the accuracy of our data crosswalks between EIA data and CEMS data. 

```log
                         unique_cc_pms ba_code
plant_id_eia subplant_id                      
56           4                      CT    SOCO
             5                      CA    SOCO
658          1                      CA    FMPP
             7                      CT    FMPP
748          3                      CA     TEC
1010         2                      CA    MISO
             10                     CT    MISO
3526         1                      CA    TXMS
             3                      CT    TXMS
3527         1                      CT    TXMS
             2                      CA    TXMS
10188        1                      CA     TEC
             2                      CT     TEC
             3                      CT     TEC
[...]
```

### Some generation and fuel data from EIA-923 is being erroneously allocated
For 7 plants, the total generation or fuel allocated to individual generators does not match the total generation or fuel associated with the plant in the raw input data. Possible causes include inconsistent reporting of prime movers or fuel types between EIA-923 tables, mid-year plant retirements, anomalous fuel totals being reported, or other unknown issues. 

```log
2024-04-03 12:47:46,436 [WARNING] oge.oge.validation:128 There are 7 plant/fuel combinations with non-zero (missing) net generation or fuel consumed that are different in allocated EIA-923 and EIA-923 Input
2024-04-03 12:47:46,436 [WARNING] oge.oge.validation:137 Percentage Difference:
2024-04-03 12:47:46,437 [WARNING] oge.oge.validation:138
                                 net_generation_mwh  fuel_consumed_mmbtu fuel_consumed_for_electricity_mmbtu ba_code
plant_id_eia energy_source_code                                                                                      
1316         NG                        8.704314e-02             0.000000                             0.000000    SWPP
10613        WDS                       1.171735e-16             0.056065                             0.558523    NBSO
54724        GEO                      -1.000000e+00            -1.000000                            -1.000000     IID
58256        MWH                      -1.000000e+00             0.000000                             0.000000     PNM
59817        SUN                      -1.000000e+00            -1.000000                            -1.000000    PACE
59825        WND                      -1.000000e+00            -1.000000                            -1.000000    CPLE
65498        WND                      -1.000000e+00            -1.000000                            -1.000000    MISO
2024-04-03 12:47:46,437 [WARNING] oge.oge.validation:141 EIA-923 Input Totals:
2024-04-03 12:47:46,439 [WARNING] oge.oge.validation:142 
                                 net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu ba_code
plant_id_eia energy_source_code                                                                                      
1316         NG                             413.588              14274.0                              14274.0    SWPP
10613        WDS                          31047.793             913514.0                              91699.0    NBSO
54724        GEO                          98804.000             337120.0                             337120.0     IID
58256        MWH                            -83.000                  0.0                                  0.0     PNM
59817        SUN                             98.000                336.0                                336.0    PACE
59825        WND                           9972.000              34024.0                              34024.0    CPLE
65498        WND                           1158.000               3950.0                               3950.0    MISO
2024-04-03 12:47:46,439 [WARNING] oge.oge.validation:148 Allocated Totals:
2024-04-03 12:47:46,441 [WARNING] oge.oge.validation:149 
                                 net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu ba_code
plant_id_eia energy_source_code                                                                                      
1316         NG                             449.588              14274.0                              14274.0    SWPP
10613        WDS                          31047.793             964730.0                             142915.0    NBSO
54724        GEO                              0.000                  0.0                                  0.0     IID
58256        MWH                                NaN                  NaN                                  NaN     PNM
59817        SUN                              0.000                  0.0                                  0.0    PACE
59825        WND                                NaN                  NaN                                  NaN    CPLE
65498        WND                                NaN                  NaN                                  NaN    MISO
```

### Primary fuel assignment is inconsistent with capacity-based assignment
In OGE, we assign a primary fuel to a plant based on the fuel that was consumed in the highest volume (by mmBTU) in a year, rather than based on the greatest nameplate capacity associated with a fuel (which is used in other contexts). Our validation check flagged 119 plants where our primary fuel assignment is inconsistent with a capacity-based assignment.
```log
2024-04-03 12:47:48,073 [WARNING] oge.oge.validation:206 There are 119 plants where the assigned primary fuel doesn't match the capacity-based primary fuel.
It is possible that these plants will categorized as a different fuel in EIA-930
2024-04-03 12:47:48,742 [WARNING] oge.oge.validation:209 
     plant_id_eia plant_primary_fuel_from_capacity_mw plant_primary_fuel plant_primary_fuel_from_capacity_mw_category plant_primary_fuel_category ba_code
0              63                                 DFO                WAT                                    petroleum                       hydro    AKMS
1              91                                 DFO                WAT                                    petroleum                       hydro    AKMS
2             160                                  NG                SUB                                  natural_gas                        coal    WALC
3             165                                  NG                SUB                                  natural_gas                        coal    SWPP
4             557                                  JF                WAT                                    petroleum                       hydro    ISNE
5             667                                 BIT                 NG                                         coal                 natural_gas     JEA
6             899                                  NG                WAT                                  natural_gas                       hydro    MISO
7             944                                  NG                WND                                  natural_gas                        wind    MISO
8             955                                 DFO                WAT                                    petroleum                       hydro    MISO
9             971                                 DFO                 NG                                    petroleum                 natural_gas    MISO
[...]
```

### Missing NOx and SO2 emission factors for Fuel Cells
Our validation check identifies that there are missing NOx and SO2 emission factors for plants with fuel cell prime movers. This is a known and expected issue - to our knowledge emission factors for fuel cells have not yet been published.


### Global extreme are detected in CEMS timeseries
We identify where gross generation, fuel consumption, and CO2 emission timeseries in the reported CEMS data may be anomalous. The check reports the number of observations where the hourly value exceeds the median value for the unit by a factor of 10 or more, as well as the mean deviation of these extreme values from the median.
```log
2024-04-03 12:51:08,337 [WARNING] oge.oge.validation:2026 Global extreme detected in CEMS time series
2024-04-03 12:51:08,343 [WARNING] oge.oge.validation:2027 
                                   gross_generation_mwh                fuel_consumed_mmbtu                   co2_mass_lb                ba_code
                                         GLOBAL_EXTREME MEAN_DEVIATION      GLOBAL_EXTREME MEAN_DEVIATION GLOBAL_EXTREME MEAN_DEVIATION        
plant_id_eia emissions_unit_id_epa                                                                                                             
315          3                                    350.0           12.0                 NaN            NaN            NaN            NaN    CISO
             4                                    297.0           12.1                 NaN            NaN            NaN            NaN    CISO
356          5                                     55.0           11.5                 NaN            NaN            NaN            NaN    CISO
             6                                    122.0           11.5                 NaN            NaN            NaN            NaN    CISO
377          5                                      NaN            NaN               113.0           10.6          157.0           10.7    LDWP
563          14B                                    NaN            NaN                 4.0           15.8            4.0           16.4    ISNE
673          S-3                                    NaN            NaN                 6.0           62.8            6.0           66.1    FMPP
874          5                                      NaN            NaN                44.0           37.9           44.0           37.7     PJM
1554         1                                      NaN            NaN                14.0           12.1           28.0           13.6     PJM
             4                                      NaN            NaN                12.0           14.2           12.0           19.2     PJM
[...]
```

### Certain subplants have inconsistent gross and net generation values
In general, a subplant's net generation should be some smaller percentage of the total gross generation. 

However, there are a number of plants that have positive net generation but zero gross generation. 
```log
2024-04-03 12:59:18,212 [WARNING] oge.oge.validation:519 There are 947 subplants at 405 plants for which there is zero gross generation associated with positive net generation.
2024-04-03 12:59:19,683 [WARNING] oge.oge.validation:522 
   plant_id_eia  subplant_id report_date  gross_generation_mwh  net_generation_mwh data_source ba_code
0             7            1  2022-11-01                   0.0              1180.0        both    SOCO
1            10            3  2022-03-01                   0.0                 3.8        both    SOCO
2            10            3  2022-04-01                   0.0               294.4        both    SOCO
3            10            3  2022-08-01                   0.0              1431.3        both    SOCO
4            10            3  2022-09-01                   0.0               523.7        both    SOCO
[...]
```
There are also a number of plants that have zero net generation associated with positive gross generation
```log
2024-04-03 12:59:19,684 [WARNING] oge.oge.validation:568 There are 76 subplants at 62 plants for which there is zero net generation associated with positive gross generation.
2024-04-03 12:59:20,952 [WARNING] oge.oge.validation:571 
   plant_id_eia  subplant_id report_date  gross_generation_mwh  net_generation_mwh data_source ba_code
0           141            3  2022-04-01                6733.0                 0.0        both     SRP
1           141            3  2022-05-01                4482.0                 0.0        both     SRP
2           141            3  2022-06-01               11353.0                 0.0        both     SRP
3           141            3  2022-07-01               16718.0                 0.0        both     SRP
4           141            3  2022-08-01               10269.0                 0.0        both     SRP
[...]
```
There are also a number of plants where the total annual net generation is substantially higher than the reported gross generation:
```log
2024-04-03 12:59:20,972 [WARNING] oge.oge.validation:620 The following plants have annual net generation that is >125% of annual gross generation:
2024-04-03 12:59:21,617 [WARNING] oge.oge.validation:623 
    plant_id_eia  gross_generation_mwh  net_generation_mwh  annual_plant_ratio ba_code
40         54096               18039.0            162310.0            8.997727    SOCO
28         10567               17734.0            134118.2            7.562772    ISNE
10          2831                1913.0              7158.4            3.741976     PJM
68         55470             1200399.0           3786786.0            3.154606    ERCO
80         65372                7427.0             17138.7            2.307621    ERCO
[...]
```
These warnings could be the result of incorrect crosswalking between subplants, inconsistent or inaccurate reporting of generation data by the plant, or an incorrect monthly allocation of annually-reported EIA-923 data by the EIA. 

### Several plants are missing default gross-to-net conversion factors
```log
2024-04-03 12:59:43,968 [WARNING] oge.oge.gross_to_net_generation:92 The following subplants are missing default GTN ratios. Using a default value of 0.97
2024-04-03 12:59:44,662 [WARNING] oge.oge.gross_to_net_generation:95 
    plant_id_eia  subplant_id ba_code
0          50852            1     PJM
1          50900            3     PJM
2          50900            7     PJM
3          55641            2    MISO
4          50733            2    MISO
5          50733            3    MISO
6          50733            4    MISO
7          50733            5    MISO
8          50733            6    MISO
9          50733            7    MISO
10         50733            8    MISO
11         50733            9    MISO
12         50900            4     PJM
13         50900            8     PJM
```

### Our calculated net generation does not match reported net generation for several plants
Currently, we trust CEMS gross generation data more than EIA-923 net generation data, so in certain instances where these are inconsistent with each other, we use default gross-to-net conversion factors to calculate net generation from the CEMS data, which means that the net generation data will not match what is reported in EIA-923.

```log
2024-04-03 12:59:48,418 [WARNING] oge.oge.validation:704 There are 5 plants where calculated annual net generation does not match EIA annual net generation.
2024-04-03 12:59:49,106 [WARNING] oge.oge.validation:707 
   plant_id_eia  net_generation_mwh_eia  net_generation_mwh_calc  pct_error                                                               gtn_method ba_code
0          1391               2979899.1              3494770.869   0.172782                               [6_default_eia_ratio, 5_annual_fleet_ratio]    MISO
1          3406                845833.6              1229468.715   0.453559  [5_annual_fleet_ratio, 6_default_eia_ratio, 4_annual_plant_shift_factor]     TVA
2          8906                828822.0              1737494.717   1.096342                               [5_annual_fleet_ratio, 6_default_eia_ratio]    NYIS
3         55075               1276484.0              1819931.568   0.425738                                                    [5_annual_fleet_ratio]    MISO
4         57865                 41887.8               753792.480  16.995514                               [5_annual_fleet_ratio, 6_default_eia_ratio]    SWPP
```

### Some data outputs are missing complete monthly or hourly values
Currently, there are a number of plants for which there may not be 12 reported monthly values or 8760 reported hourly values in the OGE outputs. These likely represent missing input data, but we are working to fill these missing timestamps and ensure that no data is being mistakenly dropped from the data pipeline.

In some of the regional data, the validation checks flagged certain fuel types for which there are more than 8760 reported values in a year. This appears to result when a BA spans multiple time zones. We are investigating this in more detail to correct this issue.

### Some plants have extreme emission factors
There are approximately 30 plants that have extremely high or low calculated emission rates. The extreme low values appear to mostly occur when a clean plant (eg nuclear) has a fossil-based backup generator that ocassionally runs, resulting in a small amount of emissions relative to the total generation. The extreme high values appear to mostly occur at natural gas and petroleum plants with a small amount of net generation. These may represent peaker plants or spinning reserves, which consume a large amount of fuel throughout the year, but generate only a small amount of net electricity. 

### Negative consumed emission rates in HST
For the small balancing area of HST (City of Homestead, FL), there are 62 hours when the calculated consumed emission rate is negative. We are currently investigating this issue, but did not correct this prior to the data release given that it affected only a small number of hours in a small balancing area. 