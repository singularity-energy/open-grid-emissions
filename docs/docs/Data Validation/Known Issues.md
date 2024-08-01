---
stoplight-id: known_issues
---
## Known Issues
While to our knowledge, the published OGE data represents the highest-quality data of its kind that is publicly available, there are nevertheless a number of known data quality issues with the published OGE data.

We maintain a [list of these known issues on GitHub](https://github.com/singularity-energy/open-grid-emissions/issues), and seek to prioritize and fix these issues over time. 

When running the OGE data pipeline, we run a growing number of [validation checks](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/validation.py) on the data at all stages of the pipeline to ensure that no unexpected transformations are occuring and the data is as complete and accurate as possible. The full results of these validation checks can be viewed in the `data_pipeline.txt` file included with the `data_quality_metrics` results, but we have also summarized warnings raised by these validation checks for the 2022 data below.

### Some subplants only contain a single component of a combined-cycle generator
When identifying subplants that include a combined cycle generator, the subplant should include both parts (i.e. the steam turbine with prime mover CA and the combustion turbine with prime mover CT). In 2022, we identified 120 subplants (out of over 33,000) that only include a single combined cycle part. This could affect the accuracy of our data crosswalks between EIA data and CEMS data. E.g., for 2022:
```log
2024-07-31 09:39:29,349 [WARNING] oge.oge.validation:378 There are 120 subplants that only contain one part of a combined cycle system.
Subplants that represent combined cycle generation should contain both CA and CT parts.
2024-07-31 09:39:30,511 [WARNING] oge.oge.validation:384 
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
For 7 plants, the total generation or fuel allocated to individual generators does not match the total generation or fuel associated with the plant in the raw input data. Possible causes include inconsistent reporting of prime movers or fuel types between EIA-923 tables, mid-year plant retirements, anomalous fuel totals being reported, or other unknown issues. E.g., for 2022:
```log
2024-07-31 09:40:25,359 [WARNING] oge.oge.validation:132 There are 7 plant/fuel combinations with non-zero (missing) net generation or fuel consumed that are different in allocated EIA-923 and EIA-923 Input
2024-07-31 09:40:25,360 [WARNING] oge.oge.validation:141 Percentage Difference:
2024-07-31 09:40:25,360 [WARNING] oge.oge.validation:142 
                                 net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu ba_code
plant_id_eia energy_source_code                                                                                      
1316         NG                            0.087043             0.000000                             0.000000    SWPP
10613        WDS                           0.000000             0.056065                             0.558523    NBSO
54724        GEO                          -1.000000            -1.000000                            -1.000000     IID
58256        MWH                          -1.000000             0.000000                             0.000000     PNM
59817        SUN                          -1.000000            -1.000000                            -1.000000    PACE
59825        WND                          -1.000000            -1.000000                            -1.000000    CPLE
65498        WND                          -1.000000            -1.000000                            -1.000000    MISO
2024-07-31 09:40:25,361 [WARNING] oge.oge.validation:145 EIA-923 Input Totals:
2024-07-31 09:40:25,364 [WARNING] oge.oge.validation:146 
                                 net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu ba_code
plant_id_eia energy_source_code                                                                                      
1316         NG                          413.588001              14274.0                              14274.0    SWPP
10613        WDS                       31047.793091             913514.0                              91699.0    NBSO
54724        GEO                       98803.999023             337120.0                             337120.0     IID
58256        MWH                         -83.000000                  0.0                                  0.0     PNM
59817        SUN                          98.000001                336.0                                336.0    PACE
59825        WND                        9972.000031              34024.0                              34024.0    CPLE
65498        WND                        1157.999996               3950.0                               3950.0    MISO
2024-07-31 09:40:25,364 [WARNING] oge.oge.validation:152 Allocated Totals:
2024-07-31 09:40:25,365 [WARNING] oge.oge.validation:153 
                                 net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu ba_code
plant_id_eia energy_source_code                                                                                      
1316         NG                          449.588001              14274.0                              14274.0    SWPP
10613        WDS                       31047.793091             964730.0                             142915.0    NBSO
54724        GEO                           0.000000                  0.0                                  0.0     IID
58256        MWH                                NaN                  NaN                                  NaN     PNM
59817        SUN                           0.000000                  0.0                                  0.0    PACE
59825        WND                                NaN                  NaN                                  NaN    CPLE
65498        WND                                NaN                  NaN                                  NaN    MISO
```

### Some fuel consumption are associated with an 'OTH' fuel type.
E.g., for 2005:
```log
2024-07-30 22:57:43,411 [WARNING] oge.oge.data_cleaning:282 
            After cleaning energy source codes, some fuel consumption is still 
            associated with an 'OTH' fuel type. This will lead to incorrect emissions 
            calculations. Check the following plants: {3992: 'MISO', 7678: 'SPS', 10004: 'TEC', 10184: 'ERCO', 10205: 'TEC', 10434: 'TEC', 10436: 'ERCO', 10639: 'MISO', 50274: 'IPCO', 50371: 'TEC', 50404: 'ERCO', 50490: 'MISO', 50509: 'CPLE', 50510: 'TEC', 50540: 'CISO', 50624: 'CISO', 50633: 'TEC', 50888: 'PJM', 50991: 'OKMS', 52006: 'MISO', 52063: 'CISO', 52064: 'CISO', 52065: 'ERCO', 54416: 'NJMS', 54806: 'SOCO', 54972: 'TXMS', 55066: 'MISO', 55120: 'MISO', 55557: 'LDWP', 56134: 'CISO', 56294: 'PJM', 10211: 'MISO'}. Assign a fuel
```


### Some plants have no subplant id
E.g., for 2005:
```log
2024-07-30 22:57:43,587 [WARNING] oge.oge.validation:487 There are 96 records for 5 plants without a subplant ID
2024-07-30 22:57:43,588 [WARNING] oge.oge.validation:490 
       plant_id_eia generator_id  subplant_id
1222            469          IC1         <NA>
1223            469          IC2         <NA>
3631           1174            1         <NA>
6019           1956            4         <NA>
7003           2218            1         <NA>
19251         50766         OE11         <NA>
19254         50766         OE14         <NA>
19258         50766         OE21         <NA>

```


### Primary fuel assignment is inconsistent with capacity-based assignment
In OGE, we assign a primary fuel to a plant based on the fuel that was consumed in the highest volume (by mmBTU) in a year, rather than based on the greatest nameplate capacity associated with a fuel (which is used in other contexts). Our validation check flagged 119 plants where our primary fuel assignment is inconsistent with a capacity-based assignment. E.g., for 2022:
```log
2024-07-31 09:40:27,095 [WARNING] oge.oge.validation:211 There are 119 plants where the assigned primary fuel doesn't match the capacity-based primary fuel.
It is possible that these plants will categorized as a different fuel in EIA-930
2024-07-31 09:40:28,088 [WARNING] oge.oge.validation:214 
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

### Missing NOx and SO2 emission factors
Our validation check identifies that there are missing NOx and SO2 emission factors for plants with fuel cell prime movers. This is a known and expected issue - to our knowledge emission factors for fuel cells have not yet been published.

For some years, NOx and SO2 emissions factors are missing for other fuel types. E.g., for 2005:
```log
2024-07-30 22:57:49,180 [WARNING] oge.oge.emissions:1051 The heat content for the following fuels is missing and NOx emissions will not be calculated for these fuel: ['RC']
2024-07-30 22:57:49,585 [WARNING] oge.oge.emissions:649 NOx emission factors are missing for the following records
2024-07-30 22:57:49,585 [WARNING] oge.oge.emissions:650 Missing factors for FC prime movers are currently expected
2024-07-30 22:57:50,460 [WARNING] oge.oge.emissions:651 
    report_date  plant_id_eia energy_source_code prime_mover_code generator_id ba_code
0    2005-01-01         50305                 RC               ST         GE10    MISO
1    2005-01-01         50305                 RC               ST         GEN6    MISO
2    2005-01-01         50305                 RC               ST         GEN7    MISO
3    2005-01-01         50305                 RC               ST         GEN8    MISO
4    2005-01-01         50305                 RC               ST         GEN9    MISO
5    2005-01-01         50481                 RC               ST         TG10     PJM
6    2005-01-01         50481                 RC               ST         TG11     PJM
7    2005-01-01         50481                 RC               ST         TG12     PJM
8    2005-01-01         50481                 RC               ST         TG13     PJM
9    2005-01-01         50481                 RC               ST         TG14     PJM
[...]
```


### Emissions control prior to 2012 has not been integrated into the data pipeline
Note that this may overestimate SO2 and NOx emissions calculated from EIA-923 data prior to 2012.

### Fuel sulfur content is missing for years 2005-2007
Then, 2008-2012 data to derive a plant, state and national fuel sulfur content average. E.g., for 2005 about 20 plants are concerned
```log
2024-07-30 22:57:59,125 [WARNING] oge.oge.emissions:2129 Sulfur content data is missing in EIA-923 for the below units.
2024-07-30 22:58:00,107 [WARNING] oge.oge.emissions:2130 
    plant_id_eia generator_id prime_mover_code energy_source_code ba_code
0          50305         GE10               ST                 RC    MISO
1          50305         GEN6               ST                 RC    MISO
2          50305         GEN7               ST                 RC    MISO
3          50305         GEN8               ST                 RC    MISO
4          50305         GEN9               ST                 RC    MISO
5          50481         TG10               ST                 RC     PJM
6          50481         TG11               ST                 RC     PJM
7          50481         TG12               ST                 RC     PJM
8          50481         TG13               ST                 RC     PJM
9          50481         TG14               ST                 RC     PJM
[...]
```


### Global extreme are detected in CEMS timeseries
We identify where gross generation, fuel consumption, and CO2 emission timeseries in the reported CEMS data may be anomalous. The check reports the number of observations where the hourly value exceeds the median value for the unit by a factor of 10 or more, as well as the mean deviation of these extreme values from the median. E.g., for 2022:
```log
2024-07-31 09:43:33,124 [WARNING] oge.oge.validation:2033 Global extreme detected in CEMS time series
2024-07-31 09:43:33,129 [WARNING] oge.oge.validation:2034 
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

However, there are a number of plants that have positive net generation but zero gross generation. E.g., for 2022:
```log
2024-07-31 09:55:19,489 [WARNING] oge.oge.validation:524 There are 946 subplants at 405 plants for which there is zero gross generation associated with positive net generation.
2024-07-31 09:55:21,242 [WARNING] oge.oge.validation:527 
   plant_id_eia  subplant_id report_date  gross_generation_mwh  net_generation_mwh data_source ba_code
0             7            1  2022-11-01                   0.0              1180.0        both    SOCO
1            10            3  2022-03-01                   0.0                 3.8        both    SOCO
2            10            3  2022-04-01                   0.0               294.4        both    SOCO
3            10            3  2022-08-01                   0.0              1431.3        both    SOCO
4            10            3  2022-09-01                   0.0               523.7        both    SOCO
[...]
```
There are also a number of plants that have zero net generation associated with positive gross generation. E.g., for 2022:
```log
2024-07-31 09:55:21,243 [WARNING] oge.oge.validation:573 There are 76 subplants at 62 plants for which there is zero net generation associated with positive gross generation.
2024-07-31 09:55:23,092 [WARNING] oge.oge.validation:576 
   plant_id_eia  subplant_id report_date  gross_generation_mwh  net_generation_mwh data_source ba_code
0           141            3  2022-04-01                6733.0                 0.0        both     SRP
1           141            3  2022-05-01                4482.0                 0.0        both     SRP
2           141            3  2022-06-01               11353.0                 0.0        both     SRP
3           141            3  2022-07-01               16718.0                 0.0        both     SRP
4           141            3  2022-08-01               10269.0                 0.0        both     SRP
[...]
```
There are also a number of plants where the total annual net generation is substantially higher than the reported gross generation. E.g., for 2022:
```log
2024-07-31 09:55:23,145 [WARNING] oge.oge.validation:628 The following plants have annual net generation that is >125% or <50% of annual gross generation:
2024-07-31 09:55:24,013 [WARNING] oge.oge.validation:631 
     plant_id_eia  gross_generation_mwh  net_generation_mwh  annual_plant_ratio ba_code
71          54096               18039.0            162310.0            8.997727    SOCO
59          10567               17734.0            134118.3            7.562778    ISNE
31           2831                1913.0              7158.4            3.741976     PJM
105         55470             1200399.0           3786786.0            3.154606    ERCO
120         65372                7427.0             17138.7            2.307621    ERCO
[...]
```
These warnings could be the result of incorrect crosswalking between subplants, inconsistent or inaccurate reporting of generation data by the plant, or an incorrect monthly allocation of annually-reported EIA-923 data by the EIA. 

### Several plants are missing default gross-to-net conversion factors
E.g., for 2022:
```log
2024-07-31 09:55:48,222 [WARNING] oge.oge.gross_to_net_generation:114 The following subplants are missing default GTN ratios. Using a default value of 0.97
2024-07-31 09:55:49,232 [WARNING] oge.oge.gross_to_net_generation:117 
    plant_id_eia  subplant_id ba_code
0           2535            1    NYIS
1           2535            2    NYIS
2           6082            1    NYIS
3          50240            4    MISO
4          50240            5    MISO
5          50240            6    MISO
6          50240            7    MISO
7          50282            1     PJM
8          50852            1     PJM
9          50900            3     PJM
10         50900            7     PJM
11         55641            2    MISO
12          6019            2     PJM
13          6019            3     PJM
14          6264            2     PJM
[...]
```

### Our calculated net generation does not match reported net generation for several plants
Currently, we trust CEMS gross generation data more than EIA-923 net generation data, so in certain instances where these are inconsistent with each other, we use default gross-to-net conversion factors to calculate net generation from the CEMS data, which means that the net generation data will not match what is reported in EIA-923. E.g., for 2022:
```log
2024-07-31 09:55:53,446 [WARNING] oge.oge.validation:712 There are 5 plants where calculated annual net generation does not match EIA annual net generation.
2024-07-31 09:55:54,483 [WARNING] oge.oge.validation:715 
   plant_id_eia  net_generation_mwh_eia  net_generation_mwh_calc  pct_error                                                             gtn_method ba_code
0          1391               2979899.1              3436941.351   0.153375                            [6_default_eia_ratio, 5_annual_fleet_ratio]    MISO
1          8906                828822.0              1713570.706   1.067477                            [5_annual_fleet_ratio, 6_default_eia_ratio]    NYIS
2         55088               4447878.7              2704114.030  -0.392044  [5_annual_fleet_ratio, 2_annual_plant_ratio, 1_annual_subplant_ratio]    MISO
3         55799                 81116.4                61396.335  -0.243108                                                 [5_annual_fleet_ratio]    MISO
4         57865                 41887.8               753792.480  16.995514                                                  [6_default_eia_ratio]    SWPP
[...]
```

### Some data outputs are missing complete monthly or hourly values
Currently, there are a number of plants for which there may not be 12 reported monthly values or 8760 reported hourly values in the OGE outputs. These likely represent missing input data, but we are working to fill these missing timestamps and ensure that no data is being mistakenly dropped from the data pipeline.

In some of the regional data, the validation checks flagged certain fuel types for which there are more than 8760 reported values in a year. This appears to result when a BA spans multiple time zones. We are investigating this in more detail to correct this issue.

### Some plants have extreme emission factors
There are approximately 30 plants that have extremely high or low calculated emission rates. The extreme low values appear to mostly occur when a clean plant (eg nuclear) has a fossil-based backup generator that ocassionally runs, resulting in a small amount of emissions relative to the total generation. The extreme high values appear to mostly occur at natural gas and petroleum plants with a small amount of net generation. These may represent peaker plants or spinning reserves, which consume a large amount of fuel throughout the year, but generate only a small amount of net electricity. 

### Negative consumed emission rates in HST
For the small balancing area of HST (City of Homestead, FL), there are 62 hours when the calculated consumed emission rate is negative. We are currently investigating this issue, but did not correct this prior to the data release given that it affected only a small number of hours in a small balancing area. 