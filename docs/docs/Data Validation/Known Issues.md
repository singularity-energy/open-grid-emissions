---
stoplight-id: known_issues
---
## Known Issues
While to our knowledge, the published OGE data represents the highest-quality data of its kind that is publicly available, there are nevertheless a number of known data quality issues with the published OGE data.

We maintain a [list of these known issues on GitHub](https://github.com/singularity-energy/open-grid-emissions/issues), and seek to prioritize and fix these issues over time. 

When running the OGE data pipeline, we run a growing number of [validation checks](https://github.com/singularity-energy/open-grid-emissions/blob/main/src/validation.py) on the data at all stages of the pipeline to ensure that no unexpected transformations are occuring and the data is as complete and accurate as possible. The full results of these validation checks can be viewed in the `data_pipeline.txt` file included with the `data_quality_metrics` results, but we have also summarized warnings raised by these validation checks for the 2022 data below.

### Some subplants only contain a single component of a combined-cycle generator
When identifying subplants that include a combined cycle generator, the subplant should include both parts (i.e. the steam turbine with prime mover CA and the combustion turbine with prime mover CT). In 2022, we identified 259 subplants (out of over 33,000) that only include a single combined cycle part. This could affect the accuracy of our data crosswalks between EIA data and CEMS data. 

```
                         unique_cc_pms
plant_id_eia subplant_id              
3            7.0                    CT
             8.0                    CA
56           3.0                    CT
             4.0                    CA
96           6.0                    CA
341          4.0                    CT
             5.0                    CT
             6.0                    CT
             7.0                    CA
             8.0                    CA
[...]
```

### Some generation and fuel data from EIA-923 is being erroneously allocated
For 8 plants, the total generation or fuel allocated to individual generators does not match the total generation or fuel associated with the plant in the raw input data. Possible causes include inconsistent reporting of prime movers or fuel types between EIA-923 tables, mid-year plant retirements, anomalous fuel totals being reported, or other unknown issues. 

```
2023-12-23 11:27:47,744 [WARNING] oge.oge.validation:79 Allocated EIA-923 doesn't match input data for plants:
2023-12-23 11:27:47,744 [WARNING] oge.oge.validation:80 Percentage Difference:
2023-12-23 11:27:47,744 [WARNING] oge.oge.validation:81 
              net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu
plant_id_eia                                                                              
1316                    0.071287             0.000000                             0.000000
10613                   0.000000             0.005637                             0.033157
50937                   0.000000            -0.121487                            -0.121487
54724                  -1.000000            -1.000000                            -1.000000
54809                   0.000000            -0.302998                            -0.797823
58256                   0.142367             0.000000                             0.000000
59817                  -0.054505            -0.054741                            -0.054741
59825                  -1.000000            -1.000000                            -1.000000
65498                  -1.000000            -1.000000                            -1.000000
2023-12-23 11:27:47,744 [WARNING] oge.oge.validation:82 EIA-923 Input Totals:
2023-12-23 11:27:47,744 [WARNING] oge.oge.validation:83 
              net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu
plant_id_eia                                                                              
1316                      505.00              17176.0                              17176.0
10613                  373515.00            9084933.0                            1544644.0
50937                       0.00               1103.0                               1103.0
54724                   98804.00             337120.0                             337120.0
54809                    9033.11            1493143.0                             567067.0
58256                     583.00               2274.0                               2274.0
59817                    1798.00               6138.0                               6138.0
59825                    9972.00              34024.0                              34024.0
65498                    1158.00               3950.0                               3950.0
2023-12-23 11:27:47,744 [WARNING] oge.oge.validation:86 Allocated Totals:
2023-12-23 11:27:47,744 [WARNING] oge.oge.validation:87 
              net_generation_mwh  fuel_consumed_mmbtu  fuel_consumed_for_electricity_mmbtu
plant_id_eia                                                                              
1316                      541.00              17176.0                              17176.0
10613                  373515.00            9136149.0                            1595860.0
50937                       0.00                969.0                                969.0
54724                       0.00                  0.0                                  0.0
54809                    9033.11            1040724.0                             114648.0
58256                     666.00               2274.0                               2274.0
59817                    1700.00               5802.0                               5802.0
59825                       0.00                  0.0                                  0.0
65498                       0.00                  0.0                                  0.0
```

### Primary fuel assignment is inconsistent with capacity-based assignment
In OGE, we assign a primary fuel to a plant based on the fuel that was consumed in the highest volume (by mmBTU) in a year, rather than based on the greatest nameplate capacity associated with a fuel (which is used in other contexts). Our validation check flagged 173 plants where our primary fuel assignment is inconsistent with a capacity-based assignment.

```
023-12-23 11:27:50,249 [WARNING] oge.oge.validation:141 There are 173 plants where the assigned primary fuel doesn't match the capacity-based primary fuel.
It is possible that these plants will categorized as a different fuel in EIA-930
2023-12-23 11:27:50,249 [WARNING] oge.oge.validation:144 
       plant_id_eia plant_primary_fuel_from_capacity_mw plant_primary_fuel plant_primary_fuel_from_capacity_mw_category plant_primary_fuel_category
23               47                                 BIT                 NG                                         coal                 natural_gas
33               63                                 DFO                WAT                                    petroleum                       hydro
49               91                                 DFO                WAT                                    petroleum                       hydro
81              160                                  NG                SUB                                  natural_gas                        coal
278             460                                  NG                DFO                                  natural_gas                   petroleum
288             477                                 BIT                 NG                                         coal                 natural_gas
305             508                                 SUB                WND                                         coal                        wind
333             552                                 DFO                WAT                                    petroleum                       hydro
336             557                                  JF                WAT                                    petroleum                       hydro
[...]
```

### Missing NOx and SO2 emission factors for Fuel Cells
Our validation check identifies that there are missing NOx and SO2 emission factors for plants with fuel cell prime movers. This is a known and expected issue - to our knowledge emission factors for fuel cells have not yet been published.

### Certain subplants have inconsistent gross and net generation values
In general, a subplant's net generation should be some smaller percentage of the total gross generation. 

However, there are a number of plants that have positive net generation but zero gross generation. 
```
2023-12-23 11:40:28,584 [WARNING] oge.oge.validation:375 There are 949 subplants at 405 plants for which there is zero gross generation associated with positive net generation.
2023-12-23 11:40:28,599 [WARNING] oge.oge.validation:378 
     plant_id_eia  subplant_id report_date  gross_generation_mwh  net_generation_mwh data_source
82              7            0  2022-11-01                   0.0              1180.0        both
144            10            3  2022-01-01                   0.0                 8.1        both
146            10            3  2022-03-01                   0.0                 3.8        both
151            10            3  2022-08-01                   0.0              1431.3        both
[...]
```

There are also a number of plants that have zero net generation associated with positive gross generation
```
2023-12-23 11:40:28,599 [WARNING] oge.oge.validation:412 There are 76 subplants at 62 plants for which there is zero net generation associated with positive gross generation.
2023-12-23 11:40:28,599 [WARNING] oge.oge.validation:415 
      plant_id_eia  subplant_id report_date  gross_generation_mwh  net_generation_mwh data_source
843            141            2  2022-04-01                6733.0                 0.0        both
844            141            2  2022-05-01                4482.0                 0.0        both
845            141            2  2022-06-01               11353.0                 0.0        both
846            141            2  2022-07-01               16718.0                 0.0        both
```
There are also a number of plants where the total annual net generation is substantially higher than the reported gross generation:
```
2023-12-23 11:40:28,646 [WARNING] oge.oge.validation:451 The following plants have annual net generation that is >125% of annual gross generation:
2023-12-23 11:40:28,646 [WARNING] oge.oge.validation:454 
      plant_id_eia  gross_generation_mwh  net_generation_mwh  annual_plant_ratio
772          54096               18039.0            162310.0            8.997727
693          10567               17734.0            134118.2            7.562772
293           2831                1913.0              7158.4            3.741976
1017         55470             1200399.0           3786786.0            3.154606
```
These warnings could be the result of incorrect crosswalking between subplants, inconsistent or inaccurate reporting of generation data by the plant, or an incorrect monthly allocation of annually-reported EIA-923 data by the EIA. 

### Several plants are missing default gross-to-net conversion factors
```
2023-12-23 11:41:07,486 [WARNING] oge.oge.gross_to_net_generation:93 The following subplants are missing default GTN ratios. Using a default value of 0.97
2023-12-23 11:41:07,486 [WARNING] oge.oge.gross_to_net_generation:96 
          plant_id_eia  subplant_id
17264257         50733            2
17338728         50900            4
23352705         55641            1
```

### Our calculated net generation does not match reported net generation for several plants
Currently, we trust CEMS gross generation data more than EIA-923 net generation data, so in certain instances where these are inconsistent with each other, we use default gross-to-net conversion factors to calculate net generation from the CEMS data, which means that the net generation data will not match what is reported in EIA-923.

```
2023-12-23 11:41:12,857 [WARNING] oge.oge.validation:528 There are 6 plants where calculated annual net generation does not match EIA annual net generation.
2023-12-23 11:41:12,872 [WARNING] oge.oge.validation:531 
              net_generation_mwh_eia  net_generation_mwh_calc  pct_error                                                               gtn_method
plant_id_eia                                                                                                                                     
3406                        845833.6              1228806.261   0.452775  [5_annual_fuel_ratio, 6_default_eia_ratio, 4_annual_plant_shift_factor]
8906                        828822.0              1736788.787   1.095491                               [5_annual_fuel_ratio, 6_default_eia_ratio]
50625                      3444415.8              6045743.480   0.755230                                                    [6_default_eia_ratio]
55075                      1276484.0              1817106.918   0.423525                                                    [5_annual_fuel_ratio]
57865                        41887.8               753792.480  16.995514                               [5_annual_fuel_ratio, 6_default_eia_ratio]
62192                      3019018.0              4452630.770   0.474861                                                    [6_default_eia_ratio]
```

### Some data outputs are missing complete monthly or hourly values
Currently, there are a number of plants for which there may not be 12 reported monthly values or 8760 reported hourly values in the OGE outputs. These likely represent missing input data, but we are working to fill these missing timestamps and ensure that no data is being mistakenly dropped from the data pipeline.

In some of the regional data, the validation checks flagged certain fuel types for which there are more than 8760 reported values in a year. This appears to result when a BA spans multiple time zones. We are investigating this in more detail to correct this issue.

### Some plants have extreme emission factors
There are approximately 30 plants that have extremely high or low calculated emission rates. The extreme low values appear to mostly occur when a clean plant (eg nuclear) has a fossil-based backup generator that ocassionally runs, resulting in a small amount of emissions relative to the total generation. The extreme high values appear to mostly occur at natural gas and petroleum plants with a small amount of net generation. These may represent peaker plants or spinning reserves, which consume a large amount of fuel throughout the year, but generate only a small amount of net electricity. 

### Negative consumed emission rates in HST
For the small balancing area of HST (City of Homestead, FL), there are 62 hours when the calculated consumed emission rate is negative. We are currently investigating this issue, but did not correct this prior to the data release given that it affected only a small number of hours in a small balancing area. 