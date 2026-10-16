# Test data

## `eia930_raw_*.csv.gz`

Two weeks of real EIA-930 data for every BA, in the format that `eia930.clean_930`
writes to `eia930_raw.csv` (after `convert_balance_data_to_gridemissions_format` and
`manual_930_adjust`). Used by `test_physics_reconciliation.py`.

- `eia930_raw_2024-11-01_2024-11-14.csv.gz`: sliced from the `eia930_raw.csv` written
  by the 2025 pipeline run. Covers the first weeks of storage reporting in several BAs.
- `eia930_raw_2018-06-24_2018-07-07.csv.gz`: crosses 2018-07-01, the date after which
  the optimization uses generation by source. Built with:

```python
from oge import eia930

df = eia930.convert_balance_data_to_gridemissions_format(2018)
df = eia930.manual_930_adjust(df)
df.loc["2018-06-24":"2018-07-07"].to_csv("eia930_raw_2018-06-24_2018-07-07.csv.gz")
```
