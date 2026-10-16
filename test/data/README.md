# Test data

## `eia930_raw_*.csv.gz`

Two weeks of real EIA-930 data for every BA, after `eia930.manual_930_adjust`, in the
wide format the original `gridemissions` package reads (one column per series, named
like `EBA.{BA}-ALL.D.H`). Used by `test_physics_reconciliation.py`, which converts them
to `EIA930Tables` for OGE.

- `eia930_raw_2024-11-01_2024-11-14.csv.gz`: covers the first weeks of storage
  reporting in several BAs.
- `eia930_raw_2018-06-24_2018-07-07.csv.gz`: crosses 2018-07-01, the date after which
  the optimization uses generation by energy source.

Both were sliced from data built with the gridemissions-format code that preceded
`oge.physics_reconciliation`. To rebuild one with the current code (the result is
identical for the same PUDL data):

```python
from test.test_physics_reconciliation import convert_tables_to_wide
from oge import eia930

tables = eia930.manual_930_adjust(eia930.load_eia930_balance_data(2018))
wide = convert_tables_to_wide(tables)
wide.loc["2018-06-24":"2018-07-07"].to_csv("eia930_raw_2018-06-24_2018-07-07.csv.gz")
```
