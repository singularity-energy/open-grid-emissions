# Test data

## `eia930_raw_*.parquet`

Two weeks of real EIA-930 data for every BA, after `eia930.manual_930_adjust`, with
one row per series and hour (columns "datetime_utc", "field",
"balancing_authority_code_eia", "other" and "value_mwh", as written by
`melt_series_to_long` in `test_physics_reconciliation.py`). Used by
`test_physics_reconciliation.py`, which converts them to `EIA930Tables` for OGE.

- `eia930_raw_2024-11-01_2024-11-14.parquet`: covers the first weeks of storage
  reporting in several BAs.

To rebuild one from PUDL data:

```python
from oge import eia930
from oge.physics_reconciliation import pivot_eia930_tables_to_series
from test.test_physics_reconciliation import melt_series_to_long

tables = eia930.manual_930_adjust(eia930.load_eia930_balance_data(2024))
series = pivot_eia930_tables_to_series(tables).loc["2024-11-01":"2024-11-14"]
melt_series_to_long(series).to_parquet(
    "eia930_raw_2024-11-01_2024-11-14.parquet", index=False, compression="zstd"
)
```

## `eia930_reconciled_*.parquet`

The reconciled data for each raw fixture, rounded to 0.001 MWh. The test checks that
`reconcile_eia930` still gives these values. After an intended change to the
reconciliation, rewrite them and review the differences before committing:

```python
import tempfile

from test.test_physics_reconciliation import FIXTURE_DATES, write_reconciled_fixture

for dates in FIXTURE_DATES:
    write_reconciled_fixture(dates, tempfile.mkdtemp())
```
