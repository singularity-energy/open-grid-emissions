# Test data

## `eia930_raw_*.parquet`

Two weeks of real EIA-930 data for every BA, after `eia930.manual_930_adjust`, with
one row per series and hour (columns "datetime_utc", "field",
"balancing_authority_code_eia", "other" and "value_mwh", as written by
`melt_series_to_long` in `test_physics_reconciliation.py`). Used by
`test_physics_reconciliation.py`, which converts them to `EIA930Tables` for OGE.

- `eia930_raw_2024-11-01_2024-11-14.parquet`: covers the first weeks of storage
  reporting in several BAs, with the energy sources that EIA-930 introduced in 2024
  (e.g. battery storage in ERCO, IPCO and NEVP, pumped storage in ISNE) as their own
  series. Rebuilt from PUDL in October 2026.
- `eia930_raw_2018-06-24_2018-07-07.parquet`: crosses 2018-07-01, the date after which
  the optimization uses generation by energy source.

To rebuild one from PUDL data:

```python
from oge import eia930
from oge.physics_reconciliation import pivot_eia930_tables_to_series
from test.test_physics_reconciliation import melt_series_to_long

tables = eia930.manual_930_adjust(eia930.load_eia930_balance_data(2018))
series = pivot_eia930_tables_to_series(tables).loc["2018-06-24":"2018-07-07"]
melt_series_to_long(series).to_parquet(
    "eia930_raw_2018-06-24_2018-07-07.parquet", index=False, compression="zstd"
)
```

## `eia930_reconciled_*.parquet`

The reconciled data for each raw fixture, rounded to 0.001 MWh, with a "stage" column
("reconciled_without_energy_sources" or "reconciled"). The test checks that
`reconcile_eia930` still gives these values. After an intended change to the
reconciliation, rewrite them and review the differences before committing:

```python
import tempfile

from test.test_physics_reconciliation import FIXTURE_DATES, write_reconciled_fixture

for dates in FIXTURE_DATES:
    write_reconciled_fixture(dates, tempfile.mkdtemp())
```
