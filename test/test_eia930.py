import logging

import numpy as np
import pandas as pd
import pytest

from oge.eia930 import (
    add_missing_operations,
    check_bas_outside_us,
    remove_months_with_zero_data,
)
from oge.physics_reconciliation import EIA930Tables, check_eia930_tables

HOURS = pd.to_datetime(["2025-06-01 06:00", "2025-06-01 07:00"], utc=True)


@pytest.fixture
def tables_missing_operations():
    """EIA-930 tables where one BA reports generation and interchange but has no
    operations data.

    "MISO" has operations data. "SIKE" only reports coal and other generation and
    interchange with "MISO" and "SWPP", like the real SIKE in 2025. In the second hour,
    SIKE's generation is blank and only its interchange with MISO is reported.
    """
    operations = pd.DataFrame(
        {
            "datetime_utc": HOURS,
            "balancing_authority_code_eia": "MISO",
            "demand_mwh": [100.0, 110.0],
            "net_generation_mwh": [90.0, 95.0],
            "interchange_mwh": [-10.0, -15.0],
        }
    )
    net_generation_by_energy_source = pd.DataFrame(
        {
            "datetime_utc": [HOURS[0], HOURS[0], HOURS[1], HOURS[1]],
            "balancing_authority_code_eia": "SIKE",
            "generation_energy_source": ["coal", "other", "coal", "other"],
            "net_generation_mwh": [10.0, 5.0, np.nan, np.nan],
        }
    )
    interchange = pd.DataFrame(
        {
            "datetime_utc": [HOURS[0], HOURS[0], HOURS[1], HOURS[1]],
            "balancing_authority_code_eia": "SIKE",
            "balancing_authority_code_adjacent_eia": ["MISO", "SWPP", "MISO", "SWPP"],
            "interchange_mwh": [12.0, 3.0, 4.0, np.nan],
        }
    )
    return EIA930Tables(operations, net_generation_by_energy_source, interchange)


def test_add_missing_operations_derives_operations(tables_missing_operations):
    """A BA with generation and interchange but no operations data should get
    operations data derived from its own reports.

    Without it, the reconciliation treats the BA as outside of the US and replaces its
    data with estimates from its neighbors' interchange, which is what happened to SIKE
    in 2025.
    - net generation is the sum of generation by energy source
    - total interchange is the sum of interchange with each adjacent BA
    - demand is net generation minus total interchange
    - hours where a value can't be computed stay blank
    """
    result = add_missing_operations(tables_missing_operations)
    sike = (
        result.operations[result.operations["balancing_authority_code_eia"] == "SIKE"]
        .set_index("datetime_utc")
        .sort_index()
    )
    # first hour: 10 + 5 MWh of generation, 12 + 3 MWh of exports
    assert sike.loc[HOURS[0], "net_generation_mwh"] == 15.0
    assert sike.loc[HOURS[0], "interchange_mwh"] == 15.0
    assert sike.loc[HOURS[0], "demand_mwh"] == 0.0
    # second hour: no generation, so net generation and demand are blank
    assert np.isnan(sike.loc[HOURS[1], "net_generation_mwh"])
    assert sike.loc[HOURS[1], "interchange_mwh"] == 4.0
    assert np.isnan(sike.loc[HOURS[1], "demand_mwh"])
    # the other tables and the existing operations data are unchanged
    pd.testing.assert_frame_equal(
        result.operations[result.operations["balancing_authority_code_eia"] == "MISO"],
        tables_missing_operations.operations,
    )
    assert result.interchange is tables_missing_operations.interchange
    assert list(result.operations.columns) == list(
        tables_missing_operations.operations.columns
    )


def test_add_missing_operations_no_missing_bas(tables_missing_operations):
    """Tables where every reporting BA has operations data should be returned as is."""
    tables = tables_missing_operations._replace(
        net_generation_by_energy_source=tables_missing_operations.net_generation_by_energy_source.iloc[
            0:0
        ],
        interchange=tables_missing_operations.interchange.iloc[0:0],
    )
    assert add_missing_operations(tables) is tables


def test_check_eia930_tables_warns_about_incomplete_bas(
    tables_missing_operations, caplog
):
    """The reconciliation should warn about BAs that report generation or interchange
    but have no operations data.

    This replaces the "Inconsistent columns" warning from the gridemissions code
    (issue #315), which also fired for every BA outside of the US.
    """
    with caplog.at_level(logging.INFO):
        check_eia930_tables(tables_missing_operations)
    warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "['SIKE']" in warnings[0]

    repaired = add_missing_operations(tables_missing_operations)
    caplog.clear()
    with caplog.at_level(logging.INFO):
        check_eia930_tables(repaired)
    # no warning once SIKE has operations data, and only SWPP is treated as outside of
    # the US because it only appears as an adjacent BA
    assert not [r for r in caplog.records if r.levelno == logging.WARNING]
    assert any("['SWPP']" in r.message for r in caplog.records)


def test_check_bas_outside_us_warns_about_us_bas(caplog):
    """BAs that only appear as adjacent BAs should only be expected outside of the US.

    The reconciliation estimates their data from their neighbors' interchange, so a US
    BA in that position (here SWPP) probably has missing data.
    """
    tables = EIA930Tables(
        operations=pd.DataFrame(
            {
                "datetime_utc": HOURS[:1],
                "balancing_authority_code_eia": ["MISO"],
                "demand_mwh": [100.0],
                "net_generation_mwh": [90.0],
                "interchange_mwh": [-10.0],
            }
        ),
        net_generation_by_energy_source=pd.DataFrame(
            columns=[
                "datetime_utc",
                "balancing_authority_code_eia",
                "generation_energy_source",
                "net_generation_mwh",
            ]
        ),
        interchange=pd.DataFrame(
            {
                "datetime_utc": [HOURS[0], HOURS[0]],
                "balancing_authority_code_eia": ["MISO", "MISO"],
                "balancing_authority_code_adjacent_eia": ["IESO", "SWPP"],
                "interchange_mwh": [-5.0, -5.0],
            }
        ),
    )
    with caplog.at_level(logging.WARNING):
        check_bas_outside_us(tables)
    warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    # SWPP is a US BA; IESO (Ontario) is expected outside of the US
    assert len(warnings) == 1
    assert "['SWPP']" in warnings[0]


def test_remove_months_with_zero_data():
    """BA-fuel-months that are zero in every hour should be removed, and the rest
    kept, including months with only some zero hours.
    """
    hours = pd.date_range("2025-01-31 22:00", periods=4, freq="h")
    report_dates = pd.to_datetime(["2025-01-01"] * 2 + ["2025-02-01"] * 2)
    eia930_data = pd.DataFrame(
        {
            "ba_code": "SWPP",
            "fuel_category_eia930": ["petroleum"] * 4 + ["solar"] * 4,
            "datetime_utc": list(hours) * 2,
            "report_date": list(report_dates) * 2,
            # petroleum: zero all January, zero in one hour of February. Solar:
            # non-zero in January, zero all February.
            "net_generation_mwh_930": [0.0, 0.0, 0.0, 5.0, 3.0, 1.0, 0.0, 0.0],
        }
    )

    result = remove_months_with_zero_data(eia930_data)

    kept = result[["fuel_category_eia930", "report_date"]].drop_duplicates()
    assert sorted(map(tuple, kept.astype(str).to_numpy())) == [
        ("petroleum", "2025-02-01"),
        ("solar", "2025-01-01"),
    ]
    # every hour of the kept months is kept
    assert len(result) == 4
