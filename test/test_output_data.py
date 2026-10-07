from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

import oge.output_data as output_data
from oge.column_checks import DATA_COLUMNS, STORAGE_DATA_COLUMNS


def capture_outputs(fleet_data, **kwargs):
    """Runs `write_power_sector_results` and returns the tables passed to
    `output_to_results`, keyed by (subfolder, file_name), along with the value of
    `skip_outputs` used for each."""
    outputs = {}

    def fake_output_to_results(
        df, year, file_name, subfolder, path_prefix, skip_outputs, **_
    ):
        outputs[(subfolder, file_name)] = (df.copy(), skip_outputs)

    with (
        patch.object(output_data, "output_to_results", fake_output_to_results),
        patch.object(output_data.load_data, "ba_timezone", return_value="Etc/GMT+5"),
    ):
        output_data.write_power_sector_results(fleet_data, 2023, "2023/", **kwargs)
    return outputs


@pytest.fixture
def monthly_fleet_data():
    """Monthly fleet data for two BAs over two months.

    BA1 and BA2 both have natural gas and solar, with every data column equal to 1 for
    BA1 and 2 for BA2. Only BA2 has storage, with negative net generation and non-zero
    charging and discharging. Storage columns are blank for all other fuels. A third
    set of records has no BA and large values (100), so any test totals that include
    these records will be obviously wrong.
    """
    records = []
    for ba, scale in [("BA1", 1.0), ("BA2", 2.0), (np.nan, 100.0)]:
        for month in ["2023-01-01", "2023-02-01"]:
            for fuel in ["natural_gas", "solar"]:
                record = {
                    "ba_code": ba,
                    "fuel_category": fuel,
                    "report_date": pd.Timestamp(month),
                }
                record.update({col: scale for col in DATA_COLUMNS})
                record["storage_charge_mwh"] = np.nan
                record["storage_discharge_mwh"] = np.nan
                records.append(record)
        if ba == "BA2":
            for month in ["2023-01-01", "2023-02-01"]:
                record = {
                    "ba_code": ba,
                    "fuel_category": "storage",
                    "report_date": pd.Timestamp(month),
                }
                record.update({col: 0.0 for col in DATA_COLUMNS})
                record["net_generation_mwh"] = -1.0
                record["storage_charge_mwh"] = 5.0
                record["storage_discharge_mwh"] = 4.0
                records.append(record)
    return pd.DataFrame(records)


@pytest.mark.parametrize("agg_level", ["monthly", "annual"])
def test_national_results_equal_sum_of_ba_results(monthly_fleet_data, agg_level):
    """The monthly and annual national files should equal the sum of the BA files.

    Checks every data and storage column for each fuel category (and month, for the
    monthly files), including the total row, and checks that the national file has the
    same columns in the same order as the BA files. Because both are written from the
    same fleet data, any difference would mean the national aggregation is using
    different data or logic than the BA files.
    """
    outputs = capture_outputs(
        monthly_fleet_data,
        skip_outputs=False,
        include_hourly=False,
        include_monthly=True,
        include_annual=True,
    )
    subfolder = f"power_sector_data/{agg_level}/"
    keys = ["fuel_category"] + (["report_date"] if agg_level == "monthly" else [])
    columns = DATA_COLUMNS + STORAGE_DATA_COLUMNS

    ba_sum = (
        pd.concat([outputs[(subfolder, ba)][0] for ba in ["BA1", "BA2"]])
        .groupby(keys)[columns]
        .sum(min_count=1)
        .sort_index()
    )
    national = outputs[(subfolder, "US")][0].set_index(keys)[columns].sort_index()

    pd.testing.assert_frame_equal(national, ba_sum)
    # the national file has the same columns, in the same order, as the BA files
    assert list(outputs[(subfolder, "US")][0].columns) == list(
        outputs[(subfolder, "BA1")][0].columns
    )


def test_national_annual_results(monthly_fleet_data):
    """The annual national file should only include data assigned to a BA, and should
    handle storage data the same way as the BA files.

    - Records that are not assigned to a BA are excluded, matching the BA files.
    - Storage columns are blank for non-storage fuels (rather than being filled with
      zero), are summed for the storage fuel, and are included in the total row.
    - The total row sums all fuel categories, including negative net generation from
      storage.
    """
    outputs = capture_outputs(
        monthly_fleet_data,
        skip_outputs=False,
        include_hourly=False,
        include_monthly=True,
        include_annual=True,
    )
    national = outputs[("power_sector_data/annual/", "US")][0].set_index(
        "fuel_category"
    )

    # records without a BA are excluded: 2 months x (BA1 + BA2)
    assert national.loc["natural_gas", "net_generation_mwh"] == 6.0
    # storage columns are blank for non-storage fuels, and included in the total
    assert national.loc["natural_gas", STORAGE_DATA_COLUMNS].isna().all()
    assert national.loc["storage", "storage_charge_mwh"] == 10.0
    assert national.loc["total", "storage_charge_mwh"] == 10.0
    assert national.loc["total", "storage_discharge_mwh"] == 8.0
    # the total row sums all fuels, including negative storage net generation
    assert national.loc["total", "net_generation_mwh"] == 6.0 + 6.0 - 2.0


def test_national_annual_written_when_skipping_outputs(monthly_fleet_data):
    """The annual national file should be written even if `skip_outputs` is True.

    `consumed.py` reads the annual national file to estimate emission rates for BAs
    outside the US, so it must always exist. All other files, including the monthly
    national file and the BA files, should still follow `skip_outputs`.
    """
    outputs = capture_outputs(
        monthly_fleet_data,
        skip_outputs=True,
        include_hourly=False,
        include_monthly=True,
        include_annual=True,
    )

    # only the annual national file is written
    assert outputs[("power_sector_data/annual/", "US")][1] is False
    assert outputs[("power_sector_data/monthly/", "US")][1] is True
    assert outputs[("power_sector_data/annual/", "BA1")][1] is True


def test_national_hourly_results():
    """The hourly national file should equal the sum of the hourly BA files, even when
    the BAs are in different timezones.

    Each hourly BA file covers that BA's local year, so BAs in different timezones
    cover different UTC hours. Here, BA1's data starts one hour earlier than BA2's.
    The national data should:
    - be aggregated by UTC hour, covering the hours of both BAs, with the first and
      last hours only including data from one BA,
    - equal the sum of the BA data in each hour, for each fuel and the total row, and
    - not include a local datetime column, since the US spans multiple timezones,
      while the BA files still include one.
    """
    datetimes = pd.date_range("2023-01-01 05:00", periods=3, freq="h", tz="UTC")
    records = []
    for ba, hours in [("BA1", datetimes[:2]), ("BA2", datetimes[1:])]:
        for hour in hours:
            for fuel in ["natural_gas", "wind"]:
                record = {
                    "ba_code": ba,
                    "fuel_category": fuel,
                    "datetime_utc": hour,
                    "report_date": pd.Timestamp("2023-01-01"),
                }
                record.update({col: 1.0 for col in DATA_COLUMNS})
                records.append(record)
    outputs = capture_outputs(
        pd.DataFrame(records),
        skip_outputs=False,
        include_hourly=True,
        include_monthly=False,
        include_annual=False,
    )
    national = outputs[("power_sector_data/hourly/", "US")][0]

    # there is no local datetime for the national data
    assert "datetime_local" not in national.columns
    assert "datetime_local" in outputs[("power_sector_data/hourly/", "BA1")][0]
    # the national data equals the sum of the BA data in each hour
    total = national[national["fuel_category"] == "total"].set_index("datetime_utc")
    assert total["net_generation_mwh"].tolist() == [2.0, 4.0, 2.0]
    wind = national[national["fuel_category"] == "wind"].set_index("datetime_utc")
    assert wind["net_generation_mwh"].tolist() == [1.0, 2.0, 1.0]
