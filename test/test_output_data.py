from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

import oge.output_data as output_data
from oge.column_checks import DATA_COLUMNS, HG_DATA_COLUMNS, STORAGE_DATA_COLUMNS


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


def test_residual_mix_annual_results(monthly_fleet_data):
    """The residual row should sum only the residual mix fuel categories.

    In the monthly fleet data, natural gas is the only residual mix fuel, so the
    residual row should equal the natural gas row, excluding solar and storage, with blank storage
    columns and the same emission rates as natural gas.
    """
    outputs = capture_outputs(
        monthly_fleet_data,
        skip_outputs=False,
        include_hourly=False,
        include_monthly=True,
        include_annual=True,
    )
    for ba in ["BA1", "BA2", "US"]:
        table = outputs[("power_sector_data/annual/", ba)][0].set_index("fuel_category")
        pd.testing.assert_series_equal(
            table.loc[
                "residual", DATA_COLUMNS + output_data.GENERATED_EMISSION_RATE_COLS
            ],
            table.loc[
                "natural_gas", DATA_COLUMNS + output_data.GENERATED_EMISSION_RATE_COLS
            ],
            check_names=False,
        )
        assert table.loc["residual", STORAGE_DATA_COLUMNS].isna().all()


def test_residual_mix_hourly_results():
    """Every hourly BA file should have a residual row for every hour.

    BA1 has coal, natural gas, waste, and wind, so its residual row should sum coal and
    natural gas, excluding waste and wind. BA2 only has wind, so its residual row
    should be zero, with zero emission rates, rather than missing.
    """
    datetimes = pd.date_range("2023-01-01 05:00", periods=2, freq="h", tz="UTC")
    records = []
    for ba, fuels in [
        ("BA1", ["coal", "natural_gas", "waste", "wind"]),
        ("BA2", ["wind"]),
    ]:
        for hour in datetimes:
            for fuel in fuels:
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

    ba1 = outputs[("power_sector_data/hourly/", "BA1")][0]
    ba1_residual = ba1[ba1["fuel_category"] == "residual"]
    assert ba1_residual["net_generation_mwh"].tolist() == [2.0, 2.0]
    assert ba1_residual["co2_mass_lb_for_electricity"].tolist() == [2.0, 2.0]
    assert ba1_residual["generated_co2_rate_lb_per_mwh_for_electricity"].tolist() == [
        1.0,
        1.0,
    ]

    ba2 = outputs[("power_sector_data/hourly/", "BA2")][0]
    ba2_residual = ba2[ba2["fuel_category"] == "residual"]
    assert len(ba2_residual) == len(datetimes)
    assert (ba2_residual[DATA_COLUMNS] == 0).all().all()
    assert (ba2_residual[output_data.GENERATED_EMISSION_RATE_COLS] == 0).all().all()

    us = outputs[("power_sector_data/hourly/", "US")][0]
    us_residual = us[us["fuel_category"] == "residual"]
    assert us_residual["net_generation_mwh"].tolist() == [2.0, 2.0]


def test_residual_row_without_residual_generation(monthly_fleet_data):
    """A BA with no fossil generation should still have a residual row of zeros in the
    monthly and annual files."""
    no_residual = monthly_fleet_data[
        monthly_fleet_data["fuel_category"] != "natural_gas"
    ]
    outputs = capture_outputs(
        no_residual,
        skip_outputs=False,
        include_hourly=False,
        include_monthly=True,
        include_annual=True,
    )
    for agg_level, n_rows in [("monthly", 2), ("annual", 1)]:
        table = outputs[(f"power_sector_data/{agg_level}/", "BA1")][0]
        residual = table[table["fuel_category"] == "residual"]
        assert len(residual) == n_rows
        assert (residual[DATA_COLUMNS] == 0).all().all()


@pytest.mark.parametrize(
    "resolution, keys",
    [("hourly", ["datetime_utc"]), ("monthly", ["report_date"]), ("annual", [])],
)
def test_fuel_category_sort_order(resolution, keys):
    """Fuel categories should be sorted alphabetically, followed by the residual and
    total rows, and by date within each fuel category."""
    datetimes = pd.date_range("2023-01-01 05:00", periods=2, freq="h", tz="UTC")
    records = []
    # records are deliberately out of order
    for hour in reversed(datetimes):
        for fuel in ["wind", "natural_gas", "solar", "coal"]:
            record = {
                "ba_code": "BA1",
                "fuel_category": fuel,
                "datetime_utc": hour,
                "report_date": pd.Timestamp("2023-01-01"),
            }
            record.update({col: 1.0 for col in DATA_COLUMNS})
            records.append(record)
    fleet_data = pd.DataFrame(records)
    if resolution != "hourly":
        fleet_data = fleet_data.drop(columns="datetime_utc")
    outputs = capture_outputs(
        fleet_data,
        skip_outputs=False,
        include_hourly=resolution == "hourly",
        include_monthly=resolution == "monthly",
        include_annual=resolution == "annual",
    )

    table = outputs[(f"power_sector_data/{resolution}/", "BA1")][0]
    expected_order = ["coal", "natural_gas", "solar", "wind", "residual", "total"]
    assert table["fuel_category"].unique().tolist() == expected_order
    pd.testing.assert_frame_equal(
        table,
        table.sort_values(keys)
        .sort_values(
            "fuel_category", key=lambda x: x.map(expected_order.index), kind="stable"
        )
        .reset_index(drop=True),
    )


@pytest.mark.parametrize("resolution", ["monthly", "annual"])
@pytest.mark.parametrize("plant_part", ["plant", "subplant"])
def test_plant_data_does_not_sum_subplant_id(plant_part, resolution):
    """Plant-level results should not have a subplant_id column (which would otherwise
    be the sum of the plant's subplant IDs), while subplant-level results keep it."""
    records = []
    for subplant_id in [1, 2]:
        for month in ["2024-01-01", "2024-02-01"]:
            record = {
                "plant_id_eia": 1,
                "subplant_id": subplant_id,
                "report_date": pd.Timestamp(month),
            }
            record.update({col: 1.0 for col in DATA_COLUMNS})
            record.update({col: np.nan for col in STORAGE_DATA_COLUMNS})
            record.update({col: np.nan for col in HG_DATA_COLUMNS})
            records.append(record)
    monthly_subplant_data = pd.DataFrame(records).astype({"subplant_id": "Int32"})
    plant_attributes = pd.DataFrame(
        {
            "plant_id_eia": [1],
            "plant_name_eia": ["Plant 1"],
            "fuel_category": ["natural_gas"],
            "plant_primary_fuel": ["NG"],
            "capacity_mw": [100.0],
            "ba_code": ["BA1"],
            "city": ["City"],
            "county": ["County"],
            "state": ["ST"],
        }
    )

    outputs = {}

    def fake_output_to_results(df, year, file_name, subfolder, **_):
        outputs[(subfolder, file_name)] = df.copy()

    with (
        patch.object(output_data, "output_to_results", fake_output_to_results),
        patch.object(output_data.validation, "check_for_anomalous_co2_factors"),
    ):
        output_data.write_plant_data_to_results(
            monthly_subplant_data,
            2024,
            path_prefix="2024/",
            plant_part=plant_part,
            resolution=resolution,
            skip_outputs=False,
            plant_attributes=plant_attributes,
        )

    table = outputs[(f"plant_data/{resolution}/", f"{plant_part}_data")]
    months_per_row = 2 if resolution == "annual" else 1
    if plant_part == "plant":
        assert "subplant_id" not in table.columns
        expected_value = 2.0 * months_per_row
    else:
        assert sorted(table["subplant_id"].unique()) == [1, 2]
        expected_value = 1.0 * months_per_row
    assert (table[DATA_COLUMNS] == expected_value).all().all()
    assert table[STORAGE_DATA_COLUMNS].isna().all().all()


def test_capacity_factor_uses_hours_in_year():
    """Capacity factors should divide net generation by nameplate capacity times the
    hours in the data year, including the extra day of a leap year, and be negative
    where net generation is negative.

    Capacity factors are blank where capacity is zero or blank, rather than infinite.
    """
    generation = pd.Series([8784.0, 8760.0, -100.0, 10.0, 10.0])
    capacity = pd.Series([1.0, 1.0, 1.0, 0.0, np.nan])
    leap_year = output_data.calculate_capacity_factor(generation, capacity, 2024)
    assert leap_year[0] == 1.0
    assert leap_year[2] == pytest.approx(-100 / 8784)
    # zero or blank capacity
    assert leap_year[[3, 4]].isna().all()
    # a year that is not a leap year has 8,760 hours
    assert output_data.calculate_capacity_factor(generation, capacity, 2023)[1] == 1.0


def test_nonbaseload_factor_follows_egrid():
    """The nonbaseload factor should follow eGRID's piecewise linear function of the
    capacity factor and its exclusions.

    - Below a capacity factor of 0.2 a plant is fully nonbaseload (1), above 0.8 it is
      fully baseload (0), and in between the factor is -5/3 * CF + 4/3.
    - Plants whose capacity factor is zero or blank have no factor, as in eGRID, since
      they have no nonbaseload generation. Negative capacity factors also have none,
      although the plant results set them to zero before this is calculated.
    - Plants whose primary fuel is not dispatched to follow load (e.g. nuclear, solar,
      hydro) have no factor, but biomass plants do, as in eGRID's code.
    """
    capacity_factor = pd.Series([0.1, 0.2, 0.5, 0.8, 0.95, 0.0, -0.1, np.nan, 0.5, 0.1])
    primary_fuel = pd.Series(
        ["NG", "NG", "BIT", "NG", "NG", "NG", "NG", "NG", "NUC", "WDS"]
    )
    factor = output_data.calculate_nonbaseload_factor(capacity_factor, primary_fuel)
    assert factor[0] == 1.0
    assert factor[1] == pytest.approx(1.0)
    assert factor[2] == pytest.approx(-5 / 3 * 0.5 + 4 / 3)
    assert factor[3] == pytest.approx(0.0)
    assert factor[4] == 0.0
    # no factor for zero, negative or blank capacity factors
    assert factor[[5, 6, 7]].isna().all()
    # nuclear is excluded, biomass is not
    assert pd.isna(factor[8])
    assert factor[9] == 1.0


def test_annual_plant_data_has_capacity_and_nonbaseload_factors():
    """The annual plant results should include each plant's capacity factor and
    nonbaseload factor next to its capacity, without the plant primary fuel used to
    calculate the nonbaseload factor.

    The eGRID-format file reads both factors from these results, so they must be
    calculated from the plant's annual generation and capacity. Plant 1 generates half
    of what it could at full capacity in a leap year. Plant 2 has negative net
    generation, so, as in eGRID's plant sheet, its capacity factor is zero rather than
    negative, and it has no nonbaseload factor.
    """
    records = []
    for plant_id, generation in [(1, 100 * 8784 / 4), (2, -50.0)]:
        for month in ["2024-01-01", "2024-02-01"]:
            record = {
                "plant_id_eia": plant_id,
                "subplant_id": 1,
                "report_date": pd.Timestamp(month),
            }
            record.update({col: 0.0 for col in DATA_COLUMNS})
            record.update({col: np.nan for col in STORAGE_DATA_COLUMNS})
            record.update({col: np.nan for col in HG_DATA_COLUMNS})
            record["net_generation_mwh"] = generation
            records.append(record)
    monthly_subplant_data = pd.DataFrame(records).astype({"subplant_id": "Int32"})
    plant_attributes = pd.DataFrame(
        {
            "plant_id_eia": [1, 2],
            "plant_name_eia": ["Plant 1", "Plant 2"],
            "fuel_category": ["natural_gas", "natural_gas"],
            "plant_primary_fuel": ["NG", "NG"],
            "capacity_mw": [100.0, 100.0],
            "ba_code": ["BA1", "BA1"],
            "city": ["City", "City"],
            "county": ["County", "County"],
            "state": ["ST", "ST"],
        }
    )
    outputs = {}

    def fake_output_to_results(df, year, file_name, subfolder, **_):
        outputs[(subfolder, file_name)] = df.copy()

    with (
        patch.object(output_data, "output_to_results", fake_output_to_results),
        patch.object(output_data.validation, "check_for_anomalous_co2_factors"),
    ):
        output_data.write_plant_data_to_results(
            monthly_subplant_data,
            2024,
            path_prefix="2024/",
            plant_part="plant",
            resolution="annual",
            skip_outputs=False,
            plant_attributes=plant_attributes,
        )

    table = outputs[("plant_data/annual/", "plant_data")]
    columns = list(table.columns)
    assert columns[
        columns.index("capacity_mw") + 1 : columns.index("capacity_mw") + 3
    ] == [
        "capacity_factor",
        "nonbaseload_factor",
    ]
    assert "plant_primary_fuel" not in table
    # half of the year's generation at full capacity
    assert table.loc[0, "capacity_factor"] == pytest.approx(0.5)
    assert table.loc[0, "nonbaseload_factor"] == pytest.approx(-5 / 3 * 0.5 + 4 / 3)
    # negative generation gives a zero capacity factor and no nonbaseload factor
    assert table.loc[1, "capacity_factor"] == 0
    assert pd.isna(table.loc[1, "nonbaseload_factor"])


def test_plant_factors_are_rounded_to_fixed_decimals():
    """Capacity and nonbaseload factors should keep four decimals however small the
    median factor is, since they are multiplied by generation and emissions.
    """
    table = pd.DataFrame(
        {
            "capacity_factor": [0.123456, 0.9],
            "nonbaseload_factor": [0.987654, 0.0],
        }
    )
    rounded = output_data.round_table(table)
    assert rounded["capacity_factor"].tolist() == [0.1235, 0.9]
    assert rounded["nonbaseload_factor"].tolist() == [0.9877, 0.0]
