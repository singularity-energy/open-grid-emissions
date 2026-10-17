import logging

import numpy as np
import pandas as pd
import pytest

import oge.consumed as consumed
from oge.consumed import FUEL_TYPE_MAP
from oge.eia930 import (
    add_missing_operations,
    check_bas_outside_us,
    combine_renamed_energy_sources,
    get_columns,
    prepare_eia930_for_shaping,
    remove_months_with_zero_data,
)
from oge.physics_reconciliation import (
    ENERGY_SOURCES,
    NET_GENERATION_BY_ENERGY_SOURCE,
    EIA930Tables,
    check_eia930_tables,
)

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


@pytest.fixture
def pudl_generation_across_rename():
    """PUDL generation by energy source for two BAs, in the hour before and the hour
    after PUDL renamed hydro to "hydro_excluding_pumped_storage" on 2024-07-01.

    As in PUDL, there is a row for every BA, hour and energy source, blank when
    nothing is reported:
        - hydro is reported under "hydro" before the rename and under
          "hydro_excluding_pumped_storage" after it
        - ISNE reports pumped storage after the rename; SOCO never does
        - neither BA reports "unknown", an original energy source
    """
    hours = pd.to_datetime(["2024-07-01 03:00", "2024-07-01 04:00"], utc=True)
    values = {
        # energy source: (ISNE before, ISNE after, SOCO before, SOCO after)
        "hydro": (100.0, np.nan, 300.0, np.nan),
        "hydro_excluding_pumped_storage": (np.nan, 110.0, np.nan, 310.0),
        "pumped_storage": (np.nan, -20.0, np.nan, np.nan),
        "unknown": (np.nan, np.nan, np.nan, np.nan),
    }
    rows = []
    for energy_source, (isne_0, isne_1, soco_0, soco_1) in values.items():
        for ba, hour, value in [
            ("ISNE", hours[0], isne_0),
            ("ISNE", hours[1], isne_1),
            ("SOCO", hours[0], soco_0),
            ("SOCO", hours[1], soco_1),
        ]:
            rows.append(
                {
                    "datetime_utc": hour,
                    "balancing_authority_code_eia": ba,
                    "generation_energy_source": energy_source,
                    "net_generation_mwh": value,
                }
            )
    data = pd.DataFrame(rows)
    data["generation_energy_source"] = data["generation_energy_source"].astype(
        "category"
    )
    return data


def test_combine_renamed_energy_sources(pudl_generation_across_rename):
    """Renamed energy sources should be combined with their original names, and new
    energy sources should be kept only for the BAs that report them.

    The new energy sources are reconciled as their own series, but a series that a BA
    never reports would only add zeros to the reconciliation.
    """
    result = combine_renamed_energy_sources(pudl_generation_across_rename)
    values = result.set_index(
        ["balancing_authority_code_eia", "generation_energy_source", "datetime_utc"]
    )["net_generation_mwh"]

    # hydro is one series across the rename
    assert values.loc[("ISNE", "hydro")].to_list() == [100.0, 110.0]
    assert values.loc[("SOCO", "hydro")].to_list() == [300.0, 310.0]
    assert "hydro_excluding_pumped_storage" not in set(
        result["generation_energy_source"]
    )
    # pumped storage is its own series for ISNE, zero before it was reported, and
    # isn't added for SOCO
    assert values.loc[("ISNE", "pumped_storage")].to_list() == [0.0, -20.0]
    assert ("SOCO", "pumped_storage") not in values.index.droplevel(2)
    # original energy sources are kept for every BA, even if never reported
    assert values.loc[("ISNE", "unknown")].to_list() == [0.0, 0.0]
    assert values.loc[("SOCO", "unknown")].to_list() == [0.0, 0.0]


def test_combine_renamed_energy_sources_warns_about_both_names(
    pudl_generation_across_rename, caplog
):
    """A BA that reports an energy source under both its original and its renamed
    name in the same hour should be reported, and its values added together.

    This doesn't happen in PUDL data from 2023-10 to 2026-01, but a value reported
    twice would otherwise be counted twice without notice.
    """
    data = pudl_generation_across_rename.copy()
    is_isne_hydro_after = (
        (data["balancing_authority_code_eia"] == "ISNE")
        & (data["generation_energy_source"] == "hydro")
        & (data["datetime_utc"] == data["datetime_utc"].max())
    )
    data.loc[is_isne_hydro_after, "net_generation_mwh"] = 5.0

    with caplog.at_level(logging.WARNING):
        result = combine_renamed_energy_sources(data)

    warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "ISNE" in warnings[0]
    # the two values are added
    isne_hydro = result[
        (result["balancing_authority_code_eia"] == "ISNE")
        & (result["generation_energy_source"] == "hydro")
    ]
    assert isne_hydro["net_generation_mwh"].to_list() == [100.0, 115.0]


def test_prepare_eia930_for_shaping_adds_new_energy_sources_before_clipping():
    """New energy sources should be added to the fuel category that included them
    before, and only then should negative values be set to zero.

    Battery charging (negative battery storage) then offsets other generation in the
    "other" category, as it did before the BA reported batteries separately, rather
    than being dropped. The scenario is one ERCO hour in which batteries charge more
    than "other" generates, and one in which both are positive.
    """
    hours = pd.to_datetime(["2025-06-01 18:00", "2025-06-01 19:00"], utc=True)
    reconciled = pd.DataFrame(
        {
            "datetime_utc": list(hours) * 4,
            "balancing_authority_code_eia": "ERCO",
            "generation_energy_source": ["other"] * 2
            + ["battery_storage"] * 2
            + ["geothermal"] * 2
            + ["solar"] * 2,
            "net_generation_mwh": [100.0, 100.0, -150.0, 50.0, 0.0, 10.0, 7.0, 8.0],
        }
    )

    result = prepare_eia930_for_shaping(reconciled, 2025).set_index(
        ["fuel_category_eia930", "datetime_utc"]
    )["net_generation_mwh_930"]

    # other + battery storage + geothermal, with the net load in the first hour set to
    # zero after adding them up
    assert result.loc["other"].to_list() == [0.0, 160.0]
    # other categories are unchanged
    assert result.loc["solar"].to_list() == [7.0, 8.0]


def test_get_columns_includes_new_energy_sources():
    """The timestamp adjustments should shift every energy source of a BA, including
    the ones EIA-930 introduced in 2024.

    `adjust_series_timestamps` shifts the generation of some BAs (e.g. PJM, CISO) by
    selecting their series with `get_columns`.
    """
    columns = pd.MultiIndex.from_tuples(
        [
            (NET_GENERATION_BY_ENERGY_SOURCE, "CISO", "solar"),
            (NET_GENERATION_BY_ENERGY_SOURCE, "CISO", "geothermal"),
            (NET_GENERATION_BY_ENERGY_SOURCE, "CISO", "battery_storage"),
            (NET_GENERATION_BY_ENERGY_SOURCE, "ERCO", "battery_storage"),
            ("net_generation", "CISO", ""),
            ("demand", "CISO", ""),
        ]
    )

    result = get_columns("CISO", columns)

    assert set(result) == {
        (NET_GENERATION_BY_ENERGY_SOURCE, "CISO", "solar"),
        (NET_GENERATION_BY_ENERGY_SOURCE, "CISO", "geothermal"),
        (NET_GENERATION_BY_ENERGY_SOURCE, "CISO", "battery_storage"),
        ("net_generation", "CISO", ""),
        ("demand", "CISO", ""),
    }


def test_fuel_type_map_covers_every_energy_source():
    """Every EIA-930 energy source should have an emission rate category for BAs
    outside of the US.

    `consumed.get_average_emission_factors` looks up each energy source in
    `ENERGY_SOURCES`, so a missing one would raise a KeyError in step 19.
    """
    assert set(ENERGY_SOURCES) <= set(FUEL_TYPE_MAP)


@pytest.mark.parametrize("has_storage_row", [True, False])
def test_get_average_emission_factors_uses_zero_for_storage(
    has_storage_row, tmp_path, monkeypatch
):
    """EIA-930 storage energy sources should get the US storage fleet's emission rate,
    which is zero, and zero if the US fleet data has no storage row.

    These rates apply to BAs outside of the US. Storage has no direct emissions, and
    emissions from the electricity used to charge it aren't tracked through storage
    yet, so the average rate of all generation would overstate them.
    """
    columns = [
        consumed.get_pollutant_rate_column_name(pol, adjustment, generated=True)
        for pol in consumed.POLLUTANTS
        for adjustment in consumed.ADJUSTMENTS
    ]
    rates = {"total": 800.0, "geothermal": 50.0, "solar": 0.0, "storage": -0.0}
    if not has_storage_row:
        del rates["storage"]
    us_fleet = pd.DataFrame(
        {column: list(rates.values()) for column in columns},
        index=pd.Index(list(rates), name="fuel_category"),
    )
    folder = tmp_path / "2025" / "power_sector_data" / "annual" / "us_units"
    folder.mkdir(parents=True)
    us_fleet.to_csv(folder / "US.csv")
    monkeypatch.setattr(consumed, "results_folder", lambda rel: str(tmp_path / rel))

    factors = consumed.get_average_emission_factors("2025")

    co2 = factors["CO2"]["for_electricity"]
    # storage energy sources have a zero rate
    for energy_source in [
        "battery_storage",
        "unknown_energy_storage",
        "other_energy_storage",
        "pumped_storage",
    ]:
        assert co2[energy_source] == 0.0
    # the other new energy sources use the rate of their own fuel
    assert co2["geothermal"] == 50.0
    assert co2["solar_w_integrated_battery_storage"] == 0.0
    # energy sources without their own fleet use the total rate
    assert co2["other"] == 800.0
