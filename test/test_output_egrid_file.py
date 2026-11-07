import pytest

from oge.output_egrid_file import (
    EGRID_CATEGORIES,
    EGRID_SHEETS,
    get_egrid_columns,
    replace_units,
)

# number of columns in each sheet of the eGRID2023 U.S. units and metric files
US_COLUMN_COUNTS = {
    "UNT": 33,
    "GEN": 17,
    "PLNT": 150,
    "ST": 169,
    "BA": 169,
    "SRL": 169,
    "NRL": 169,
    "US": 167,
    "GGL": 6,
}
METRIC_COLUMN_COUNTS = {
    "UNT": 33,
    "GEN": 19,
    "PLNT": 187,
    "ST": 255,
    "BA": 255,
    "SRL": 255,
    "NRL": 255,
    "US": 253,
    "GGL": 9,
}


def headers(columns):
    """Returns the header of each column, keyed by its code."""
    return {code: header for code, header, _ in columns}


@pytest.mark.parametrize("sheet", EGRID_SHEETS)
@pytest.mark.parametrize("metric", [False, True])
def test_column_counts_and_unique_codes(sheet, metric):
    columns = get_egrid_columns(sheet, 2023, metric=metric)
    expected = (METRIC_COLUMN_COUNTS if metric else US_COLUMN_COUNTS)[sheet]
    assert len(columns) == expected
    codes = [code for code, _, _ in columns]
    assert len(set(codes)) == len(codes)


def test_sequence_number_code_uses_data_year():
    assert get_egrid_columns("PLNT", 2025)[0] == (
        "SEQPLT25",
        "Plant file sequence number",
        None,
    )
    assert get_egrid_columns("UNT", 2023)[0][0] == "SEQUNT23"


def test_regional_sheets_share_data_columns():
    us = get_egrid_columns("US", 2023)
    assert us[:3] == [
        ("YEAR", "Data year", None),
        ("USNAMEPCAP", "U.S. nameplate capacity (MW)", None),
        ("USHTIAN", "U.S. annual heat input from combustion (MMBtu)", 1),
    ]
    ba = headers(get_egrid_columns("BA", 2023))
    assert ba["BACODE"] == "Balancing authority code"
    assert ba["BANGENAN"] == "BA annual net generation (MWh)"
    srl = headers(get_egrid_columns("SRL", 2023))
    assert srl["SRCO2RTA"] == (
        "eGRID subregion annual CO2 total output emission rate (lb/MWh)"
    )


def test_metric_columns_convert_units_and_add_gj_columns():
    metric = get_egrid_columns("BA", 2023, metric=True)
    codes = [code for code, _, _ in metric]
    headers_by_code = headers(metric)
    # columns in MWh keep their units and are followed by a column in GJ
    assert codes[codes.index("BANGENAN") + 1] == "BANGENAN2"
    assert headers_by_code["BANGENAN"] == "BA annual net generation (MWh)"
    assert headers_by_code["BANGENAN2"] == "BA annual net generation (GJ)"
    # rates per MWh are converted to kg and followed by a rate per GJ
    assert (
        headers_by_code["BANOXRTA"]
        == "BA annual NOx total output emission rate (kg/MWh)"
    )
    assert (
        headers_by_code["BANOXRTA2"]
        == "BA annual NOx total output emission rate (kg/GJ)"
    )
    # other units are converted without adding a column
    assert headers_by_code["BACO2AN"] == "BA annual CO2 emissions (metric tons)"
    assert headers_by_code["BAHTIAN"] == "BA annual heat input from combustion (GJ)"
    assert headers_by_code["BANOXRA"] == "BA annual NOx input emission rate (kg/GJ)"
    assert headers_by_code["BAHGAN"] == "BA annual Hg emissions (kg)"
    assert "BAHTIAN2" not in headers_by_code
    # units without a metric equivalent are unchanged
    assert headers_by_code["BANAMEPCAP"] == "BA nameplate capacity (MW)"
    assert headers_by_code["BACLPR"] == "BA coal generation percent (resource mix)"


def test_headers_fix_egrid_inconsistencies():
    plant = headers(get_egrid_columns("PLNT", 2023))
    state = headers(get_egrid_columns("ST", 2023))
    # the other unknown/purchased total is labeled renewable in eGRID's plant sheet
    # and nonrenewable in its regional sheets
    assert plant["PLGENATO"] == (
        "Plant annual total nonrenewables other unknown/purchased net generation (MWh)"
    )
    assert state["STGENATO"] == (
        "State annual total nonrenewables other unknown/purchased net generation (MWh)"
    )
    assert state["STNBNOX"] == (
        "State annual NOx nonbaseload output emission rate (lb/MWh)"
    )
    assert plant["CHPCO2E"] == (
        "Plant annual CO2 equivalent emissions CHP adjustment value (tons)"
    )
    assert EGRID_CATEGORIES[14].endswith("per MMBtu)")


def test_replace_units_only_replaces_listed_units():
    units = {"MWh": "GJ"}
    assert replace_units("Direct use (MWh)", units) == "Direct use (GJ)"
    assert replace_units("Total disposition (MWh) without exports", units) == (
        "Total disposition (GJ) without exports"
    )
    assert replace_units("coal percent (resource mix)", units) == (
        "coal percent (resource mix)"
    )
    assert replace_units("nameplate capacity (MW)", units) == (
        "nameplate capacity (MW)"
    )


@pytest.mark.parametrize("sheet", EGRID_SHEETS)
@pytest.mark.parametrize("metric", [False, True])
def test_each_category_is_one_block_of_columns(sheet, metric):
    categories = [category for _, _, category in get_egrid_columns(sheet, 2023, metric)]
    assert set(categories) - {None} <= set(EGRID_CATEGORIES)
    blocks = [
        category
        for i, category in enumerate(categories)
        if category is not None and (i == 0 or categories[i - 1] != category)
    ]
    assert len(blocks) == len(set(blocks))


def test_categories_by_sheet():
    def category_counts(sheet):
        counts = {}
        for _, _, category in get_egrid_columns(sheet, 2023):
            counts[category] = counts.get(category, 0) + 1
        return counts

    assert category_counts("UNT") == {None: 19, 2: 14}
    assert category_counts("GEN") == {None: 15, 1: 2}
    assert category_counts("PLNT") == {
        None: 37,
        1: 15,
        2: 22,
        3: 16,
        4: 8,
        5: 8,
        6: 8,
        7: 11,
        8: 4,
        9: 3,
        10: 11,
        11: 4,
        12: 3,
    }
    assert category_counts("BA") == {
        None: 4,
        1: 15,
        4: 8,
        5: 8,
        6: 8,
        7: 11,
        8: 4,
        9: 3,
        10: 11,
        11: 4,
        12: 3,
        13: 30,
        14: 30,
        15: 8,
        16: 11,
        17: 11,
    }
    assert category_counts("GGL") == {None: 6}


def test_gj_columns_have_the_category_of_their_column():
    metric = get_egrid_columns("PLNT", 2023, metric=True)
    categories = {code: category for code, _, category in metric}
    assert categories["PLNGENAN"] == categories["PLNGENAN2"] == 1
    assert categories["PLGENACL"] == categories["PLGENACL2"] == 7
