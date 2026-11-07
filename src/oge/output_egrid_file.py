"""Writes OGE annual results in the layout of EPA's eGRID data file.

The eGRID data file is an Excel workbook with one sheet for each level of aggregation
(unit, generator, plant, state, balancing authority, eGRID subregion, NERC region, U.S.
and grid gross loss). The first row of each sheet describes each column and the second
row holds the column codes. EPA publishes the file in U.S. units and in metric units.

This module defines the columns of each sheet in eGRID's order, using the headers of
the eGRID2023 U.S. units file, and the category of each column, which eGRID uses to
color the headers. The headers and category names fix eGRID's typos and
inconsistencies, such as "adjustement", "non-baseload" in some headers and
"nonbaseload" in others, and a fuel-type input emission rate category described as
"per MWh". Columns that OGE does not calculate are left blank.
"""

import re

# sheets of the eGRID data file, in workbook order. Each sheet name ends with the
# two-digit data year, e.g. "PLNT23"
EGRID_SHEETS = ["UNT", "GEN", "PLNT", "ST", "BA", "SRL", "NRL", "US", "GGL"]

# categories of columns, from the contents sheet of the eGRID data file
EGRID_CATEGORIES = {
    1: "Annual Values (generation, emissions, and heat input)",
    2: "Unadjusted Annual Values (emissions and heat input)",
    3: "Adjustment Values (emissions, heat input, heat rate)",
    4: "Output Emission Rates (emissions per MWh)",
    5: "Input Emission Rates (emissions per MMBtu)",
    6: "Combustion Output Rates (emissions per MWh)",
    7: "Generation by Fuel Type (MWh)",
    8: "Renewable and Non-Renewable Generation (MWh)",
    9: "Combustion and Non-Combustion Generation (MWh)",
    10: "Resource Mix (percentages)",
    11: "Renewable and Non-Renewable Resource Mix (percentages)",
    12: "Combustion and Non-Combustion Resource Mix (percentages)",
    13: "Output Emission Rates by Fuel Type (emissions by fuel type per MWh)",
    14: "Input Emission Rates by Fuel Type (emissions by fuel type per MMBtu)",
    15: "Nonbaseload Output Emission Rates (emissions per MWh)",
    16: "Nonbaseload Generation by Fuel Type (MWh)",
    17: "Nonbaseload Resource Mix (percentages)",
}

# units that eGRID's metric file converts, from U.S. units to metric units
METRIC_UNITS = {
    "tons": "metric tons",
    "lbs": "kg",
    "MMBtu": "GJ",
    "lb/MWh": "kg/MWh",
    "lb/MMBtu": "kg/GJ",
    "feet": "meters",
}
# eGRID's metric file also reports each column in these U.S. units in GJ, in a new
# column after it whose code ends in "2"
GJ_UNITS = {"MWh": "GJ", "lb/MWh": "kg/GJ"}


def get_egrid_columns(
    sheet: str, year: int, metric: bool = False
) -> list[tuple[str, str, int | None]]:
    """Returns the code, header and category of each column of a sheet of the eGRID
    data file.

    Args:
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        year (int): data year, used in the sequence number column codes.
        metric (bool, optional): whether to return the columns of the metric units
            file. Defaults to False, for the U.S. units file.

    Returns:
        list[tuple[str, str, int | None]]: (code, header, category) of each column, in
            eGRID's order. The category is a key of `EGRID_CATEGORIES`, or None for
            columns that are not in a category.
    """
    if sheet in REGION_SHEETS:
        prefix, label, id_columns = REGION_SHEETS[sheet]
        blocks = [(None, id_columns)] + [
            (category, [(prefix + code, f"{label} {header}") for code, header in block])
            for category, block in REGION_DATA_COLUMNS
        ]
    else:
        blocks = {
            "UNT": UNIT_COLUMNS,
            "GEN": GENERATOR_COLUMNS,
            "PLNT": PLANT_COLUMNS,
            "GGL": GRID_LOSS_COLUMNS,
        }[sheet]
    columns = [
        (code.format(yy=f"{year % 100:02d}"), header, category)
        for category, block in blocks
        for code, header in block
    ]
    if not metric:
        return columns
    metric_columns = []
    for code, header, category in columns:
        metric_columns.append((code, replace_units(header, METRIC_UNITS), category))
        if replace_units(header, GJ_UNITS) != header:
            metric_columns.append(
                (f"{code}2", replace_units(header, GJ_UNITS), category)
            )
    return metric_columns


def replace_units(header: str, units: dict[str, str]) -> str:
    """Replaces the units in parentheses in a column header, e.g. "(tons)".

    Args:
        header (str): column header.
        units (dict[str, str]): new unit for each unit to replace.

    Returns:
        str: the header with its units replaced. Headers with none of the units in
            `units` are returned unchanged.
    """
    return re.sub(
        r"\(([^()]+)\)",
        lambda match: f"({units.get(match.group(1), match.group(1))})",
        header,
    )


# Columns of each sheet, as blocks of (category, [(code, header), ...]) in eGRID's
# order. Columns that are not in a category have the category None. "{yy}" in a code
# is replaced by the two-digit data year
UNIT_COLUMNS = [
    # identifying and descriptive columns
    (
        None,
        [
            ("SEQUNT{yy}", "Unit file sequence number"),
            ("YEAR", "Data year"),
            ("PSTATABB", "Plant state abbreviation"),
            ("PNAME", "Plant name"),
            ("ORISPL", "DOE/EIA ORIS plant or facility code"),
            ("UNITID", "Unit ID"),
            ("PRMVR", "Prime mover"),
            ("UNTOPST", "Unit operational status"),
            ("CAPDFLAG", "CAPD program flag"),
            ("PRGCODE", "Program code(s)"),
            ("BOTFIRTY", "Unit bottom and firing type"),
            ("NUMGEN", "Number of associated generators"),
            ("FUELU1", "Unit primary fuel"),
            ("HRSOP", "Unit operating hours"),
        ],
    ),
    # 2) Unadjusted Annual Values (emissions and heat input)
    (
        2,
        [
            ("HTIAN", "Unit unadjusted annual heat input (MMBtu)"),
            ("HTIOZ", "Unit unadjusted ozone season heat input (MMBtu)"),
            ("NOXAN", "Unit unadjusted annual NOx emissions (tons)"),
            ("NOXOZ", "Unit unadjusted ozone season NOx emissions (tons)"),
            ("SO2AN", "Unit unadjusted annual SO2 emissions (tons)"),
            ("CO2AN", "Unit unadjusted annual CO2 emissions (tons)"),
            ("HGAN", "Unit unadjusted annual Hg emissions (lbs)"),
            ("HTIANSRC", "Unit unadjusted annual heat input source"),
            ("HTIOZSRC", "Unit unadjusted ozone season heat input source"),
            ("NOXANSRC", "Unit unadjusted annual NOx emissions source"),
            ("NOXOZSRC", "Unit unadjusted ozone season NOx emissions source"),
            ("SO2SRC", "Unit unadjusted annual SO2 emissions source"),
            ("CO2SRC", "Unit unadjusted annual CO2 emissions source"),
            ("HGSRC", "Unit unadjusted annual Hg emissions source"),
        ],
    ),
    # not in a category
    (
        None,
        [
            ("SO2CTLDV", "Unit SO2 (scrubber) first control device"),
            ("NOXCTLDV", "Unit NOx first control device"),
            ("HGCTLDV", "Unit Hg activated carbon injection system flag"),
            ("UNTYRONL", "Unit year on-line"),
            ("STACKHT", "Stack height (feet)"),
        ],
    ),
]

GENERATOR_COLUMNS = [
    # identifying and descriptive columns
    (
        None,
        [
            ("SEQGEN{yy}", "Generator file sequence number"),
            ("YEAR", "Data year"),
            ("PSTATABB", "Plant state abbreviation"),
            ("PNAME", "Plant name"),
            ("ORISPL", "DOE/EIA ORIS plant or facility code"),
            ("GENID", "Generator ID"),
            ("NUMBLR", "Number of associated boilers"),
            ("GENSTAT", "Generator status"),
            ("PRMVR", "Generator prime mover type"),
            ("FUELG1", "Generator primary fuel"),
            ("NAMEPCAP", "Generator nameplate capacity (MW)"),
            ("CFACT", "Generator capacity factor"),
        ],
    ),
    # 1) Annual Values (generation, emissions, and heat input)
    (
        1,
        [
            ("GENNTAN", "Generator annual net generation (MWh)"),
            ("GENNTOZ", "Generator ozone season net generation (MWh)"),
        ],
    ),
    # not in a category
    (
        None,
        [
            ("GENERSRC", "Generation data source"),
            ("GENYRONL", "Generator year on-line"),
            ("GENYRRET", "Generator planned or actual retirement year"),
        ],
    ),
]

PLANT_COLUMNS = [
    # identifying and descriptive columns
    (
        None,
        [
            ("SEQPLT{yy}", "Plant file sequence number"),
            ("YEAR", "Data year"),
            ("PSTATABB", "Plant state abbreviation"),
            ("PNAME", "Plant name"),
            ("ORISPL", "DOE/EIA ORIS plant or facility code"),
            ("OPRNAME", "Plant transmission or distribution system owner name"),
            ("OPRCODE", "Plant transmission or distribution system owner ID"),
            ("UTLSRVNM", "Utility name"),
            ("UTLSRVID", "Utility ID"),
            ("SECTOR", "Plant-level sector"),
            ("BANAME", "Balancing authority name"),
            ("BACODE", "Balancing authority code"),
            ("NERC", "NERC region acronym"),
            ("SUBRGN", "eGRID subregion acronym"),
            ("SRNAME", "eGRID subregion name"),
            ("ISORTO", "Plant associated ISO/RTO territory"),
            ("FIPSST", "Plant FIPS state code"),
            ("FIPSCNTY", "Plant FIPS county code"),
            ("CNTYNAME", "Plant county name"),
            ("LAT", "Plant latitude"),
            ("LON", "Plant longitude"),
            ("CAPDFLAG", "CAPD program flag"),
            ("NUMUNT", "Number of units"),
            ("NUMGEN", "Number of generators"),
            ("PLPRMFL", "Plant primary fuel"),
            ("PLFUELCT", "Plant primary fuel category"),
            (
                "COALFLAG",
                "Flag indicating if the plant burned or generated any amount of coal",
            ),
            ("CAPFAC", "Plant capacity factor"),
            ("NAMEPCAP", "Plant nameplate capacity (MW)"),
            ("NBFACTOR", "Nonbaseload factor"),
            ("RMBMFLAG", "Biogas/ biomass plant adjustment flag"),
            ("CHPFLAG", "Combined heat and power (CHP) plant adjustment flag"),
            ("USETHRMO", "CHP plant useful thermal output (MMBtu)"),
            ("PWRTOHT", "CHP plant power to heat ratio"),
            ("ELCALLOC", "CHP plant electric allocation factor"),
            ("PSFLAG", "Plant pumped storage flag"),
        ],
    ),
    # 1) Annual Values (generation, emissions, and heat input)
    (
        1,
        [
            ("PLHTIAN", "Plant annual heat input from combustion (MMBtu)"),
            ("PLHTIOZ", "Plant ozone season heat input from combustion (MMBtu)"),
            ("PLHTIANT", "Plant total annual heat input (MMBtu)"),
            ("PLHTIOZT", "Plant total ozone season heat input (MMBtu)"),
            ("PLNGENAN", "Plant annual net generation (MWh)"),
            ("PLNGENOZ", "Plant ozone season net generation (MWh)"),
            ("PLNGENNB", "Plant annual nonbaseload generation (MWh)"),
            ("PLNOXAN", "Plant annual NOx emissions (tons)"),
            ("PLNOXOZ", "Plant ozone season NOx emissions (tons)"),
            ("PLSO2AN", "Plant annual SO2 emissions (tons)"),
            ("PLCO2AN", "Plant annual CO2 emissions (tons)"),
            ("PLCH4AN", "Plant annual CH4 emissions (lbs)"),
            ("PLN2OAN", "Plant annual N2O emissions (lbs)"),
            ("PLCO2EQA", "Plant annual CO2 equivalent emissions (tons)"),
            ("PLHGAN", "Plant annual Hg emissions (lbs)"),
        ],
    ),
    # 4) Output Emission Rates (emissions per MWh)
    (
        4,
        [
            ("PLNOXRTA", "Plant annual NOx total output emission rate (lb/MWh)"),
            ("PLNOXRTO", "Plant ozone season NOx total output emission rate (lb/MWh)"),
            ("PLSO2RTA", "Plant annual SO2 total output emission rate (lb/MWh)"),
            ("PLCO2RTA", "Plant annual CO2 total output emission rate (lb/MWh)"),
            ("PLCH4RTA", "Plant annual CH4 total output emission rate (lb/MWh)"),
            ("PLN2ORTA", "Plant annual N2O total output emission rate (lb/MWh)"),
            (
                "PLC2ERTA",
                "Plant annual CO2 equivalent total output emission rate (lb/MWh)",
            ),
            ("PLHGRTA", "Plant annual Hg total output emission rate (lb/MWh)"),
        ],
    ),
    # 5) Input Emission Rates (emissions per MMBtu)
    (
        5,
        [
            ("PLNOXRA", "Plant annual NOx input emission rate (lb/MMBtu)"),
            ("PLNOXRO", "Plant ozone season NOx input emission rate (lb/MMBtu)"),
            ("PLSO2RA", "Plant annual SO2 input emission rate (lb/MMBtu)"),
            ("PLCO2RA", "Plant annual CO2 input emission rate (lb/MMBtu)"),
            ("PLCH4RA", "Plant annual CH4 input emission rate (lb/MMBtu)"),
            ("PLN2ORA", "Plant annual N2O input emission rate (lb/MMBtu)"),
            ("PLC2ERA", "Plant annual CO2 equivalent input emission rate (lb/MMBtu)"),
            ("PLHGRA", "Plant annual Hg input emission rate (lb/MMBtu)"),
        ],
    ),
    # 6) Combustion Output Rates (emissions per MWh)
    (
        6,
        [
            ("PLNOXCRT", "Plant annual NOx combustion output emission rate (lb/MWh)"),
            (
                "PLNOXCRO",
                "Plant ozone season NOx combustion output emission rate (lb/MWh)",
            ),
            ("PLSO2CRT", "Plant annual SO2 combustion output emission rate (lb/MWh)"),
            ("PLCO2CRT", "Plant annual CO2 combustion output emission rate (lb/MWh)"),
            ("PLCH4CRT", "Plant annual CH4 combustion output emission rate (lb/MWh)"),
            ("PLN2OCRT", "Plant annual N2O combustion output emission rate (lb/MWh)"),
            (
                "PLC2ECRT",
                "Plant annual CO2 equivalent combustion output emission rate (lb/MWh)",
            ),
            ("PLHGCRT", "Plant annual Hg combustion output emission rate (lb/MWh)"),
        ],
    ),
    # 2) Unadjusted Annual Values (emissions and heat input)
    (
        2,
        [
            ("UNNOX", "Plant unadjusted annual NOx emissions (tons)"),
            ("UNNOXOZ", "Plant unadjusted ozone season NOx emissions (tons)"),
            ("UNSO2", "Plant unadjusted annual SO2 emissions (tons)"),
            ("UNCO2", "Plant unadjusted annual CO2 emissions (tons)"),
            ("UNCH4", "Plant unadjusted annual CH4 emissions (lbs)"),
            ("UNN2O", "Plant unadjusted annual N2O emissions (lbs)"),
            ("UNCO2E", "Plant unadjusted annual CO2 equivalent emissions (tons)"),
            ("UNHG", "Plant unadjusted annual Hg emissions (lbs)"),
            ("UNHTI", "Plant unadjusted annual heat input from combustion (MMBtu)"),
            (
                "UNHTIOZ",
                "Plant unadjusted ozone season heat input from combustion (MMBtu)",
            ),
            ("UNHTIT", "Plant unadjusted total annual heat input (MMBtu)"),
            ("UNHTIOZT", "Plant unadjusted total ozone season heat input (MMBtu)"),
            ("UNNOXSRC", "Plant unadjusted annual NOx emissions source"),
            ("UNNOZSRC", "Plant unadjusted ozone season NOx emissions source"),
            ("UNSO2SRC", "Plant unadjusted annual SO2 emissions source"),
            ("UNCO2SRC", "Plant unadjusted annual CO2 emissions source"),
            ("UNCH4SRC", "Plant unadjusted annual CH4 emissions source"),
            ("UNN2OSRC", "Plant unadjusted annual N2O emissions source"),
            ("UNC2ESRC", "Plant unadjusted annual CO2 equivalent emissions source"),
            ("UNHGSRC", "Plant unadjusted annual Hg emissions source"),
            ("UNHTISRC", "Plant unadjusted annual heat input source"),
            ("UNHOZSRC", "Plant unadjusted ozone season heat input source"),
        ],
    ),
    # 3) Adjustment Values (emissions, heat input, heat rate)
    (
        3,
        [
            ("BIONOX", "Plant annual NOx biomass emissions (tons)"),
            ("BIONOXOZ", "Plant ozone season NOx biomass emissions (tons)"),
            ("BIOSO2", "Plant annual SO2 biomass emissions (tons)"),
            ("BIOCO2", "Plant annual CO2 biomass emissions (tons)"),
            ("BIOCH4", "Plant annual CH4 biomass emissions (lbs)"),
            ("BION2O", "Plant annual N2O biomass emissions (lbs)"),
            ("BIOCO2E", "Plant annual CO2 equivalent biomass emissions (tons)"),
            ("CHPCHTI", "Plant combustion heat input CHP adjustment value (MMBtu)"),
            (
                "CHPCHTIOZ",
                "Plant combustion annual ozone season heat input CHP adjustment value (MMBtu)",
            ),
            ("CHPNOX", "Plant annual NOx emissions CHP adjustment value (tons)"),
            (
                "CHPNOXOZ",
                "Plant ozone season NOx emissions CHP adjustment value (tons)",
            ),
            ("CHPSO2", "Plant annual SO2 emissions CHP adjustment value (tons)"),
            ("CHPCO2", "Plant annual CO2 emissions CHP adjustment value (tons)"),
            ("CHPCH4", "Plant annual CH4 emissions CHP adjustment value (lbs)"),
            ("CHPN2O", "Plant annual N2O emissions CHP adjustment value (lbs)"),
            (
                "CHPCO2E",
                "Plant annual CO2 equivalent emissions CHP adjustment value (tons)",
            ),
        ],
    ),
    # not in a category
    (
        None,
        [
            ("PLHTRT", "Plant nominal heat rate (Btu/kWh)"),
        ],
    ),
    # 7) Generation by Fuel Type (MWh)
    (
        7,
        [
            ("PLGENACL", "Plant annual coal net generation (MWh)"),
            ("PLGENAOL", "Plant annual oil net generation (MWh)"),
            ("PLGENAGS", "Plant annual gas net generation (MWh)"),
            ("PLGENANC", "Plant annual nuclear net generation (MWh)"),
            ("PLGENAHY", "Plant annual hydro net generation (MWh)"),
            ("PLGENABM", "Plant annual biomass net generation (MWh)"),
            ("PLGENAWI", "Plant annual wind net generation (MWh)"),
            ("PLGENASO", "Plant annual solar net generation (MWh)"),
            ("PLGENAGT", "Plant annual geothermal net generation (MWh)"),
            ("PLGENAOF", "Plant annual other fossil net generation (MWh)"),
            (
                "PLGENAOP",
                "Plant annual other unknown/purchased fuel net generation (MWh)",
            ),
        ],
    ),
    # 8) Renewable and Non-Renewable Generation (MWh)
    (
        8,
        [
            ("PLGENATN", "Plant annual total nonrenewables net generation (MWh)"),
            ("PLGENATR", "Plant annual total renewables net generation (MWh)"),
            (
                "PLGENATO",
                "Plant annual total nonrenewables other unknown/purchased net generation (MWh)",
            ),
            ("PLGENATH", "Plant annual total nonhydro renewables net generation (MWh)"),
        ],
    ),
    # 9) Combustion and Non-Combustion Generation (MWh)
    (
        9,
        [
            ("PLGENACY", "Plant annual total combustion net generation (MWh)"),
            ("PLGENACN", "Plant annual total noncombustion net generation (MWh)"),
            (
                "PLGENACO",
                "Plant annual total noncombustion other unknown/purchased net generation (MWh)",
            ),
        ],
    ),
    # 10) Resource Mix (percentages)
    (
        10,
        [
            ("PLCLPR", "Plant coal generation percent (resource mix)"),
            ("PLOLPR", "Plant oil generation percent (resource mix)"),
            ("PLGSPR", "Plant gas generation percent (resource mix)"),
            ("PLNCPR", "Plant nuclear generation percent (resource mix)"),
            ("PLHYPR", "Plant hydro generation percent (resource mix)"),
            ("PLBMPR", "Plant biomass generation percent (resource mix)"),
            ("PLWIPR", "Plant wind generation percent (resource mix)"),
            ("PLSOPR", "Plant solar generation percent (resource mix)"),
            ("PLGTPR", "Plant geothermal generation percent (resource mix)"),
            ("PLOFPR", "Plant other fossil generation percent (resource mix)"),
            (
                "PLOPPR",
                "Plant other unknown / purchased fuel generation percent (resource mix)",
            ),
        ],
    ),
    # 11) Renewable and Non-Renewable Resource Mix (percentages)
    (
        11,
        [
            ("PLTNPR", "Plant total nonrenewables generation percent (resource mix)"),
            ("PLTRPR", "Plant total renewables generation percent (resource mix)"),
            (
                "PLTOPR",
                "Plant total nonrenewables other unknown/purchased generation percent (resource mix)",
            ),
            (
                "PLTHPR",
                "Plant total nonhydro renewables generation percent (resource mix)",
            ),
        ],
    ),
    # 12) Combustion and Non-Combustion Resource Mix (percentages)
    (
        12,
        [
            ("PLCYPR", "Plant total combustion generation percent (resource mix)"),
            ("PLCNPR", "Plant total noncombustion generation percent (resource mix)"),
            (
                "PLCOPR",
                "Plant total noncombustion other unknown/purchased generation percent (resource mix)",
            ),
        ],
    ),
]

# The state, BA, eGRID subregion, NERC region and U.S. sheets have their own
# identifying columns followed by the same data columns. The data column codes start
# with the sheet's prefix and the headers with the sheet's label, e.g. "BANGENAN",
# "BA annual net generation (MWh)"
REGION_DATA_COLUMNS = [
    # not in a category
    (
        None,
        [
            ("NAMEPCAP", "nameplate capacity (MW)"),
        ],
    ),
    # 1) Annual Values (generation, emissions, and heat input)
    (
        1,
        [
            ("HTIAN", "annual heat input from combustion (MMBtu)"),
            ("HTIOZ", "ozone season heat input from combustion (MMBtu)"),
            ("HTIANT", "total annual heat input (MMBtu)"),
            ("HTIOZT", "total ozone season heat input (MMBtu)"),
            ("NGENAN", "annual net generation (MWh)"),
            ("NGENOZ", "ozone season net generation (MWh)"),
            ("NGENNB", "annual nonbaseload generation (MWh)"),
            ("NOXAN", "annual NOx emissions (tons)"),
            ("NOXOZ", "ozone season NOx emissions (tons)"),
            ("SO2AN", "annual SO2 emissions (tons)"),
            ("CO2AN", "annual CO2 emissions (tons)"),
            ("CH4AN", "annual CH4 emissions (lbs)"),
            ("N2OAN", "annual N2O emissions (lbs)"),
            ("CO2EQA", "annual CO2 equivalent emissions (tons)"),
            ("HGAN", "annual Hg emissions (lbs)"),
        ],
    ),
    # 4) Output Emission Rates (emissions per MWh)
    (
        4,
        [
            ("NOXRTA", "annual NOx total output emission rate (lb/MWh)"),
            ("NOXRTO", "ozone season NOx total output emission rate (lb/MWh)"),
            ("SO2RTA", "annual SO2 total output emission rate (lb/MWh)"),
            ("CO2RTA", "annual CO2 total output emission rate (lb/MWh)"),
            ("CH4RTA", "annual CH4 total output emission rate (lb/MWh)"),
            ("N2ORTA", "annual N2O total output emission rate (lb/MWh)"),
            ("C2ERTA", "annual CO2 equivalent total output emission rate (lb/MWh)"),
            ("HGRTA", "annual Hg total output emission rate (lb/MWh)"),
        ],
    ),
    # 5) Input Emission Rates (emissions per MMBtu)
    (
        5,
        [
            ("NOXRA", "annual NOx input emission rate (lb/MMBtu)"),
            ("NOXRO", "ozone season NOx input emission rate (lb/MMBtu)"),
            ("SO2RA", "annual SO2 input emission rate (lb/MMBtu)"),
            ("CO2RA", "annual CO2 input emission rate (lb/MMBtu)"),
            ("CH4RA", "annual CH4 input emission rate (lb/MMBtu)"),
            ("N2ORA", "annual N2O input emission rate (lb/MMBtu)"),
            ("C2ERA", "annual CO2 equivalent input emission rate (lb/MMBtu)"),
            ("HGRA", "annual Hg input emission rate (lb/MMBtu)"),
        ],
    ),
    # 6) Combustion Output Rates (emissions per MWh)
    (
        6,
        [
            ("NOXCRT", "annual NOx combustion output emission rate (lb/MWh)"),
            ("NOXCRO", "ozone season NOx combustion output emission rate (lb/MWh)"),
            ("SO2CRT", "annual SO2 combustion output emission rate (lb/MWh)"),
            ("CO2CRT", "annual CO2 combustion output emission rate (lb/MWh)"),
            ("CH4CRT", "annual CH4 combustion output emission rate (lb/MWh)"),
            ("N2OCRT", "annual N2O combustion output emission rate (lb/MWh)"),
            (
                "C2ECRT",
                "annual CO2 equivalent combustion output emission rate (lb/MWh)",
            ),
            ("HGCRT", "annual Hg combustion output emission rate (lb/MWh)"),
        ],
    ),
    # 13) Output Emission Rates by Fuel Type (emissions by fuel type per MWh)
    (
        13,
        [
            ("CNOXRT", "annual NOx coal output emission rate (lb/MWh)"),
            ("ONOXRT", "annual NOx oil output emission rate (lb/MWh)"),
            ("GNOXRT", "annual NOx gas output emission rate (lb/MWh)"),
            ("FSNXRT", "annual NOx fossil fuel output emission rate (lb/MWh)"),
            ("CNXORT", "ozone season NOx coal output emission rate (lb/MWh)"),
            ("ONXORT", "ozone season NOx oil output emission rate (lb/MWh)"),
            ("GNXORT", "ozone season NOx gas output emission rate (lb/MWh)"),
            ("FSNORT", "ozone season NOx fossil fuel output emission rate (lb/MWh)"),
            ("CSO2RT", "annual SO2 coal output emission rate (lb/MWh)"),
            ("OSO2RT", "annual SO2 oil output emission rate (lb/MWh)"),
            ("GSO2RT", "annual SO2 gas output emission rate (lb/MWh)"),
            ("FSS2RT", "annual SO2 fossil fuel output emission rate (lb/MWh)"),
            ("CCO2RT", "annual CO2 coal output emission rate (lb/MWh)"),
            ("OCO2RT", "annual CO2 oil output emission rate (lb/MWh)"),
            ("GCO2RT", "annual CO2 gas output emission rate (lb/MWh)"),
            ("FSC2RT", "annual CO2 fossil fuel output emission rate (lb/MWh)"),
            ("CCH4RT", "annual CH4 coal output emission rate (lb/MWh)"),
            ("OCH4RT", "annual CH4 oil output emission rate (lb/MWh)"),
            ("GCH4RT", "annual CH4 gas output emission rate (lb/MWh)"),
            ("FCH4RT", "annual CH4 fossil fuel output emission rate (lb/MWh)"),
            ("CN2ORT", "annual N2O coal output emission rate (lb/MWh)"),
            ("ON2ORT", "annual N2O oil output emission rate (lb/MWh)"),
            ("GN2ORT", "annual N2O gas output emission rate (lb/MWh)"),
            ("FN2ORT", "annual N2O fossil fuel output emission rate (lb/MWh)"),
            ("CC2ERT", "annual CO2 equivalent coal output emission rate (lb/MWh)"),
            ("OC2ERT", "annual CO2 equivalent oil output emission rate (lb/MWh)"),
            ("GC2ERT", "annual CO2 equivalent gas output emission rate (lb/MWh)"),
            (
                "FSC2ERT",
                "annual CO2 equivalent fossil fuel output emission rate (lb/MWh)",
            ),
            ("CHGRT", "annual Hg coal output emission rate (lb/MWh)"),
            ("FSHGRT", "annual Hg fossil fuel output emission rate (lb/MWh)"),
        ],
    ),
    # 14) Input Emission Rates by Fuel Type (emissions by fuel type per MMBtu)
    (
        14,
        [
            ("CNOXR", "annual NOx coal input emission rate (lb/MMBtu)"),
            ("ONOXR", "annual NOx oil input emission rate (lb/MMBtu)"),
            ("GNOXR", "annual NOx gas input emission rate (lb/MMBtu)"),
            ("FSNXR", "annual NOx fossil fuel input emission rate (lb/MMBtu)"),
            ("CNXOR", "ozone season NOx coal input emission rate (lb/MMBtu)"),
            ("ONXOR", "ozone season NOx oil input emission rate (lb/MMBtu)"),
            ("GNXOR", "ozone season NOx gas input emission rate (lb/MMBtu)"),
            ("FSNOR", "ozone season NOx fossil fuel input emission rate (lb/MMBtu)"),
            ("CSO2R", "annual SO2 coal input emission rate (lb/MMBtu)"),
            ("OSO2R", "annual SO2 oil input emission rate (lb/MMBtu)"),
            ("GSO2R", "annual SO2 gas input emission rate (lb/MMBtu)"),
            ("FSS2R", "annual SO2 fossil fuel input emission rate (lb/MMBtu)"),
            ("CCO2R", "annual CO2 coal input emission rate (lb/MMBtu)"),
            ("OCO2R", "annual CO2 oil input emission rate (lb/MMBtu)"),
            ("GCO2R", "annual CO2 gas input emission rate (lb/MMBtu)"),
            ("FSC2R", "annual CO2 fossil fuel input emission rate (lb/MMBtu)"),
            ("CCH4R", "annual CH4 coal input emission rate (lb/MMBtu)"),
            ("OCH4R", "annual CH4 oil input emission rate (lb/MMBtu)"),
            ("GCH4R", "annual CH4 gas input emission rate (lb/MMBtu)"),
            ("FCH4R", "annual CH4 fossil fuel input emission rate (lb/MMBtu)"),
            ("CN2OR", "annual N2O coal input emission rate (lb/MMBtu)"),
            ("ON2OR", "annual N2O oil input emission rate (lb/MMBtu)"),
            ("GN2OR", "annual N2O gas input emission rate (lb/MMBtu)"),
            ("FN2OR", "annual N2O fossil fuel input emission rate (lb/MMBtu)"),
            ("CC2ER", "annual CO2 equivalent coal input emission rate (lb/MMBtu)"),
            ("OC2ER", "annual CO2 equivalent oil input emission rate (lb/MMBtu)"),
            ("GC2ER", "annual CO2 equivalent gas input emission rate (lb/MMBtu)"),
            (
                "FSC2ER",
                "annual CO2 equivalent fossil fuel input emission rate (lb/MMBtu)",
            ),
            ("CHGR", "annual Hg coal input emission rate (lb/MMBtu)"),
            ("FSHGR", "annual Hg fossil fuel input emission rate (lb/MMBtu)"),
        ],
    ),
    # 15) Nonbaseload Output Emission Rates (emissions per MWh)
    (
        15,
        [
            ("NBNOX", "annual NOx nonbaseload output emission rate (lb/MWh)"),
            ("NBNXO", "ozone season NOx nonbaseload output emission rate (lb/MWh)"),
            ("NBSO2", "annual SO2 nonbaseload output emission rate (lb/MWh)"),
            ("NBCO2", "annual CO2 nonbaseload output emission rate (lb/MWh)"),
            ("NBCH4", "annual CH4 nonbaseload output emission rate (lb/MWh)"),
            ("NBN2O", "annual N2O nonbaseload output emission rate (lb/MWh)"),
            (
                "NBC2E",
                "annual CO2 equivalent nonbaseload output emission rate (lb/MWh)",
            ),
            ("NBHG", "annual Hg nonbaseload output emission rate (lb/MWh)"),
        ],
    ),
    # 7) Generation by Fuel Type (MWh)
    (
        7,
        [
            ("GENACL", "annual coal net generation (MWh)"),
            ("GENAOL", "annual oil net generation (MWh)"),
            ("GENAGS", "annual gas net generation (MWh)"),
            ("GENANC", "annual nuclear net generation (MWh)"),
            ("GENAHY", "annual hydro net generation (MWh)"),
            ("GENABM", "annual biomass net generation (MWh)"),
            ("GENAWI", "annual wind net generation (MWh)"),
            ("GENASO", "annual solar net generation (MWh)"),
            ("GENAGT", "annual geothermal net generation (MWh)"),
            ("GENAOF", "annual other fossil net generation (MWh)"),
            ("GENAOP", "annual other unknown/purchased fuel net generation (MWh)"),
        ],
    ),
    # 8) Renewable and Non-Renewable Generation (MWh)
    (
        8,
        [
            ("GENATN", "annual total nonrenewables net generation (MWh)"),
            ("GENATR", "annual total renewables net generation (MWh)"),
            (
                "GENATO",
                "annual total nonrenewables other unknown/purchased net generation (MWh)",
            ),
            ("GENATH", "annual total nonhydro renewables net generation (MWh)"),
        ],
    ),
    # 9) Combustion and Non-Combustion Generation (MWh)
    (
        9,
        [
            ("GENACY", "annual total combustion net generation (MWh)"),
            ("GENACN", "annual total noncombustion net generation (MWh)"),
            (
                "GENACO",
                "annual total noncombustion other unknown/purchased net generation (MWh)",
            ),
        ],
    ),
    # 10) Resource Mix (percentages)
    (
        10,
        [
            ("CLPR", "coal generation percent (resource mix)"),
            ("OLPR", "oil generation percent (resource mix)"),
            ("GSPR", "gas generation percent (resource mix)"),
            ("NCPR", "nuclear generation percent (resource mix)"),
            ("HYPR", "hydro generation percent (resource mix)"),
            ("BMPR", "biomass generation percent (resource mix)"),
            ("WIPR", "wind generation percent (resource mix)"),
            ("SOPR", "solar generation percent (resource mix)"),
            ("GTPR", "geothermal generation percent (resource mix)"),
            ("OFPR", "other fossil generation percent (resource mix)"),
            ("OPPR", "other unknown/purchased fuel generation percent (resource mix)"),
        ],
    ),
    # 11) Renewable and Non-Renewable Resource Mix (percentages)
    (
        11,
        [
            ("TNPR", "total nonrenewables generation percent (resource mix)"),
            ("TRPR", "total renewables generation percent (resource mix)"),
            (
                "TOPR",
                "total nonrenewables other unknown/purchased generation percent (resource mix)",
            ),
            ("THPR", "total nonhydro renewables generation percent (resource mix)"),
        ],
    ),
    # 12) Combustion and Non-Combustion Resource Mix (percentages)
    (
        12,
        [
            ("CYPR", "total combustion generation percent (resource mix)"),
            ("CNPR", "total noncombustion generation percent (resource mix)"),
            (
                "COPR",
                "total noncombustion other unknown/purchased generation percent (resource mix)",
            ),
        ],
    ),
    # 16) Nonbaseload Generation by Fuel Type (MWh)
    (
        16,
        [
            ("NBGNCL", "annual nonbaseload coal net generation (MWh)"),
            ("NBGNOL", "annual nonbaseload oil net generation (MWh)"),
            ("NBGNGS", "annual nonbaseload gas net generation (MWh)"),
            ("NBGNNC", "annual nonbaseload nuclear net generation (MWh)"),
            ("NBGNHY", "annual nonbaseload hydro net generation (MWh)"),
            ("NBGNBM", "annual nonbaseload biomass net generation (MWh)"),
            ("NBGNWI", "annual nonbaseload wind net generation (MWh)"),
            ("NBGNSO", "annual nonbaseload solar net generation (MWh)"),
            ("NBGNGT", "annual nonbaseload geothermal net generation (MWh)"),
            ("NBGNOF", "annual nonbaseload other fossil net generation (MWh)"),
            (
                "NBGNOP",
                "annual nonbaseload other unknown/purchased fuel net generation (MWh)",
            ),
        ],
    ),
    # 17) Nonbaseload Resource Mix (percentages)
    (
        17,
        [
            ("NBCLPR", "nonbaseload coal generation percent (resource mix)"),
            ("NBOLPR", "nonbaseload oil generation percent (resource mix)"),
            ("NBGSPR", "nonbaseload gas generation percent (resource mix)"),
            ("NBNCPR", "nonbaseload nuclear generation percent (resource mix)"),
            ("NBHYPR", "nonbaseload hydro generation percent (resource mix)"),
            ("NBBMPR", "nonbaseload biomass generation percent (resource mix)"),
            ("NBWIPR", "nonbaseload wind generation percent (resource mix)"),
            ("NBSOPR", "nonbaseload solar generation percent (resource mix)"),
            ("NBGTPR", "nonbaseload geothermal generation percent (resource mix)"),
            ("NBOFPR", "nonbaseload other fossil generation percent (resource mix)"),
            (
                "NBOPPR",
                "nonbaseload other unknown/purchased fuel generation percent (resource mix)",
            ),
        ],
    ),
]

# prefix, label and identifying columns of each regional sheet
REGION_SHEETS = {
    "ST": (
        "ST",
        "State",
        [
            ("YEAR", "Data year"),
            ("PSTATABB", "State abbreviation"),
            ("FIPSST", "FIPS state code"),
        ],
    ),
    "BA": (
        "BA",
        "BA",
        [
            ("YEAR", "Data year"),
            ("BANAME", "Balancing authority name"),
            ("BACODE", "Balancing authority code"),
        ],
    ),
    "SRL": (
        "SR",
        "eGRID subregion",
        [
            ("YEAR", "Data year"),
            ("SUBRGN", "eGRID subregion acronym"),
            ("SRNAME", "eGRID subregion name"),
        ],
    ),
    "NRL": (
        "NR",
        "NERC region",
        [
            ("YEAR", "Data year"),
            ("NERC", "NERC region acronym"),
            ("NERCNAME", "NERC region name"),
        ],
    ),
    "US": ("US", "U.S.", [("YEAR", "Data year")]),
}

GRID_LOSS_COLUMNS = [
    (
        None,
        [
            ("YEAR", "Data year"),
            (
                "REGION",
                "One of the three interconnect power grids in the U.S. (plus Alaska, Hawaii, and the entire U.S.)",
            ),
            ("ESTLOSS", "Estimated losses (MWh)"),
            ("TOTDISP", "Total disposition without exports (MWh)"),
            ("DIRCTUSE", "Direct use (MWh)"),
            (
                "GGRSLOSS",
                "Grid gross loss [Estimated losses/(Total disposition without exports - Direct use)]",
            ),
        ],
    ),
]
