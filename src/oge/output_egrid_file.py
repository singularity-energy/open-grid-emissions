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
"per MWh".

The workbooks are built from OGE's annual results, using OGE's plants, balancing
authorities and fuel categories. Columns that OGE does not calculate yet are left blank
and their headers are grayed out.
"""

import os
import re
from datetime import date
from importlib.metadata import version

import numpy as np
import pandas as pd
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell

from oge.constants import ConversionFactors
from oge.filepaths import outputs_folder, results_folder
from oge.load_data import load_ba_reference, load_campd_facility, load_pudl_table
from oge.logging_util import get_logger

logger = get_logger(__name__)

# sheets of the eGRID data file, in workbook order. Each sheet name ends with the
# two-digit data year, e.g. "PLNT23"
EGRID_SHEETS = ["UNT", "GEN", "PLNT", "ST", "BA", "SRL", "NRL", "US", "GGL"]

# description of each sheet, from the contents sheet of the eGRID data file
EGRID_SHEET_DESCRIPTIONS = {
    "UNT": "Unit",
    "GEN": "Generator",
    "PLNT": "Plant",
    "ST": "State",
    "BA": "Balancing authority area",
    "SRL": "eGRID subregion",
    "NRL": "NERC region",
    "US": "U.S.",
    "GGL": "Grid Gross Loss (%)",
}

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

# header fill color of each category, from the eGRID data file. Columns that are not
# in a category use the color of the None key
CATEGORY_COLORS = {
    None: "#F2F2F2",
    1: "#F2DCDB",
    2: "#E6B8B7",
    3: "#D58785",
    4: "#EBF1DE",
    5: "#C4D79B",
    6: "#76933C",
    7: "#DAEEF3",
    8: "#92CDDC",
    9: "#31869B",
    10: "#FDE9D9",
    11: "#FABF8F",
    12: "#E26B0A",
    13: "#FFFFCC",
    14: "#FFFF99",
    15: "#E4DFEC",
    16: "#B1A0C7",
    17: "#60497A",
}
# categories whose header fill is dark enough to need white text
DARK_CATEGORIES = [6, 9, 17]

# units that eGRID's metric file converts, as {U.S. unit: (metric unit, factor)}
METRIC_UNITS = {
    "tons": (
        "metric tons",
        ConversionFactors.short_ton_to_lbs * ConversionFactors.lb_to_kg / 1000,
    ),
    "lbs": ("kg", ConversionFactors.lb_to_kg),
    "MMBtu": ("GJ", ConversionFactors.mmbtu_to_GJ),
    "lb/MWh": ("kg/MWh", ConversionFactors.lb_to_kg),
    "lb/MMBtu": ("kg/GJ", ConversionFactors.lb_to_kg / ConversionFactors.mmbtu_to_GJ),
    "feet": ("meters", ConversionFactors.ft_to_m),
}
# eGRID's metric file also reports each column in these U.S. units in GJ, in a new
# column after it whose code ends in "2", as {U.S. unit: (GJ unit, factor)}
GJ_UNITS = {
    "MWh": ("GJ", ConversionFactors.mwh_to_GJ),
    "lb/MWh": ("kg/GJ", ConversionFactors.lb_to_kg / ConversionFactors.mwh_to_GJ),
}

# Layout of each sheet, from the eGRID data file: the number of columns at the left
# that stay visible when scrolling (None for no frozen panes) and the height of the
# header row
FROZEN_COLUMNS = {
    "UNT": 6,
    "GEN": 6,
    "PLNT": 5,
    "ST": 3,
    "BA": 3,
    "SRL": 3,
    "NRL": 3,
    "US": None,
    "GGL": None,
}
HEADER_HEIGHTS = {"UNT": 60.75, "GEN": 60.75, "GGL": 60.75}
DEFAULT_HEADER_HEIGHT = 67.5

# Column widths, from the eGRID data file. Columns not listed here have the default
# width. Widths of the regional sheets' data columns are keyed by the code without the
# sheet's prefix
DEFAULT_COLUMN_WIDTH = 13
COLUMN_WIDTHS = {
    "UNT": {"PNAME": 35.45, "PRGCODE": 17.73},
    "GEN": {"PNAME": 35.45, "CFACT": 13.82, "GENERSRC": 30.18, "GENYRRET": 15.82},
    "PLNT": {
        "PNAME": 35.45,
        "OPRNAME": 32.82,
        "UTLSRVNM": 35.73,
        "SECTOR": 15.54,
        "BANAME": 32.73,
        "SRNAME": 18,
        "CAPFAC": 13.82,
        "NAMEPCAP": 14,
        "PWRTOHT": 14,
        "UNCO2E": 15.73,
        "CHPCHTIOZ": 15.18,
        "BION2O": 14.73,
        "PLGENATO": 17.27,
        "PLGENACO": 17.27,
        "PLTOPR": 17,
        "PLCOPR": 17,
    },
    "ST": {"PSTATABB": 14.73, "FIPSST": 19.18},
    "BA": {"BANAME": 73.82, "BACODE": 19.18},
    "SRL": {"SUBRGN": 14.73, "SRNAME": 19.18},
    "NRL": {"NERC": 14.73, "NERCNAME": 30.18},
    "US": {"YEAR": 14.82},
    "GGL": {"YEAR": 10.73, "REGION": 22, "ESTLOSS": 11.82, "GGRSLOSS": 23.73},
}
REGION_DATA_WIDTHS = {
    "NAMEPCAP": 14.73,
    "GENATO": 17.27,
    "GENATH": 14.73,
    "GENACO": 17.27,
    "CLPR": 14.73,
    "TOPR": 17.27,
    "THPR": 14.73,
    "COPR": 17.27,
    "NBGNOL": 14.73,
    "NBGNBM": 15.73,
    "NBGNWI": 14.73,
    "NBGNOP": 16.54,
    "NBCLPR": 14.73,
    "NBBMPR": 15.73,
    "NBWIPR": 14.73,
    "NBOPPR": 16.27,
}

# Number formats, from the eGRID data file. Quantities are whole numbers with negative
# values in parentheses, emission rates have three decimals, resource mixes are
# percentages, and the exceptions below are keyed by the code on the unit, generator
# and plant sheets
QUANTITY_FORMAT = "#,##0;\\ \\(#,##0\\)"
RATE_FORMAT = "#,##0.000"
PERCENT_FORMAT = "0.0%"
NUMBER_FORMATS = {
    "NAMEPCAP": "#,##0.0",
    "HRSOP": "#,##0.0",
    "NUMUNT": "#,##0",
    "NUMGEN": "#,##0",
    "PWRTOHT": "#,##0",
    "CAPFAC": "#,##0.0000",
    "NBFACTOR": "#,##0.0000",
    "ELCALLOC": "#,##0.0000",
    "CFACT": "#,##0.000;\\ \\(#,##0.000\\)",
}
# the unit sheet reports emissions with one decimal, and the unit sheet's heat input,
# the plant sheet's unadjusted values and the grid gross loss sheet are whole numbers
# without parentheses
UNIT_MASS_FORMAT = "#,##0.0"
# plant Hg masses and the CH4 and N2O emission rates by fuel type are small, so they
# have more decimals
HG_MASS_FORMAT = "#,##0.000"
SMALL_RATE_FORMAT = "#,##0.0000"

# links on the contents sheet
OGE_REPOSITORY_URL = "https://github.com/singularity-energy/open-grid-emissions"
OGE_DOCUMENTATION_URL = "https://docs.singularity.energy/docs/open-grid-emissions"
OGE_ISSUES_URL = "https://github.com/singularity-energy/open-grid-emissions/issues"

# months of the ozone season, May to September
OZONE_SEASON_MONTHS = [5, 6, 7, 8, 9]
# OGE fuel categories in each of eGRID's fuel groups. The groups use OGE's categories,
# which are defined differently from eGRID's (see the eGRID comparison doc): OGE's
# natural gas includes other gases such as blast furnace gas, its waste fills eGRID's
# other fossil group, and storage is counted with other unknown/purchased fuel
EGRID_FUEL_GROUPS = {
    "coal": ["coal"],
    "oil": ["petroleum"],
    "gas": ["natural_gas"],
    "nuclear": ["nuclear"],
    "hydro": ["hydro"],
    "biomass": ["biomass"],
    "wind": ["wind"],
    "solar": ["solar"],
    "geothermal": ["geothermal"],
    "other_fossil": ["waste"],
    "other": ["other", "storage"],
}
EGRID_FUEL_GROUPS |= {
    "nonrenewable": [
        category
        for group in ["coal", "oil", "gas", "other_fossil", "nuclear"]
        for category in EGRID_FUEL_GROUPS[group]
    ],
    "renewable": [
        category
        for group in ["biomass", "hydro", "wind", "solar", "geothermal"]
        for category in EGRID_FUEL_GROUPS[group]
    ],
    "nonhydro_renewable": [
        category
        for group in ["biomass", "wind", "solar", "geothermal"]
        for category in EGRID_FUEL_GROUPS[group]
    ],
    "combustion": [
        category
        for group in ["coal", "oil", "gas", "other_fossil", "biomass"]
        for category in EGRID_FUEL_GROUPS[group]
    ],
    "noncombustion": [
        category
        for group in ["nuclear", "hydro", "wind", "solar", "geothermal"]
        for category in EGRID_FUEL_GROUPS[group]
    ],
    "fossil": [
        category
        for group in ["coal", "oil", "gas", "other_fossil"]
        for category in EGRID_FUEL_GROUPS[group]
    ],
}
# eGRID fuel group of each generation by fuel column, keyed by the code without the
# sheet's prefix. eGRID's other unknown/purchased fuel generation is counted as
# nonrenewable and noncombustion, so it is in three columns
GENERATION_BY_FUEL_CODES = {
    "GENACL": "coal",  # coal net generation
    "GENAOL": "oil",  # oil net generation
    "GENAGS": "gas",  # gas net generation
    "GENANC": "nuclear",  # nuclear net generation
    "GENAHY": "hydro",  # hydro net generation
    "GENABM": "biomass",  # biomass net generation
    "GENAWI": "wind",  # wind net generation
    "GENASO": "solar",  # solar net generation
    "GENAGT": "geothermal",  # geothermal net generation
    "GENAOF": "other_fossil",  # other fossil net generation
    "GENAOP": "other",  # other unknown/purchased fuel net generation
    "GENATN": "nonrenewable",  # total nonrenewables net generation
    "GENATR": "renewable",  # total renewables net generation
    "GENATO": "other",  # other unknown/purchased, with the nonrenewables
    "GENATH": "nonhydro_renewable",  # total nonhydro renewables net generation
    "GENACY": "combustion",  # total combustion net generation
    "GENACN": "noncombustion",  # total noncombustion net generation
    "GENACO": "other",  # other unknown/purchased, with the noncombustion
}

# Emission rates are calculated from the columns of each sheet. The mass column of
# each pollutant's rates, keyed by the pollutant's code in the rate columns, e.g. the
# CO2 equivalent output rate "PLC2ERTA" divides "PLCO2EQA" by "PLNGENAN"
EMISSION_RATE_MASSES = {
    "NOX": "NOXAN",  # annual NOx emissions
    "SO2": "SO2AN",  # annual SO2 emissions
    "CO2": "CO2AN",  # annual CO2 emissions
    "CH4": "CH4AN",  # annual CH4 emissions
    "N2O": "N2OAN",  # annual N2O emissions
    "C2E": "CO2EQA",  # annual CO2 equivalent emissions
}
# code of the ozone season combustion net generation, which is not an eGRID column
OZONE_COMBUSTION_GENERATION = "GENACYOZ"
# the denominator of each type of emission rate: output rates divide by net
# generation, input rates by combustion heat input, and combustion output rates by
# combustion net generation. The NOx rates also have an ozone season version, which
# divides ozone season NOx by the ozone season denominator. Each entry is (annual
# denominator, ozone season denominator, code of the ozone season NOx rate)
EMISSION_RATE_DENOMINATORS = {
    # output rate: net generation, ozone season net generation
    "RTA": ("NGENAN", "NGENOZ", "RTO"),
    # input rate: combustion heat input, ozone season combustion heat input
    "RA": ("HTIAN", "HTIOZ", "RO"),
    # combustion output rate: combustion net generation, and its ozone season value
    "CRT": ("GENACY", OZONE_COMBUSTION_GENERATION, "CRO"),
}
# generation column of each resource mix column, keyed by the codes without the
# sheet's prefix. Shares are of the sum of the positive generation of the fuel columns
RESOURCE_MIX_GENERATION = {
    "CLPR": "GENACL",  # coal
    "OLPR": "GENAOL",  # oil
    "GSPR": "GENAGS",  # gas
    "NCPR": "GENANC",  # nuclear
    "HYPR": "GENAHY",  # hydro
    "BMPR": "GENABM",  # biomass
    "WIPR": "GENAWI",  # wind
    "SOPR": "GENASO",  # solar
    "GTPR": "GENAGT",  # geothermal
    "OFPR": "GENAOF",  # other fossil
    "OPPR": "GENAOP",  # other unknown/purchased fuel
    "TNPR": "GENATN",  # total nonrenewables
    "TRPR": "GENATR",  # total renewables
    "TOPR": "GENATO",  # other unknown/purchased, with the nonrenewables
    "THPR": "GENATH",  # total nonhydro renewables
    "CYPR": "GENACY",  # total combustion
    "CNPR": "GENACN",  # total noncombustion
    "COPR": "GENACO",  # other unknown/purchased, with the noncombustion
}
# Emission rates by fuel type divide the emissions of a fuel group's rows in the power
# sector results by their net generation (output rates) or heat input (input rates).
# OGE column of each pollutant's mass, keyed by the pollutant's code
FUEL_TYPE_RATE_MASSES = {
    "NOX": "nox_mass_lb_for_electricity_adjusted",
    "SO2": "so2_mass_lb_for_electricity_adjusted",
    "CO2": "co2_mass_lb_for_electricity_adjusted",
    "CH4": "ch4_mass_lb_for_electricity_adjusted",
    "N2O": "n2o_mass_lb_for_electricity_adjusted",
    "C2E": "co2e_mass_lb_for_electricity_adjusted",
}
# output rate code of each pollutant and fuel group, without the sheet's prefix. The
# input rate code is the same without the final "T", e.g. "CNOXR". "NOXOZ" is the
# ozone season NOx rate, which divides by ozone season generation or heat input
FUEL_TYPE_OUTPUT_RATE_CODES = {
    "NOX": {"coal": "CNOXRT", "oil": "ONOXRT", "gas": "GNOXRT", "fossil": "FSNXRT"},
    "NOXOZ": {"coal": "CNXORT", "oil": "ONXORT", "gas": "GNXORT", "fossil": "FSNORT"},
    "SO2": {"coal": "CSO2RT", "oil": "OSO2RT", "gas": "GSO2RT", "fossil": "FSS2RT"},
    "CO2": {"coal": "CCO2RT", "oil": "OCO2RT", "gas": "GCO2RT", "fossil": "FSC2RT"},
    "CH4": {"coal": "CCH4RT", "oil": "OCH4RT", "gas": "GCH4RT", "fossil": "FCH4RT"},
    "N2O": {"coal": "CN2ORT", "oil": "ON2ORT", "gas": "GN2ORT", "fossil": "FN2ORT"},
    "C2E": {"coal": "CC2ERT", "oil": "OC2ERT", "gas": "GC2ERT", "fossil": "FSC2ERT"},
}

# generation columns of the individual fuels, whose positive generation sums to the
# total of the resource mix
RESOURCE_MIX_FUELS = [
    "GENACL",  # coal
    "GENAOL",  # oil
    "GENAGS",  # gas
    "GENANC",  # nuclear
    "GENAHY",  # hydro
    "GENABM",  # biomass
    "GENAWI",  # wind
    "GENASO",  # solar
    "GENAGT",  # geothermal
    "GENAOF",  # other fossil
    "GENAOP",  # other unknown/purchased fuel
]

# OGE column for each column of the plant sheet that holds the same quantity. Plant
# attributes come from OGE's plant attributes, or from EIA-860 when OGE does not have
# them. OGE masses are in lb, and are converted to the units of each eGRID column
PLANT_DATA_SOURCES = {
    "PSTATABB": "state",
    "PNAME": "plant_name_eia",
    "ORISPL": "plant_id_eia",
    "OPRNAME": "transmission_distribution_owner_name",
    "OPRCODE": "transmission_distribution_owner_id",
    "UTLSRVNM": "utility_name_eia",
    "UTLSRVID": "utility_id_eia",
    "SECTOR": "sector",
    "BANAME": "ba_name",
    "BACODE": "ba_code",
    "NERC": "nerc_region",
    "FIPSST": "state_id_fips",
    "FIPSCNTY": "county_id_fips",
    "CNTYNAME": "county",
    "LAT": "latitude",
    "LON": "longitude",
    "PLPRMFL": "plant_primary_fuel",
    "NAMEPCAP": "capacity_mw",
    "PLHTIAN": "combustion_fuel_consumed_for_electricity_mmbtu",
    "PLHTIOZ": "combustion_fuel_consumed_for_electricity_mmbtu_ozone_season",
    "PLHTIANT": "fuel_consumed_for_electricity_mmbtu",
    "PLHTIOZT": "fuel_consumed_for_electricity_mmbtu_ozone_season",
    "PLNGENAN": "net_generation_mwh",
    "PLNGENOZ": "net_generation_mwh_ozone_season",
    "PLNOXOZ": "nox_mass_lb_for_electricity_adjusted_ozone_season",
    "PLCO2EQA": "co2e_mass_lb_for_electricity_adjusted",
    # not an eGRID column: used for the ozone season combustion output rate
    f"PL{OZONE_COMBUSTION_GENERATION}": "net_generation_mwh_combustion_ozone_season",
    "PLNOXAN": "nox_mass_lb_for_electricity_adjusted",
    "PLSO2AN": "so2_mass_lb_for_electricity_adjusted",
    "PLCO2AN": "co2_mass_lb_for_electricity_adjusted",
    "PLCH4AN": "ch4_mass_lb_for_electricity_adjusted",
    "PLN2OAN": "n2o_mass_lb_for_electricity_adjusted",
    "PLHGAN": "hg_mass_lb_for_electricity",
    "UNNOX": "nox_mass_lb",
    "UNNOXOZ": "nox_mass_lb_ozone_season",
    "UNSO2": "so2_mass_lb",
    "UNCO2": "co2_mass_lb",
    "UNCH4": "ch4_mass_lb",
    "UNN2O": "n2o_mass_lb",
    "UNHG": "hg_mass_lb",
    "UNCO2E": "co2e_mass_lb",
    "UNHTI": "combustion_fuel_consumed_mmbtu",
    "UNHTIOZ": "combustion_fuel_consumed_mmbtu_ozone_season",
    "UNHTIT": "fuel_consumed_mmbtu",
    "UNHTIOZT": "fuel_consumed_mmbtu_ozone_season",
} | {
    f"PL{code}": f"net_generation_mwh_{group}"
    for code, group in GENERATION_BY_FUEL_CODES.items()
}
# OGE column for each column of the unit sheet that holds the same quantity, from the
# annual totals of OGE's cleaned hourly CEMS data
UNIT_DATA_SOURCES = {
    "PSTATABB": "state",
    "PNAME": "plant_name_eia",
    "ORISPL": "plant_id_eia",
    "UNITID": "emissions_unit_id_epa",
    "FUELU1": "energy_source_code",
    "HRSOP": "operating_hours",
    "HTIAN": "fuel_consumed_mmbtu",
    "HTIOZ": "fuel_consumed_mmbtu_ozone_season",
    "NOXAN": "nox_mass_lb",
    "NOXOZ": "nox_mass_lb_ozone_season",
    "SO2AN": "so2_mass_lb",
    "CO2AN": "co2_mass_lb",
    "HGAN": "hg_mass_lb",
    # unit attributes from the CAMPD facility file, as eGRID codes
    "PRMVR": "prime_mover_code",
    "UNTOPST": "operating_status",
    "PRGCODE": "program_code_epa",
    "BOTFIRTY": "firing_type",
    "SO2CTLDV": "so2_controls",
    "NOXCTLDV": "nox_controls",
    "HGCTLDV": "hg_controls_flag",
    "UNTYRONL": "operating_year",
}
# OGE column for each column of the generator sheet that holds the same quantity, from
# the annual totals of OGE's EIA-923 data allocated to generators
GENERATOR_DATA_SOURCES = {
    "PSTATABB": "state",
    "PNAME": "plant_name_eia",
    "ORISPL": "plant_id_eia",
    "GENID": "generator_id",
    "PRMVR": "prime_mover_code",
    "GENNTAN": "net_generation_mwh",
    "GENNTOZ": "net_generation_mwh_ozone_season",
    "GENYRONL": "generator_operating_year",
    "GENYRRET": "generator_retirement_year",
    # generator attributes from EIA-860
    "GENSTAT": "operational_status_code",
    "FUELG1": "energy_source_code_1",
    "NAMEPCAP": "capacity_mw",
}

# eGRID codes for the unit attributes in the CAMPD facility file, from eGRID's
# crosswalks. Operating statuses are EIA status codes
CAMPD_OPERATING_STATUS_CODES = {
    "Operating": "OP",
    "Retired": "RE",
    "Long-term Cold Storage": "OS",
    "Future": "P",
}
CAMPD_UNIT_TYPE_CODES = {
    "Arch-fired boiler": "ARCH",
    "Bubbling fluidized bed boiler": "FLUIDIZED",
    "Cell burner boiler": "CELL",
    "Cement Kiln": "KILN",
    "Circulating fluidized bed boiler": "FLUIDIZED",
    "Combined cycle": "CC",
    "Combustion turbine": "CT",
    "Cyclone boiler": "CYCLONE",
    "Dry bottom turbo-fired boiler": "TURBO",
    "Dry bottom vertically-fired boiler": "VERTICAL",
    "Dry bottom wall-fired boiler": "WALL",
    "Integrated gasification combined cycle": "IGC",
    "Internal combustion engine": "ENGINE",
    "Other boiler": "OTHER BOILER",
    "Other turbine": "OTHER TURBINE",
    "Pressurized fluidized bed boiler": "FLUIDIZED",
    "Process Heater": "PROCESS",
    "Stoker": "STOKER",
    "Tangentially-fired": "TANGENTIAL",
    "Wet bottom turbo-fired boiler": "TURBO",
    "Wet bottom vertically-fired boiler": "VERTICAL",
    "Wet bottom wall-fired boiler": "WALL",
}
# EIA prime mover of each CAMPD unit type, as eGRID assigns it. One correction: eGRID's
# code groups arch-fired boilers with combustion turbines (GT), which looks like an
# error since they are boilers, so they are steam turbines (ST) here like other boilers
CAMPD_UNIT_TYPE_PRIME_MOVERS = {
    "Arch-fired boiler": "ST",
    "Bubbling fluidized bed boiler": "ST",
    "Cell burner boiler": "ST",
    "Cement Kiln": "ST",
    "Circulating fluidized bed boiler": "ST",
    "Combined cycle": "CT",
    "Combustion turbine": "GT",
    "Cyclone boiler": "ST",
    "Dry bottom turbo-fired boiler": "ST",
    "Dry bottom vertically-fired boiler": "ST",
    "Dry bottom wall-fired boiler": "ST",
    "Integrated gasification combined cycle": "ST",
    "Other boiler": "ST",
    "Other turbine": "OT",
    "Pressurized fluidized bed boiler": "ST",
    "Process Heater": "ST",
    "Stoker": "ST",
    "Tangentially-fired": "ST",
    "Wet bottom turbo-fired boiler": "ST",
    "Wet bottom vertically-fired boiler": "ST",
    "Wet bottom wall-fired boiler": "ST",
}
CAMPD_SO2_CONTROL_CODES = {
    "Dry Lime FGD": "DL",
    "Dry Sorbent Injection": "DSI",
    "Dual Alkali": "DA",
    "Fluidized Bed Limestone Injection": "FBL",
    "Other": "O",
    "Sodium Based": "SB",
    "Wet Lime FGD": "WL",
    "Wet Limestone": "WLS",
}
CAMPD_NOX_CONTROL_CODES = {
    "Ammonia Injection": "NH3",
    "Combustion Modification/Fuel Reburning": "CM",
    "Dry Low NOx Burners": "LNB",
    "Low NOx Burner Technology (Dry Bottom only)": "LNB",
    "Low NOx Burner Technology w/ Closed-coupled OFA": "LNC1",
    "Low NOx Burner Technology w/ Separated OFA": "LNC2",
    "Low NOx Burner Technology w/ Closed-coupled/Separated OFA": "LNC3",
    "Low NOx Burner Technology w/ Overfire Air": "LNBO",
    "Low NOx Cell Burner": "LNCB",
    "Other": "O",
    "Overfire Air": "OFA",
    "Selective Catalytic Reduction": "SCR",
    "Selective Non-catalytic Reduction": "SNCR",
    "Steam Injection": "STM",
    "Water Injection": "H2O",
}

# CAMPD notes when a status, unit type or control started or ended in the year, e.g.
# "Water Injection (Retired Jun 01, 2024)"
CAMPD_DATE_NOTE = r"\s*\((Began|Ended|Retired|Started)[^)]*\)"

# eGRID plant sector of each EIA sector ID
EIA_SECTORS = {
    1: "Electric Utility",
    2: "IPP Non-CHP",
    3: "IPP CHP",
    4: "Commercial Non-CHP",
    5: "Commercial CHP",
    6: "Industrial Non-CHP",
    7: "Industrial CHP",
}
# OGE fuel category or eGRID fuel group, column and season of each data column of the
# BA and U.S. sheets that holds the same quantity, keyed by the code without the
# sheet's prefix. A fuel group sums its fuel categories in EGRID_FUEL_GROUPS, and the
# ozone season sums May to September from the monthly results
REGION_DATA_SOURCES = {
    "HTIAN": ("combustion", "fuel_consumed_for_electricity_mmbtu", "annual"),
    "HTIOZ": ("combustion", "fuel_consumed_for_electricity_mmbtu", "ozone"),
    "HTIANT": ("total", "fuel_consumed_for_electricity_mmbtu", "annual"),
    "HTIOZT": ("total", "fuel_consumed_for_electricity_mmbtu", "ozone"),
    "NGENAN": ("total", "net_generation_mwh", "annual"),
    "NGENOZ": ("total", "net_generation_mwh", "ozone"),
    "NOXAN": ("total", "nox_mass_lb_for_electricity_adjusted", "annual"),
    "NOXOZ": ("total", "nox_mass_lb_for_electricity_adjusted", "ozone"),
    "SO2AN": ("total", "so2_mass_lb_for_electricity_adjusted", "annual"),
    "CO2AN": ("total", "co2_mass_lb_for_electricity_adjusted", "annual"),
    "CH4AN": ("total", "ch4_mass_lb_for_electricity_adjusted", "annual"),
    "N2OAN": ("total", "n2o_mass_lb_for_electricity_adjusted", "annual"),
    "CO2EQA": ("total", "co2e_mass_lb_for_electricity_adjusted", "annual"),
    # not an eGRID column: used for the ozone season combustion output rate
    OZONE_COMBUSTION_GENERATION: ("combustion", "net_generation_mwh", "ozone"),
} | {
    code: (group, "net_generation_mwh", "annual")
    for code, group in GENERATION_BY_FUEL_CODES.items()
}


def write_egrid_files(year: int, path_prefix: str, skip_outputs: bool = False):
    """Writes OGE's annual results in the layout of the eGRID data file, in U.S. units
    and in metric units, to the results/{path_prefix}eGRID/ folder.

    Args:
        year (int): data year.
        path_prefix (str): prefix of the results folder of the year, e.g. "2024/".
        skip_outputs (bool, optional): if True, the files are not written. Defaults to
            False.
    """
    if skip_outputs:
        return
    tables = build_egrid_tables(year, path_prefix)
    metric_tables = {
        sheet: convert_to_metric(table, sheet, year) for sheet, table in tables.items()
    }
    os.makedirs(results_folder(f"{path_prefix}eGRID"), exist_ok=True)
    write_egrid_workbook(
        tables,
        year,
        results_folder(f"{path_prefix}eGRID/oge_egrid{year}_data.xlsx"),
        metric=False,
    )
    write_egrid_workbook(
        metric_tables,
        year,
        results_folder(f"{path_prefix}eGRID/oge_egrid{year}_data_metric.xlsx"),
        metric=True,
    )


def build_egrid_tables(year: int, path_prefix: str) -> dict[str, pd.DataFrame]:
    """Builds each sheet of the eGRID data file in U.S. units from OGE's annual results
    and intermediate outputs.

    Sheets that OGE does not calculate yet have no rows.

    Args:
        year (int): data year.
        path_prefix (str): prefix of the results folder of the year, e.g. "2024/".

    Returns:
        dict[str, pd.DataFrame]: table of each sheet, keyed by its name in
            `EGRID_SHEETS`, with the columns of `get_egrid_columns()`.
    """
    annual_folder = f"{path_prefix}power_sector_data/annual/us_units"
    ba_codes = sorted(
        file.removesuffix(".csv")
        for file in os.listdir(results_folder(annual_folder))
        if file.endswith(".csv") and file != "US.csv"
    )
    tables = {
        "UNT": build_unit_table(year, path_prefix),
        "GEN": build_generator_table(year, path_prefix),
        "PLNT": build_plant_table(year, path_prefix),
        "BA": build_region_table(year, path_prefix, "BA", ba_codes),
        "US": build_region_table(year, path_prefix, "US", ["US"]),
    }
    return {
        sheet: tables.get(sheet, pd.DataFrame()).reindex(
            columns=[code for code, _, _ in get_egrid_columns(sheet, year)]
        )
        for sheet in EGRID_SHEETS
    }


def build_unit_table(year: int, path_prefix: str) -> pd.DataFrame:
    """Builds the unit sheet of the eGRID data file from OGE's cleaned hourly CEMS data
    and the unit attributes in the CAMPD facility file.

    The sheet has the units that report to CEMS. Hg is the sum of the hours with Hg
    data, as in eGRID, and is blank for units that do not report Hg.

    Args:
        year (int): data year.
        path_prefix (str): prefix of the outputs and results folders of the year, e.g.
            "2024/".

    Returns:
        pd.DataFrame: one row for each CEMS unit, in U.S. units, with the columns of
            `UNIT_DATA_SOURCES`, the sequence number, the year and the CAPD flag.
    """
    keys = ["plant_id_eia", "emissions_unit_id_epa"]
    data_columns = [
        "fuel_consumed_mmbtu",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb",
        "hg_mass_lb",
    ]
    cems = pd.read_csv(
        outputs_folder(f"{path_prefix}cems_cleaned_{year}.csv.zip"),
        usecols=keys
        + data_columns
        + ["energy_source_code", "operating_time_hours", "report_date"],
        dtype={"emissions_unit_id_epa": "str", "energy_source_code": "str"},
        parse_dates=["report_date"],
    )
    ozone_season = cems["report_date"].dt.month.isin(OZONE_SEASON_MONTHS)
    for column in ["fuel_consumed_mmbtu", "nox_mass_lb"]:
        cems[f"{column}_ozone_season"] = cems[column].where(ozone_season, 0)
        data_columns.append(f"{column}_ozone_season")
    units = (
        cems.groupby(keys)
        .agg(
            energy_source_code=("energy_source_code", "first"),
            # eGRID counts the hours in which the unit operated, as CAMPD does
            operating_hours=("operating_time_hours", lambda hours: (hours > 0).sum()),
            **{
                column: (column, lambda values: values.sum(min_count=1))
                for column in data_columns
            },
        )
        .reset_index()
    )
    units = units.merge(
        load_plant_names(path_prefix), on="plant_id_eia", how="left"
    ).merge(load_campd_unit_attributes(year), on=keys, how="left")
    table = build_sheet_rows(
        units, UNIT_DATA_SOURCES, year, "SEQUNT", ["state", "plant_name_eia", *keys]
    )
    # every unit in OGE's CEMS data reports to EPA's Clean Air Markets Division
    table["CAPDFLAG"] = "Yes"
    return convert_mass_to_tons(table, "UNT", year)


def build_generator_table(year: int, path_prefix: str) -> pd.DataFrame:
    """Builds the generator sheet of the eGRID data file from OGE's EIA-923 data
    allocated to generators and the generator attributes in EIA-860.

    Args:
        year (int): data year.
        path_prefix (str): prefix of the outputs and results folders of the year, e.g.
            "2024/".

    Returns:
        pd.DataFrame: one row for each generator in OGE's allocated EIA-923 data, with
            the columns of `GENERATOR_DATA_SOURCES`, the sequence number and the year.
    """
    keys = ["plant_id_eia", "generator_id"]
    allocated = pd.read_csv(
        outputs_folder(f"{path_prefix}eia923_allocated_{year}.csv.zip"),
        usecols=keys + ["prime_mover_code", "net_generation_mwh", "report_date"],
        dtype={"generator_id": "str"},
        parse_dates=["report_date"],
    )
    allocated["net_generation_mwh_ozone_season"] = allocated[
        "net_generation_mwh"
    ].where(allocated["report_date"].dt.month.isin(OZONE_SEASON_MONTHS), 0)
    generators = (
        allocated.groupby(keys)
        .agg(
            prime_mover_code=("prime_mover_code", "first"),
            **{
                column: (column, lambda values: values.sum(min_count=1))
                for column in ["net_generation_mwh", "net_generation_mwh_ozone_season"]
            },
        )
        .reset_index()
    )
    dates = (
        pd.read_csv(
            outputs_folder(f"{path_prefix}subplant_crosswalk_{year}.csv.zip"),
            usecols=keys + ["generator_operating_date", "generator_retirement_date"],
            dtype={"generator_id": "str"},
            parse_dates=["generator_operating_date", "generator_retirement_date"],
        )
        .drop_duplicates(subset=keys)
        .assign(
            generator_operating_year=lambda df: df[
                "generator_operating_date"
            ].dt.year.astype("Int64"),
            generator_retirement_year=lambda df: df[
                "generator_retirement_date"
            ].dt.year.astype("Int64"),
        )
    )
    generators = generators.merge(
        dates[keys + ["generator_operating_year", "generator_retirement_year"]],
        on=keys,
        how="left",
    ).merge(load_plant_names(path_prefix), on="plant_id_eia", how="left")
    eia860_generators = load_pudl_table(
        "core_eia860__scd_generators",
        year=year,
        columns=keys
        + ["operational_status_code", "energy_source_code_1", "capacity_mw"],
    )
    eia860_generators["operational_status_code"] = eia860_generators[
        "operational_status_code"
    ].str.upper()
    generators = generators.merge(eia860_generators, on=keys, how="left")
    return build_sheet_rows(
        generators,
        GENERATOR_DATA_SOURCES,
        year,
        "SEQGEN",
        ["state", "plant_name_eia", *keys],
    )


def load_campd_unit_attributes(year: int) -> pd.DataFrame:
    """Loads the attributes of each unit in the CAMPD facility file, translated to
    eGRID's codes.

    Args:
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia, emissions_unit_id_epa, prime mover, operating
            status, program codes, firing type, SO2 and NOx control codes, Hg control flag and year of
            commercial operation of each unit.
    """
    facility = load_campd_facility(year).drop_duplicates(
        subset=["plant_id_eia", "emissions_unit_id_epa"]
    )
    return pd.DataFrame(
        {
            "plant_id_eia": facility["plant_id_eia"],
            "emissions_unit_id_epa": facility["emissions_unit_id_epa"],
            "prime_mover_code": facility["unit_type_epa"].map(
                CAMPD_UNIT_TYPE_PRIME_MOVERS
            ),
            "operating_status": translate_campd_values(
                facility["operating_status_epa"], CAMPD_OPERATING_STATUS_CODES
            ),
            "program_code_epa": facility["program_code_epa"],
            "firing_type": translate_campd_values(
                facility["unit_type_epa"], CAMPD_UNIT_TYPE_CODES
            ),
            "so2_controls": translate_campd_values(
                facility["so2_controls_epa"], CAMPD_SO2_CONTROL_CODES
            ),
            "nox_controls": translate_campd_values(
                facility["nox_controls_epa"], CAMPD_NOX_CONTROL_CODES
            ),
            "hg_controls_flag": facility["hg_controls_epa"]
            .notna()
            .map({True: "Yes", False: pd.NA}),
            "operating_year": facility["commercial_operation_date_epa"].dt.year.astype(
                "Int64"
            ),
        }
    )


def translate_campd_values(values: pd.Series, codes: dict[str, str]) -> pd.Series:
    """Translates CAMPD descriptions to eGRID codes.

    CAMPD lists several values separated by "|" and notes when each started or ended
    in parentheses, e.g. "Water Injection (Retired Jun 01, 2024)|Selective Catalytic
    Reduction". Each value is translated, and the codes are joined by commas as in
    eGRID, e.g. "H2O,SCR". Values without a code are kept as they are.

    Args:
        values (pd.Series): CAMPD descriptions.
        codes (dict[str, str]): eGRID code of each CAMPD description.

    Returns:
        pd.Series: eGRID codes, blank where the description is blank.
    """
    return values.str.split("|").map(
        lambda descriptions: ",".join(
            codes.get(description, description)
            for description in (
                re.sub(CAMPD_DATE_NOTE, "", description).strip()
                for description in descriptions
            )
        ),
        na_action="ignore",
    )


def load_eia860_plant_attributes(year: int) -> pd.DataFrame:
    """Loads the attributes of each plant reported in EIA-860.

    Args:
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia, name, state, county, location, state and county
            FIPS codes, transmission or distribution system owner, utility, eGRID
            sector and NERC region of each plant. The county FIPS code is the three
            digits within the state, as in eGRID.
    """
    plants = load_pudl_table(
        "core_eia__entity_plants",
        columns=[
            "plant_id_eia",
            "plant_name_eia",
            "state",
            "county",
            "latitude",
            "longitude",
            "state_id_fips",
            "county_id_fips",
        ],
    ).merge(
        load_pudl_table(
            "core_eia860__scd_plants",
            year=year,
            columns=[
                "plant_id_eia",
                "transmission_distribution_owner_name",
                "transmission_distribution_owner_id",
                "utility_id_eia",
                "sector_id_eia",
                "nerc_region",
            ],
        ).merge(
            load_pudl_table("core_eia__entity_utilities"),
            on="utility_id_eia",
            how="left",
        ),
        on="plant_id_eia",
        how="left",
    )
    plants["county_id_fips"] = plants["county_id_fips"].str[-3:]
    plants["sector"] = plants["sector_id_eia"].map(EIA_SECTORS)
    return plants.drop(columns="sector_id_eia")


def load_plant_names(path_prefix: str) -> pd.DataFrame:
    """Loads the state and name of each plant from OGE's plant static attributes.

    Args:
        path_prefix (str): prefix of the results folder of the year, e.g. "2024/".

    Returns:
        pd.DataFrame: plant_id_eia, state and plant_name_eia of each plant.
    """
    return pd.read_csv(
        results_folder(f"{path_prefix}plant_data/plant_static_attributes.csv"),
        usecols=["plant_id_eia", "state", "plant_name_eia"],
    )


def build_sheet_rows(
    data: pd.DataFrame,
    sources: dict[str, str],
    year: int,
    sequence_code: str,
    sort_columns: list[str],
) -> pd.DataFrame:
    """Builds the rows of an eGRID sheet from OGE data, sorted and numbered as in eGRID.

    Args:
        data (pd.DataFrame): OGE data with one row for each row of the sheet.
        sources (dict[str, str]): OGE column for each eGRID column.
        year (int): data year.
        sequence_code (str): code of the sequence number column without the year,
            e.g. "SEQUNT".
        sort_columns (list[str]): OGE columns to sort the rows by.

    Returns:
        pd.DataFrame: the sequence number, year and columns of `sources`.
    """
    data = data.sort_values(sort_columns)
    table = pd.DataFrame(
        {code: data[column].to_numpy() for code, column in sources.items()}
    )
    table.insert(0, f"{sequence_code}{year % 100:02d}", range(1, len(table) + 1))
    table.insert(1, "YEAR", year)
    return table


def build_plant_table(year: int, path_prefix: str) -> pd.DataFrame:
    """Builds the plant sheet of the eGRID data file from OGE's annual plant data, and
    the generation by fuel, combustion heat input and ozone season totals of OGE's
    monthly subplant data.

    Args:
        year (int): data year.
        path_prefix (str): prefix of the results folder of the year, e.g. "2024/".

    Returns:
        pd.DataFrame: one row for each plant in OGE's results, in U.S. units, with the
            columns of `PLANT_DATA_SOURCES` and the sequence number and year.
    """
    plant_data = pd.read_csv(
        results_folder(f"{path_prefix}plant_data/annual/us_units/plant_data.csv")
    )
    plant_attributes = pd.read_csv(
        results_folder(f"{path_prefix}plant_data/plant_static_attributes.csv"),
        usecols=["plant_id_eia", "latitude", "longitude", "plant_primary_fuel"],
    )
    plants = (
        plant_data.merge(plant_attributes, on="plant_id_eia", how="left")
        .merge(load_ba_reference()[["ba_code", "ba_name"]], on="ba_code", how="left")
        .merge(
            load_plant_totals_from_subplants(path_prefix),
            on="plant_id_eia",
            how="left",
        )
        .merge(
            load_eia860_plant_attributes(year),
            on="plant_id_eia",
            how="left",
            suffixes=("", "_eia860"),
        )
    )
    # use the EIA-860 attributes where OGE does not have them
    for column in ["state", "plant_name_eia", "county", "latitude", "longitude"]:
        plants[column] = plants[column].fillna(plants[f"{column}_eia860"])
    table = build_sheet_rows(
        plants,
        PLANT_DATA_SOURCES,
        year,
        "SEQPLT",
        ["state", "plant_name_eia", "plant_id_eia"],
    )
    table = convert_mass_to_tons(table, "PLNT", year)
    return add_emission_rates_and_resource_mix(table, "PLNT", year)


def load_plant_totals_from_subplants(path_prefix: str) -> pd.DataFrame:
    """Sums OGE's monthly subplant data to each plant's generation by eGRID fuel group,
    combustion heat input and ozone season totals.

    Each subplant's data counts toward the fuel groups in `EGRID_FUEL_GROUPS` that
    include its fuel category. Plants without subplants in a fuel group have zero
    generation or heat input in that group.

    Args:
        path_prefix (str): prefix of the results folder of the year, e.g. "2024/".

    Returns:
        pd.DataFrame: plant_id_eia, the net generation of each fuel group of
            `GENERATION_BY_FUEL_CODES` (e.g. "net_generation_mwh_coal"), the
            combustion heat input (e.g. "combustion_fuel_consumed_mmbtu"), and the
            ozone season totals of heat input, net generation, NOx and combustion heat
            input (with "_ozone_season" added to their column names).
    """
    heat_input_columns = ["fuel_consumed_mmbtu", "fuel_consumed_for_electricity_mmbtu"]
    ozone_season_columns = heat_input_columns + [
        "net_generation_mwh",
        "nox_mass_lb",
        "nox_mass_lb_for_electricity_adjusted",
    ]
    subplants = pd.read_csv(
        results_folder(f"{path_prefix}plant_data/monthly/us_units/subplant_data.csv"),
        usecols=["plant_id_eia", "subplant_id", "report_date"] + ozone_season_columns,
        parse_dates=["report_date"],
    ).merge(
        pd.read_csv(
            results_folder(f"{path_prefix}plant_data/subplant_attributes.csv"),
            usecols=["plant_id_eia", "subplant_id", "fuel_category"],
        ),
        on=["plant_id_eia", "subplant_id"],
        how="left",
        validate="m:1",
    )
    generation_columns = []
    for group in dict.fromkeys(GENERATION_BY_FUEL_CODES.values()):
        in_group = subplants["fuel_category"].isin(EGRID_FUEL_GROUPS[group])
        subplants[f"net_generation_mwh_{group}"] = subplants[
            "net_generation_mwh"
        ].where(in_group, 0)
        generation_columns.append(f"net_generation_mwh_{group}")
    combustion = subplants["fuel_category"].isin(EGRID_FUEL_GROUPS["combustion"])
    for column in heat_input_columns:
        subplants[f"combustion_{column}"] = subplants[column].where(combustion, 0)
    ozone_season = subplants["report_date"].dt.month.isin(OZONE_SEASON_MONTHS)
    seasonal_columns = (
        ozone_season_columns
        + [f"combustion_{column}" for column in heat_input_columns]
        + ["net_generation_mwh_combustion"]
    )
    for column in seasonal_columns:
        subplants[f"{column}_ozone_season"] = subplants[column].where(ozone_season, 0)
    return (
        subplants.groupby("plant_id_eia")[
            generation_columns
            + [f"combustion_{column}" for column in heat_input_columns]
            + [f"{column}_ozone_season" for column in seasonal_columns]
        ]
        .sum(min_count=1)
        .reset_index()
    )


def build_region_table(
    year: int, path_prefix: str, sheet: str, ba_codes: list[str]
) -> pd.DataFrame:
    """Builds the BA or U.S. sheet of the eGRID data file from OGE's annual and monthly
    power sector results, including the emission rates by fuel type.

    Generation and heat input from fuel groups that a BA does not have are zero.

    Args:
        year (int): data year.
        path_prefix (str): prefix of the results folder of the year, e.g. "2024/".
        sheet (str): "BA" or "US".
        ba_codes (list[str]): name of each power sector results file to read, i.e.
            each BA code, or "US".

    Returns:
        pd.DataFrame: one row for each file, in U.S. units, with the year, the BA name
            and code for the BA sheet, and the columns of `REGION_DATA_SOURCES`.
    """
    prefix = REGION_SHEETS[sheet][0]
    rows = []
    for ba_code in ba_codes:
        folder = f"{path_prefix}power_sector_data/{{resolution}}/us_units/{ba_code}.csv"
        monthly = pd.read_csv(
            results_folder(folder.format(resolution="monthly")),
            parse_dates=["report_date"],
        )
        ba_data = {
            "annual": pd.read_csv(
                results_folder(folder.format(resolution="annual"))
            ).set_index("fuel_category"),
            "ozone": monthly[monthly["report_date"].dt.month.isin(OZONE_SEASON_MONTHS)]
            .groupby("fuel_category")
            .sum(numeric_only=True, min_count=1),
        }
        row = {"YEAR": year}
        for code, (fuel_category, column, season) in REGION_DATA_SOURCES.items():
            fuel_categories = EGRID_FUEL_GROUPS.get(fuel_category, [fuel_category])
            data = ba_data[season]
            if data.index.isin(fuel_categories).any():
                row[prefix + code] = data.loc[
                    data.index.isin(fuel_categories), column
                ].sum(min_count=1)
            elif fuel_category != "total":
                # the BA has none of the fuel categories
                row[prefix + code] = 0
        row |= {
            prefix + code: rate
            for code, rate in calculate_fuel_type_rates(
                ba_data["annual"], ba_data["ozone"]
            ).items()
        }
        if sheet == "BA":
            row["BACODE"] = ba_code
        rows.append(row)
    table = pd.DataFrame(rows)
    if sheet == "BA":
        table = (
            table.merge(
                load_ba_reference()[["ba_code", "ba_name"]].rename(
                    columns={"ba_code": "BACODE", "ba_name": "BANAME"}
                ),
                on="BACODE",
                how="left",
            )
            .sort_values("BANAME")
            .reset_index(drop=True)
        )
    table = convert_mass_to_tons(table, sheet, year)
    return add_emission_rates_and_resource_mix(table, sheet, year)


def add_emission_rates_and_resource_mix(
    table: pd.DataFrame, sheet: str, year: int
) -> pd.DataFrame:
    """Adds the emission rates and resource mix to a plant, BA or U.S. sheet, from
    the sheet's emissions, heat input and generation.

    Rates are in lb per MWh or MMBtu, and are blank where the denominator is zero or
    negative, as in eGRID. Each fuel's share of the resource mix is its positive
    generation divided by the sum of the positive generation of every fuel, so that
    the shares of each row add up to 100%.

    Args:
        table (pd.DataFrame): table of the sheet with masses in the units of eGRID's
            columns.
        sheet (str): "PLNT", "BA" or "US".
        year (int): data year.

    Returns:
        pd.DataFrame: the table with the emission rate and resource mix columns.
    """
    table = table.copy()
    prefix = "PL" if sheet == "PLNT" else REGION_SHEETS[sheet][0]
    units = {
        code: get_unit(header) for code, header, _ in get_egrid_columns(sheet, year)
    }

    def mass_lb(code: str) -> pd.Series:
        factor = ConversionFactors.short_ton_to_lbs if units[code] == "tons" else 1
        return table[code] * factor

    def rate(mass: pd.Series, denominator: pd.Series) -> pd.Series:
        return mass / denominator.where(denominator > 0)

    for rate_code, (
        annual,
        ozone_season,
        ozone_rate_code,
    ) in EMISSION_RATE_DENOMINATORS.items():
        for pollutant, mass_code in EMISSION_RATE_MASSES.items():
            table[f"{prefix}{pollutant}{rate_code}"] = rate(
                mass_lb(prefix + mass_code), table[prefix + annual]
            )
        table[f"{prefix}NOX{ozone_rate_code}"] = rate(
            mass_lb(f"{prefix}NOXOZ"), table[prefix + ozone_season]
        )

    total = (
        table[[prefix + code for code in RESOURCE_MIX_FUELS]]
        .clip(lower=0)
        .sum(axis=1, min_count=1)
    )
    for mix_code, generation_code in RESOURCE_MIX_GENERATION.items():
        table[prefix + mix_code] = table[prefix + generation_code].clip(
            lower=0
        ) / total.where(total > 0)
    return table.drop(columns=f"{prefix}{OZONE_COMBUSTION_GENERATION}")


def calculate_fuel_type_rates(
    annual: pd.DataFrame, ozone_season: pd.DataFrame
) -> dict[str, float]:
    """Calculates the emission rates by fuel type of a BA from its power sector results.

    Each rate divides the emissions of the BA's fuel category rows in a fuel group by
    their net generation (output rates, lb/MWh) or heat input (input rates,
    lb/MMBtu). Rates are blank where the BA has none of the fuel group's categories or
    the denominator is zero or negative.

    Args:
        annual (pd.DataFrame): annual power sector results, indexed by fuel category.
        ozone_season (pd.DataFrame): ozone season sums of the monthly power sector
            results, indexed by fuel category.

    Returns:
        dict[str, float]: each rate, keyed by its code without the sheet's prefix.
    """
    rates = {}
    for pollutant, codes in FUEL_TYPE_OUTPUT_RATE_CODES.items():
        data = ozone_season if pollutant == "NOXOZ" else annual
        mass_column = FUEL_TYPE_RATE_MASSES[pollutant.removesuffix("OZ")]
        for group, output_code in codes.items():
            rows = data.loc[data.index.isin(EGRID_FUEL_GROUPS[group])]
            mass = rows[mass_column].sum(min_count=1)
            for code, denominator_column in [
                (output_code, "net_generation_mwh"),
                (output_code.removesuffix("T"), "fuel_consumed_for_electricity_mmbtu"),
            ]:
                denominator = rows[denominator_column].sum(min_count=1)
                rates[code] = mass / denominator if denominator > 0 else np.nan
    return rates


def convert_mass_to_tons(table: pd.DataFrame, sheet: str, year: int) -> pd.DataFrame:
    """Converts OGE masses in lb to short tons, in the columns of an eGRID sheet whose
    headers are in tons.

    Args:
        table (pd.DataFrame): table of an eGRID sheet, with masses in lb.
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        year (int): data year.

    Returns:
        pd.DataFrame: the table with masses in the units of the eGRID sheet.
    """
    table = table.copy()
    for code, header, _ in get_egrid_columns(sheet, year):
        if code in table and get_unit(header) == "tons":
            table[code] = table[code] / ConversionFactors.short_ton_to_lbs
    return table


def convert_to_metric(table: pd.DataFrame, sheet: str, year: int) -> pd.DataFrame:
    """Converts a table of an eGRID sheet from U.S. units to metric units.

    Args:
        table (pd.DataFrame): table of an eGRID sheet in U.S. units, with the columns
            of `get_egrid_columns()`.
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        year (int): data year.

    Returns:
        pd.DataFrame: the table in metric units, with the columns of
            `get_egrid_columns(metric=True)`.
    """
    metric = {}
    for code, header, _ in get_egrid_columns(sheet, year):
        unit = get_unit(header)
        if unit in METRIC_UNITS:
            metric[code] = table[code] * METRIC_UNITS[unit][1]
        else:
            metric[code] = table[code]
        if unit in GJ_UNITS:
            metric[f"{code}2"] = table[code] * GJ_UNITS[unit][1]
    return pd.DataFrame(metric, index=table.index)


def write_egrid_workbook(
    tables: dict[str, pd.DataFrame], year: int, path: str, metric: bool
):
    """Writes the tables of each sheet to an Excel workbook in the layout of the eGRID
    data file, after a contents sheet.

    The first row of each sheet holds the column headers and the second row the column
    codes, colored by category. The headers of columns with no data are grayed out.

    Args:
        tables (dict[str, pd.DataFrame]): table of each sheet, keyed by its name in
            `EGRID_SHEETS`, with the columns of `get_egrid_columns()`.
        year (int): data year.
        path (str): path of the workbook.
        metric (bool): whether the tables are in metric units.
    """
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        formats = WorkbookFormats(writer.book)
        write_contents_sheet(writer.book, formats, year, metric)
        for sheet in EGRID_SHEETS:
            sheet_name = f"{sheet}{year % 100:02d}"
            table = tables[sheet]
            table.to_excel(
                writer, sheet_name=sheet_name, startrow=2, header=False, index=False
            )
            worksheet = writer.sheets[sheet_name]
            number_formats = get_number_formats(sheet, year, metric)
            for i, (code, header, category) in enumerate(
                get_egrid_columns(sheet, year, metric)
            ):
                calculated = table[code].notna().any()
                worksheet.write(0, i, header, formats.header(category, calculated))
                worksheet.write(
                    1, i, code, formats.header(category, calculated, code_row=True)
                )
                worksheet.set_column(
                    i,
                    i,
                    get_column_width(sheet, code),
                    formats.get({"num_format": number_formats[code]}),
                )
            # filter on the column codes
            worksheet.autofilter(1, 0, len(table) + 1, len(table.columns) - 1)
            worksheet.set_row(0, HEADER_HEIGHTS.get(sheet, DEFAULT_HEADER_HEIGHT))
            if FROZEN_COLUMNS[sheet] is not None:
                worksheet.freeze_panes(2, FROZEN_COLUMNS[sheet])


class WorkbookFormats:
    """Creates the formats of an eGRID workbook, reusing each distinct format.

    Every format uses the Arial font at eGRID's size unless the properties set another.
    """

    FONT = {"font_name": "Arial", "font_size": 8.5}

    def __init__(self, workbook: xlsxwriter.Workbook):
        self.workbook = workbook
        self.formats = {}

    def get(self, properties: dict) -> xlsxwriter.format.Format:
        """Returns the format with these properties, creating it if needed."""
        properties = {**self.FONT, **properties}
        key = tuple(sorted(properties.items()))
        if key not in self.formats:
            self.formats[key] = self.workbook.add_format(properties)
        return self.formats[key]

    def header(
        self, category: int | None, calculated: bool = True, code_row: bool = False
    ) -> xlsxwriter.format.Format:
        """Returns the format of a column header or code, colored by the column's
        category. The headers of columns that OGE does not calculate are grayed out.
        """
        properties = {
            "bold": True,
            "text_wrap": True,
            "border": 1,
            "border_color": "#000000",
        }
        if not code_row:
            properties |= {"align": "center", "valign": "vcenter"}
        if not calculated:
            properties |= {"bg_color": "#FFFFFF", "font_color": "#A6A6A6"}
        else:
            properties |= {
                "bg_color": CATEGORY_COLORS[category],
                "font_color": "#FFFFFF" if category in DARK_CATEGORIES else "#000000",
            }
        return self.get(properties)


def get_number_formats(sheet: str, year: int, metric: bool = False) -> dict[str, str]:
    """Returns the Excel number format of each column of a sheet of the eGRID data
    file, following the formats of the eGRID data file.

    The GJ columns of the metric file have the format of the column they follow.

    Args:
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        year (int): data year.
        metric (bool, optional): whether to return the formats of the metric units
            file. Defaults to False.

    Returns:
        dict[str, str]: number format of each column code. Columns that are not
            numbers have the "General" format.
    """
    formats = {}
    for code, header, category in get_egrid_columns(sheet, year):
        units = re.findall(r"\(([^()]+)\)", header)
        unit = units[-1] if units else None
        if code in NUMBER_FORMATS and sheet not in REGION_SHEETS:
            number_format = NUMBER_FORMATS[code]
        elif unit == "resource mix" or code == "GGRSLOSS":
            number_format = PERCENT_FORMAT
        elif unit not in [*METRIC_UNITS, *GJ_UNITS, "MW", "Btu/kWh"] or unit == "feet":
            number_format = "General"
        elif "/" in unit and unit != "Btu/kWh":
            small = category in [13, 14] and ("CH4" in code or "N2O" in code)
            number_format = SMALL_RATE_FORMAT if small else RATE_FORMAT
        elif sheet == "UNT":
            number_format = UNIT_MASS_FORMAT if unit in ["tons", "lbs"] else "#,##0"
        elif sheet == "PLNT" and "Hg" in header:
            number_format = HG_MASS_FORMAT
        elif sheet == "GGL" or (sheet == "PLNT" and category == 2):
            number_format = "#,##0"
        else:
            number_format = QUANTITY_FORMAT
        formats[code] = number_format
        if metric and get_unit(header) in GJ_UNITS:
            formats[f"{code}2"] = number_format
    return formats


def get_column_width(sheet: str, code: str) -> float:
    """Returns the width of a column of a sheet of the eGRID data file.

    Args:
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        code (str): column code.

    Returns:
        float: column width, in Excel's character units.
    """
    if code in COLUMN_WIDTHS[sheet]:
        return COLUMN_WIDTHS[sheet][code]
    if sheet in REGION_SHEETS:
        return REGION_DATA_WIDTHS.get(
            code.removeprefix(REGION_SHEETS[sheet][0]).removesuffix("2"),
            DEFAULT_COLUMN_WIDTH,
        )
    return DEFAULT_COLUMN_WIDTH


def write_contents_sheet(
    workbook: xlsxwriter.Workbook,
    formats: WorkbookFormats,
    year: int,
    metric: bool,
):
    """Adds a contents sheet to an eGRID workbook, laid out like the contents sheet of
    the eGRID data file.

    The sheet has a table of contents, links to OGE's documentation, a legend of the
    header colors with a link to each category on each sheet, notes and conversion
    factors. Each section is a gray box with a thick outline. The table of contents
    uses the eGRID data file's text, except that it says the file was produced with
    OGE.

    Args:
        workbook (xlsxwriter.Workbook): workbook to add the sheet to.
        formats (WorkbookFormats): formats of the workbook.
        year (int): data year.
        metric (bool): whether the workbook is in metric units.
    """
    worksheet = workbook.add_worksheet("Contents")
    yy = f"{year % 100:02d}"
    link = {"underline": 1, "font_color": "#0563C1"}

    title = {"font_size": 14, "bold": True, "align": "center", "text_wrap": True}
    heading = {"font_size": 12, "bold": True}
    large = {"font_size": 12}

    # table of contents
    sheets = [
        [
            (1, 1, i, {**large, "align": "center"}),
            (
                2,
                3,
                f"internal:'{sheet}{yy}'!A1",
                {**large, **link, "string": f"{sheet}{yy}"},
            ),
            (4, 10, f"{EGRID_SHEET_DESCRIPTIONS[sheet]} year {year} data", large),
        ]
        for i, sheet in enumerate(EGRID_SHEETS, start=1)
    ]
    write_contents_box(
        worksheet,
        formats,
        1,
        [
            [(1, 10, "Table of Contents", title)],
            [(1, 10, "for", title)],
            [(1, 10, f"eGRID{year}_Data{'_Metric' if metric else ''}.xlsx", title)],
            [
                (
                    1,
                    10,
                    f"eGRID{year} {describe_sheets()} Data Files",
                    title,
                )
            ],
            [
                (
                    1,
                    10,
                    date.today().strftime("%B %d, %Y"),
                    {**heading, "align": "center"},
                )
            ],
            [
                (1, 4, "Produced with", {**heading, "align": "right", "valign": "top"}),
                (
                    5,
                    10,
                    OGE_REPOSITORY_URL,
                    {
                        **large,
                        **link,
                        "valign": "top",
                        "string": f"Open Grid Emissions {version('oge')}.",
                    },
                ),
            ],
            [
                (1, 1, "Sheet", {**large, "underline": 1, "align": "center"}),
                (2, 3, "Name", {**large, "underline": 1}),
                (4, 10, "Description", {**large, "underline": 1}),
            ],
            *sheets,
            [(2, 10, "", {})],
        ],
        last_col=10,
    )
    worksheet.set_row(4, 43.9)
    for row in [5, 6]:
        worksheet.set_row(row, 25.9)

    # documentation and feedback
    feedback_row = 1 + len(sheets) + 9
    write_contents_box(
        worksheet,
        formats,
        feedback_row,
        [
            [(1, 10, "Documentation and Feedback", heading)],
            [(1, 10, OGE_DOCUMENTATION_URL, {**link, "string": "OGE documentation"})],
            [
                (
                    1,
                    10,
                    OGE_ISSUES_URL,
                    {**link, "string": "Report an issue or ask a question"},
                )
            ],
        ],
        last_col=10,
    )

    # color coding legend, with a link to the first column of each category on each
    # sheet
    legend_row = feedback_row + 4
    first_columns = {}
    for col, sheet in enumerate(EGRID_SHEETS[:-1], start=9):
        for i, (_, _, category) in enumerate(get_egrid_columns(sheet, year, metric)):
            if category is not None and (category, col) not in first_columns:
                first_columns[(category, col)] = (sheet, i)
    legend = [
        [
            (1, 7, "Category", {**heading, "bottom": 2}),
            (8, 8, "Color", {**heading, "align": "center", "bottom": 2}),
            (
                9,
                16,
                "Link to sheet and category",
                {**heading, "align": "center", "bottom": 2},
            ),
        ]
    ]
    thin = {"left": 1, "right": 1, "top": 1, "bottom": 1}
    for category, name in EGRID_CATEGORIES.items():
        row = [
            (1, 7, f"{category}) {name}", thin),
            (8, 8, "", {**thin, "bg_color": CATEGORY_COLORS[category]}),
        ]
        for col in range(9, 17):
            if (category, col) in first_columns:
                sheet, i = first_columns[(category, col)]
                cell = xl_rowcol_to_cell(0, i)
                row.append(
                    (
                        col,
                        col,
                        f"internal:'{sheet}{yy}'!{cell}",
                        {**thin, **link, "align": "center", "string": sheet},
                    )
                )
            else:
                row.append((col, col, "", thin))
        legend.append(row)
    legend.append(
        [
            (1, 7, "Not calculated by OGE yet", thin),
            (8, 8, "", {**thin, "bg_color": "#FFFFFF"}),
            *[(col, col, "", thin) for col in range(9, 17)],
        ]
    )
    write_contents_box(
        worksheet,
        formats,
        legend_row,
        [[(1, 16, "Color Coding Legend", {**heading, "align": "center"})], *legend],
        last_col=16,
    )

    # notes
    notes_row = legend_row + len(legend) + 2
    notes = [
        "Values in parentheses are negative numbers.",
        "Rows are the plants and balancing authorities in OGE's results. Plants are "
        "identified by their EIA plant ID (ORISPL).",
        "Columns that OGE does not calculate yet are blank, and their headers are "
        "grayed out.",
    ]
    write_contents_box(
        worksheet,
        formats,
        notes_row,
        [[(1, 10, "Notes", heading)], *[[(1, 10, note, {})] for note in notes]],
        last_col=10,
    )

    # conversion factors
    factors_row = notes_row + len(notes) + 2
    factors = [
        ("1 megawatt-hour (MWh)", "1,000 kilowatt-hour (kWh)"),
        ("1 short ton", "2,000 pounds (lb)"),
        ("1.10231 short ton", "1 metric ton"),
        ("2.2046 lb", "1 kilogram (kg)"),
        ("0.9478 MMBtu", "1 Gigajoule (GJ)"),
        ("0.27778 MWh", "1 GJ"),
    ]
    write_contents_box(
        worksheet,
        formats,
        factors_row,
        [
            [(1, 5, "Conversion Factors", heading)],
            *[[(1, 2, left, {}), (3, 4, right, {})] for left, right in factors],
        ],
        last_col=5,
    )

    for col, width in {
        0: 1.73,
        1: 15,
        2: 11,
        4: 12.54,
        5: 18.73,
        6: 11,
        9: 5.73,
    }.items():
        worksheet.set_column(col, col, width)
    worksheet.set_column(10, 16, 6)


def describe_sheets() -> str:
    """Returns a list of the sheets of the eGRID data file for the contents sheet's
    title, e.g. "Unit, Generator, ..., and Grid Gross Loss (%)", in title case.
    """
    descriptions = [
        " ".join(
            word.capitalize() if word.islower() else word
            for word in description.split()
        )
        for description in EGRID_SHEET_DESCRIPTIONS.values()
    ]
    return f"{', '.join(descriptions[:-1])}, and {descriptions[-1]}"


def write_contents_box(
    worksheet: xlsxwriter.worksheet.Worksheet,
    formats: WorkbookFormats,
    first_row: int,
    rows: list[list[tuple]],
    last_col: int,
):
    """Writes a section of the contents sheet as a gray box with a thick outline,
    starting at column B.

    Args:
        worksheet (xlsxwriter.worksheet.Worksheet): contents sheet.
        formats (WorkbookFormats): formats of the workbook.
        first_row (int): first row of the box, from 0.
        rows (list[list[tuple]]): cells of each row of the box, as (first column, last
            column, value, format properties), merged across their columns. Values
            that start with "internal:" or "https:" are written as links, whose text
            is the "string" property.
        last_col (int): last column of the box, from 0.
    """
    last_row = first_row + len(rows) - 1

    def cell_format(row: int, first_col: int, cell_last_col: int, properties: dict):
        properties = {
            "font_size": 11,
            "bg_color": CATEGORY_COLORS[None],
            "valign": "vcenter",
            "border_color": "#000000",
            **properties,
        }
        for side, on_edge in {
            "left": first_col == 1,
            "right": cell_last_col == last_col,
            "top": row == first_row,
            "bottom": row == last_row,
        }.items():
            if on_edge:
                properties[side] = 2
        return formats.get(properties)

    for row, cells in enumerate(rows, start=first_row):
        for col in range(1, last_col + 1):
            worksheet.write_blank(row, col, None, cell_format(row, col, col, {}))
        for first_col, cell_last_col, value, properties in cells:
            properties = dict(properties)
            string = properties.pop("string", None)
            merged_format = cell_format(row, first_col, cell_last_col, properties)
            if cell_last_col > first_col:
                worksheet.merge_range(
                    row, first_col, row, cell_last_col, "", merged_format
                )
            if isinstance(value, str) and value.startswith(("internal:", "https:")):
                worksheet.write_url(row, first_col, value, merged_format, string=string)
            else:
                worksheet.write(row, first_col, value, merged_format)


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
    metric_units = {unit: name for unit, (name, _) in METRIC_UNITS.items()}
    gj_units = {unit: name for unit, (name, _) in GJ_UNITS.items()}
    metric_columns = []
    for code, header, category in columns:
        metric_columns.append((code, replace_units(header, metric_units), category))
        if get_unit(header) in GJ_UNITS:
            metric_columns.append(
                (f"{code}2", replace_units(header, gj_units), category)
            )
    return metric_columns


def get_unit(header: str) -> str | None:
    """Returns the U.S. unit of a column header that eGRID converts in its metric file,
    e.g. "tons" for "Plant annual CO2 emissions (tons)".

    Args:
        header (str): column header.

    Returns:
        str | None: the unit in parentheses that is a key of `METRIC_UNITS` or
            `GJ_UNITS`, or None if there is none.
    """
    for unit in re.findall(r"\(([^()]+)\)", header):
        if unit in METRIC_UNITS or unit in GJ_UNITS:
            return unit
    return None


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
