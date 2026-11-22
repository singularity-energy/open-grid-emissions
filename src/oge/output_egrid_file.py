"""Writes OGE annual results in the layout of EPA's eGRID data file.

The tables of each sheet are built from OGE's annual results and intermediate outputs,
using OGE's plants, balancing authorities and fuel categories. Columns that OGE does not
calculate yet are left blank. The sheets and columns of the eGRID data file are defined
in `oge.egrid_schema`, and the workbooks are formatted and written by
`oge.egrid_workbook`.
"""

import os
import re

import numpy as np
import pandas as pd

from oge.constants import ConversionFactors
from oge.egrid_schema import (
    EGRID_SHEETS,
    EMISSION_RATE_DENOMINATORS,
    EMISSION_RATE_MASSES,
    FUEL_TYPE_OUTPUT_RATE_CODES,
    GENERATION_BY_FUEL_CODES,
    GJ_UNITS,
    METRIC_UNITS,
    NONBASELOAD_GENERATION_CODES,
    NONBASELOAD_OZONE_NOX_RATE,
    NONBASELOAD_RATE_PREFIX,
    NONBASELOAD_RESOURCE_MIX_GENERATION,
    OZONE_COMBUSTION_GENERATION,
    OZONE_SEASON_MONTHS,
    PLANT_SOURCE_VALUES,
    REGION_SHEETS,
    RESOURCE_MIX_FUELS,
    RESOURCE_MIX_GENERATION,
    RESOURCE_MIX_GROUP_FUELS,
    get_egrid_columns,
    get_unit,
)
from oge.egrid_workbook import write_egrid_workbook
from oge.filepaths import outputs_folder, results_folder
from oge.helpers import assign_fuel_category_to_esc
from oge.load_data import (
    load_ba_reference,
    load_campd_facility,
    load_pudl_table,
    load_raw_eia860_stack_heights,
    load_unit_to_boiler_associations,
)
from oge.logging_util import get_logger
from oge.output_data import calculate_capacity_factor

logger = get_logger(__name__)

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
# eGRID label of each OGE data source, used in the source columns of the unit and
# plant sheets
DATA_SOURCE_LABELS = {"CEMS": "EPA/CAPD", "EIA": "EIA"}

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
    # tallest operating stack of the unit's boilers, from EIA-860
    "STACKHT": "stack_height_ft",
    # generators of the unit in OGE's subplant crosswalk
    "NUMGEN": "generator_count",
    # sources of the emissions, from the CEMS measurement codes, and of the heat input
    "HTIANSRC": "fuel_consumed_mmbtu_source",
    "HTIOZSRC": "fuel_consumed_mmbtu_ozone_season_source",
    "NOXANSRC": "nox_mass_lb_source",
    "NOXOZSRC": "nox_mass_lb_ozone_season_source",
    "SO2SRC": "so2_mass_lb_source",
    "CO2SRC": "co2_mass_lb_source",
    "HGSRC": "hg_mass_lb_source",
}
# column of the CEMS measurement codes of each emissions column in the cleaned CEMS
# data, which are the sources of the unit sheet's emissions
UNIT_MEASUREMENT_CODES = {
    "nox_mass_lb": "nox_mass_measurement_code",
    "nox_mass_lb_ozone_season": "nox_mass_measurement_code",
    "so2_mass_lb": "so2_mass_measurement_code",
    "co2_mass_lb": "co2_mass_measurement_code",
    "hg_mass_lb": "hg_mass_measurement_code",
}
# EIA-860 statuses of the stacks whose heights are used for each unit's stack height.
# As in eGRID, retired, cancelled, planned and out-of-service stacks are left out
OPERATING_STACK_STATUS_CODES = ["OP", "SB"]
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
    "Internal combustion engine": "IC",
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
    "ISORTO": "iso_rto",
    "PLFUELCT": "fuel_category_egrid",
    "COALFLAG": "coal_flag",
    "PSFLAG": "pumped_storage_flag",
    "NAMEPCAP": "capacity_mw",
    "CAPFAC": "capacity_factor",
    "NBFACTOR": "nonbaseload_factor",
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
    # adjustment values: the biomass adjustment is OGE's unadjusted value minus its
    # biomass-adjusted value, and the CHP adjustment is the biomass-adjusted value
    # minus the biomass- and CHP-adjusted value
    "BIONOX": "nox_mass_lb_biomass_adjustment",
    "BIONOXOZ": "nox_mass_lb_biomass_adjustment_ozone_season",
    "BIOSO2": "so2_mass_lb_biomass_adjustment",
    "BIOCO2": "co2_mass_lb_biomass_adjustment",
    "BIOCH4": "ch4_mass_lb_biomass_adjustment",
    "BION2O": "n2o_mass_lb_biomass_adjustment",
    "BIOCO2E": "co2e_mass_lb_biomass_adjustment",
    "CHPCHTI": "combustion_fuel_consumed_mmbtu_chp_adjustment",
    "CHPCHTIOZ": "combustion_fuel_consumed_mmbtu_chp_adjustment_ozone_season",
    "CHPNOX": "nox_mass_lb_chp_adjustment",
    "CHPNOXOZ": "nox_mass_lb_chp_adjustment_ozone_season",
    "CHPSO2": "so2_mass_lb_chp_adjustment",
    "CHPCO2": "co2_mass_lb_chp_adjustment",
    "CHPCH4": "ch4_mass_lb_chp_adjustment",
    "CHPN2O": "n2o_mass_lb_chp_adjustment",
    "CHPCO2E": "co2e_mass_lb_chp_adjustment",
    "RMBMFLAG": "biomass_adjustment_flag",
    "CHPFLAG": "chp_adjustment_flag",
    # sources of the unadjusted values, from OGE's plant metadata
    "UNNOXSRC": "nox_mass_lb_source",
    "UNNOZSRC": "nox_mass_lb_ozone_season_source",
    "UNSO2SRC": "so2_mass_lb_source",
    "UNCO2SRC": "co2_mass_lb_source",
    "UNCH4SRC": "ch4_mass_lb_source",
    "UNN2OSRC": "n2o_mass_lb_source",
    "UNC2ESRC": "co2e_mass_lb_source",
    "UNHGSRC": "hg_mass_lb_source",
    "UNHTISRC": "fuel_consumed_mmbtu_source",
    "UNHOZSRC": "fuel_consumed_mmbtu_ozone_season_source",
} | {
    f"PL{code}": f"net_generation_mwh_{group}"
    for code, group in GENERATION_BY_FUEL_CODES.items()
}
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
# eGRID's ISO/RTO of each OGE BA that is an ISO or RTO
BA_ISO_RTO = {
    "CISO": "CAISO",
    "ERCO": "ERCOT",
    "ISNE": "ISONE",
    "MISO": "MISO",
    "NYIS": "NYISO",
    "PJM": "PJM",
    "SWPP": "SPP",
}
# eGRID's name of each NERC region that EIA-860 reports for plants since 2005, including
# historical regions that were later merged or dissolved. The names of ECAR, MAAC, MAIN
# and MAPP are from the eGRID 1996-2000 Technical Support Document, and the others from
# eGRID's NERC region sheets. Alaska, Hawaii and Puerto Rico are not in a NERC region,
# so, as in eGRID, their plants' region is their state (since 2020 for Alaska and
# Hawaii, which were ASCC and HICC before)
NERC_REGION_NAMES = {
    "AK": "Alaska",
    "ASCC": "Alaska Systems Coordinating Council",  # dissolved 2019 -> AK
    "ECAR": "East Central Area Reliability Coordination Agreement",  # dissolved 2005 -> RF
    "ERCOT": "Electric Reliability Council of Texas",  # legacy name for TRE
    "FRCC": "Florida Reliability Coordinating Council",  # dissolved 2019 -> SERC
    "HI": "Hawaii",
    "HICC": "Hawaiian Islands Coordinating Council",  # dissolved 2019 -> HI
    "MAAC": "Mid-Atlantic Area Council",  # dissolved 2005 -> RF
    "MAIN": "Mid-America Interconnected Network",  # dissolved 2005
    "MAPP": "Mid-Continent Area Power Pool",  # dissolved 2004 -> MRO
    "MRO": "Midwest Reliability Organization",
    "NPCC": "Northeast Power Coordinating Council",
    "PR": "Puerto Rico",
    "RFC": "Reliability First Corporation",  # aka RF / "ReliabilityFirst"
    "SERC": "SERC Reliability Corporation",
    "SPP": "Southwest Power Pool",  # legacy
    "TRE": "Texas Regional Entity",
    "WECC": "Western Electricity Coordinating Council",
}
# plant static attribute and identifying column of the regional sheets that are summed
# from plants
PLANT_REGION_ATTRIBUTES = {"ST": ("state", "PSTATABB"), "NRL": ("nerc_region", "NERC")}
# eGRID's plant primary fuel category of each OGE plant fuel category. As in the fuel
# columns, other fossil is OGE's waste and other is OGE's other and storage (DAT-1067)
PLANT_FUEL_CATEGORIES = {
    "coal": "COAL",
    "petroleum": "OIL",
    "natural_gas": "GAS",
    "nuclear": "NUCLEAR",
    "hydro": "HYDRO",
    "biomass": "BIOMASS",
    "wind": "WIND",
    "solar": "SOLAR",
    "geothermal": "GEOTHERMAL",
    "waste": "OFSL",
    "other": "OTHF",
    "storage": "OTHF",
}
# prime mover of pumped storage subplants
PUMPED_STORAGE_PRIME_MOVER = "PS"
# pollutants whose plant adjustment values are filled, as in OGE's column names
ADJUSTED_POLLUTANTS = ["nox", "so2", "co2", "ch4", "n2o", "co2e"]

# OGE fuel category or eGRID fuel group, column and season of each data column of the
# state, BA, NERC region and U.S. sheets that holds the same quantity, keyed by the code
# without the sheet's prefix. A fuel group sums its fuel categories in
# EGRID_FUEL_GROUPS, and the ozone season sums May to September from the monthly results
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


########################################################################################
# Entry points
########################################################################################


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
        "ST": build_plant_region_table(year, path_prefix, "ST"),
        "NRL": build_plant_region_table(year, path_prefix, "NRL"),
        "BA": build_region_table(year, path_prefix, "BA", ba_codes),
        "US": build_region_table(year, path_prefix, "US", ["US"]),
    }
    tables = add_regional_capacity(tables)
    tables = add_regional_nonbaseload_values(tables, year)
    tables = add_plant_unit_values(tables)
    return {
        sheet: tables.get(sheet, pd.DataFrame()).reindex(
            columns=[code for code, _, _ in get_egrid_columns(sheet, year)]
        )
        for sheet in EGRID_SHEETS
    }


########################################################################################
# Unit sheet
########################################################################################


def build_unit_table(year: int, path_prefix: str) -> pd.DataFrame:
    """Builds the unit sheet of the eGRID data file from OGE's cleaned hourly CEMS data
    and the unit attributes in the CAMPD facility file.

    The sources of the emissions are the CEMS measurement codes of each unit, e.g.
    "Measured; Substitute", listed from the code covering the most emissions to the
    least.

    The sheet has the units that report to CEMS. Hg is the sum of the hours with Hg
    data, as in eGRID, and is blank for units that do not report Hg. Stack heights come
    from EIA-860 (see `load_unit_stack_heights`).

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
        + list(UNIT_MEASUREMENT_CODES.values())
        + ["energy_source_code", "operating_time_hours", "report_date"],
        dtype={"emissions_unit_id_epa": "str", "energy_source_code": "str"}
        | {code_column: "str" for code_column in UNIT_MEASUREMENT_CODES.values()},
        parse_dates=["report_date"],
    )
    # ozone season NOx, used to list the measurement codes of the ozone season
    cems["nox_mass_lb_ozone_season"] = cems["nox_mass_lb"].where(
        cems["report_date"].dt.month.isin(OZONE_SEASON_MONTHS)
    )
    # hours in which a unit did not operate have no heat input or NOx, so their blanks
    # are zeros in the ozone season totals, which are otherwise blank when unknown
    seasonal_columns = ["fuel_consumed_mmbtu", "nox_mass_lb"]
    seasonal_data = cems.copy()[
        keys + ["report_date", "operating_time_hours"] + seasonal_columns
    ]
    not_operating = seasonal_data["operating_time_hours"] == 0
    seasonal_data.loc[not_operating, seasonal_columns] = seasonal_data.loc[
        not_operating, seasonal_columns
    ].fillna(0)
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
        .merge(
            sum_ozone_season_values(seasonal_data, keys, seasonal_columns),
            on=keys,
            how="left",
            validate="1:1",
        )
    )
    # the sources of each unit's heat input and emissions. Heat input from CEMS has no
    # measurement code, so its source is CEMS
    for column in ["fuel_consumed_mmbtu", "fuel_consumed_mmbtu_ozone_season"]:
        units[f"{column}_source"] = (
            units[column].gt(0).map({True: DATA_SOURCE_LABELS["CEMS"]}).astype("object")
        )
    for mass, code_column in UNIT_MEASUREMENT_CODES.items():
        units = units.merge(
            list_measurement_codes(cems, keys, mass, code_column), on=keys, how="left"
        )
    units = (
        units.merge(load_plant_names(path_prefix), on="plant_id_eia", how="left")
        .merge(load_campd_unit_attributes(year), on=keys, how="left")
        .merge(count_unit_generators(path_prefix, year), on=keys, how="left")
        .merge(load_unit_stack_heights(year), on=keys, how="left", validate="1:1")
    )
    table = build_sheet_rows(
        units, UNIT_DATA_SOURCES, year, "SEQUNT", ["state", "plant_name_eia", *keys]
    )
    # every unit in OGE's CEMS data reports to EPA's Clean Air Markets Division
    table["CAPDFLAG"] = "Yes"
    return convert_mass_to_tons(table, "UNT", year)


def list_measurement_codes(
    cems: pd.DataFrame, keys: list[str], mass: str, code_column: str
) -> pd.DataFrame:
    """Lists the CEMS measurement codes of each unit's emissions of a pollutant, from
    the code covering the most emissions to the least, e.g. "Measured; Substitute".

    Args:
        cems (pd.DataFrame): hourly CEMS data with the keys, mass and code columns.
        keys (list[str]): columns identifying each unit.
        mass (str): emissions column.
        code_column (str): column of the measurement codes of the emissions.

    Returns:
        pd.DataFrame: the keys and the codes, in a column named after the mass column
            with "_source". Units without positive emissions are left out.
    """
    return (
        cems.loc[cems[mass] > 0]
        .groupby(keys + [code_column], observed=True)[mass]
        .sum()
        .reset_index()
        .sort_values(mass, ascending=False)
        .groupby(keys)[code_column]
        .agg(lambda codes: "; ".join(codes.astype(str)))
        .rename(f"{mass}_source")
        .reset_index()
    )


def count_unit_generators(path_prefix: str, year: int) -> pd.DataFrame:
    """Counts the generators of each CEMS unit in OGE's subplant crosswalk.

    Args:
        path_prefix (str): prefix of the outputs folder of the year, e.g. "2024/".
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia, emissions_unit_id_epa and generator_count.
    """
    keys = ["plant_id_eia", "emissions_unit_id_epa"]
    crosswalk = pd.read_csv(
        outputs_folder(f"{path_prefix}subplant_crosswalk_{year}.csv.zip"),
        usecols=keys + ["generator_id"],
        dtype={"emissions_unit_id_epa": "str", "generator_id": "str"},
    ).dropna()
    return (
        crosswalk.groupby(keys)["generator_id"]
        .nunique()
        .rename("generator_count")
        .reset_index()
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


def load_unit_stack_heights(year: int) -> pd.DataFrame:
    """Finds the stack height of each CEMS unit from the stacks and flues of its EIA-860
    boilers.

    Only operating and standby stacks are used, leaving out retired, cancelled, planned
    and out-of-service stacks, as eGRID does. Each unit's stack height comes from:

    1. the tallest stack of the boiler whose ID is the unit's ID, as in eGRID.
    2. the stacks of the boilers of the unit's generators, from OGE's subplant
       crosswalk and the EIA-860 boiler-generator associations, if they all have the
       same height. This fills units whose IDs are not boiler IDs, such as combustion
       turbines exhausting through a heat recovery steam generator (HRSG). For combined
       cycles, these associations link each turbine to every boiler of the block,
       not only its own HRSG, so units whose boilers have stacks of different heights
       are left blank rather than given another turbine's stack.

    Stack heights are blank before 2013, when the raw EIA-860 does not report them in
    this layout.

    Args:
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia, emissions_unit_id_epa and stack_height_ft of each
            unit that has an operating stack with a reported height.
    """
    keys = ["plant_id_eia", "emissions_unit_id_epa"]
    stacks = load_raw_eia860_stack_heights(year)
    stacks = stacks.loc[
        stacks["operational_status_code"].isin(OPERATING_STACK_STATUS_CODES),
        ["plant_id_eia", "stack_flue_id_eia", "stack_height_ft"],
    ]
    boiler_heights = (
        load_pudl_table(
            "core_eia860__assn_boiler_stack_flue",
            year,
            columns=["plant_id_eia", "boiler_id", "stack_flue_id_eia"],
        )
        .merge(
            stacks,
            on=["plant_id_eia", "stack_flue_id_eia"],
            how="inner",
            validate="m:1",
        )
        .groupby(["plant_id_eia", "boiler_id"], dropna=False)["stack_height_ft"]
        .max()
        .reset_index()
    )
    # 1. the boiler with the unit's ID
    heights_by_id = boiler_heights.rename(
        columns={"boiler_id": "emissions_unit_id_epa"}
    )
    # 2. the boilers of the unit's generators, only where their stacks all have the same
    # height, since a combined cycle turbine is linked to every boiler of its block
    heights_by_generator = (
        load_unit_to_boiler_associations(year)
        .dropna(subset=keys + ["boiler_id"])
        .drop_duplicates()
        .merge(boiler_heights, on=["plant_id_eia", "boiler_id"], how="inner")
        .groupby(keys, dropna=False)["stack_height_ft"]
        .agg(["min", "max"])
        .reset_index()
    )
    heights_by_generator = heights_by_generator.loc[
        heights_by_generator["min"] == heights_by_generator["max"], keys + ["max"]
    ].rename(columns={"max": "stack_height_ft"})
    heights = heights_by_id.merge(
        heights_by_generator,
        on=keys,
        how="outer",
        suffixes=("", "_from_generators"),
        validate="1:1",
    )
    heights["stack_height_ft"] = heights["stack_height_ft"].fillna(
        heights["stack_height_ft_from_generators"]
    )
    return heights[keys + ["stack_height_ft"]]


########################################################################################
# Generator sheet
########################################################################################


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
    # months before a generator started operating or after it retired have no
    # generation, so their blanks are zeros in the ozone season totals, which are
    # otherwise blank when unknown
    seasonal_data = fill_zeros_when_not_operating(
        allocated.merge(
            dates[keys + ["generator_operating_date", "generator_retirement_date"]],
            on=keys,
            how="left",
            validate="m:1",
        ),
        ["net_generation_mwh"],
    )
    generators = (
        allocated.groupby(keys)
        .agg(
            prime_mover_code=("prime_mover_code", "first"),
            net_generation_mwh=(
                "net_generation_mwh",
                lambda values: values.sum(min_count=1),
            ),
        )
        .reset_index()
        .merge(
            sum_ozone_season_values(seasonal_data, keys, ["net_generation_mwh"]),
            on=keys,
            how="left",
            validate="1:1",
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
    table = build_sheet_rows(
        generators,
        GENERATOR_DATA_SOURCES,
        year,
        "SEQGEN",
        ["state", "plant_name_eia", *keys],
    )
    table["CFACT"] = calculate_capacity_factor(
        table["GENNTAN"], table["NAMEPCAP"], year
    )
    return table


########################################################################################
# Plant sheet
########################################################################################


def build_plant_table(year: int, path_prefix: str) -> pd.DataFrame:
    """Builds the plant sheet of the eGRID data file from OGE's annual plant data, the
    generation by fuel, combustion heat input, ozone season totals and data sources of
    OGE's monthly subplant data, and the biomass and CHP adjustments between OGE's
    unadjusted and adjusted values.

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
        usecols=[
            "plant_id_eia",
            "latitude",
            "longitude",
            "plant_primary_fuel",
            "fuel_category",
            "nerc_region",
        ],
    ).rename(columns={"fuel_category": "plant_fuel_category"})
    plants = (
        plant_data.merge(plant_attributes, on="plant_id_eia", how="left")
        .merge(load_ba_reference()[["ba_code", "ba_name"]], on="ba_code", how="left")
        .merge(
            load_plant_totals_from_subplants(path_prefix, year),
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
    plants["iso_rto"] = plants["ba_code"].map(BA_ISO_RTO)
    plants["fuel_category_egrid"] = plants["plant_fuel_category"].map(
        PLANT_FUEL_CATEGORIES
    )
    plants = plants.merge(
        load_coal_flags(path_prefix, year), on="plant_id_eia", how="left"
    )
    plants = add_adjustment_values(plants)
    table = build_sheet_rows(
        plants,
        PLANT_DATA_SOURCES,
        year,
        "SEQPLT",
        ["state", "plant_name_eia", "plant_id_eia"],
    )
    for source_code, value_code in PLANT_SOURCE_VALUES.items():
        table[source_code] = table[source_code].where(table[value_code].notna())
    table = convert_mass_to_tons(table, "PLNT", year)
    table = add_emission_rates_and_resource_mix(table, "PLNT", year)
    # the regional nonbaseload values are summed from this column, so that each
    # region's nonbaseload generation matches the nonbaseload factors of its plants
    table["PLNGENNB"] = table["PLNGENAN"] * table["NBFACTOR"]
    # combustion heat input per kWh of combustion generation, blank where either is
    # zero or negative
    table["PLHTRT"] = (
        table["PLHTIAN"].where(table["PLHTIAN"] > 0)
        * 1000
        / table["PLGENACY"].where(table["PLGENACY"] > 0)
    )
    return table


def load_plant_totals_from_subplants(path_prefix: str, year: int) -> pd.DataFrame:
    """Sums OGE's monthly subplant data to each plant's generation by eGRID fuel group,
    combustion heat input and ozone season totals, and lists the data sources of each
    plant's heat input and emissions.

    Each subplant's data counts toward the fuel groups in `EGRID_FUEL_GROUPS` that
    include its fuel category. Plants without subplants in a fuel group have zero
    generation or heat input in that group.

    Args:
        path_prefix (str): prefix of the outputs and results folders of the year, e.g.
            "2024/".
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia, the net generation of each fuel group of
            `GENERATION_BY_FUEL_CODES` (e.g. "net_generation_mwh_coal"), the
            combustion heat input (e.g. "combustion_fuel_consumed_mmbtu"), the ozone
            season totals of heat input, net generation, NOx and combustion heat input
            (with "_ozone_season" added to their column names), the data sources of
            the heat input and emissions (e.g. "nox_mass_lb_source"), and "Yes" in
            pumped_storage_flag for plants with a pumped storage subplant.
    """
    heat_input_columns = ["fuel_consumed_mmbtu", "fuel_consumed_for_electricity_mmbtu"]
    ozone_season_columns = heat_input_columns + [
        "net_generation_mwh",
        "nox_mass_lb",
        "nox_mass_lb_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
    ]
    source_columns = [
        f"{pollutant}_mass_lb"
        for pollutant in ["so2", "co2", "ch4", "n2o", "co2e", "hg"]
    ]
    keys = ["plant_id_eia", "subplant_id", "report_date"]
    subplants = pd.read_csv(
        results_folder(f"{path_prefix}plant_data/monthly/us_units/subplant_data.csv"),
        usecols=keys + ozone_season_columns + source_columns,
        parse_dates=["report_date"],
    ).merge(
        pd.read_csv(
            results_folder(f"{path_prefix}plant_data/subplant_attributes.csv"),
            usecols=[
                "plant_id_eia",
                "subplant_id",
                "fuel_category",
                "subplant_primary_prime_mover_code",
            ],
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
    seasonal_columns = (
        ozone_season_columns
        + [f"combustion_{column}" for column in heat_input_columns]
        + ["net_generation_mwh_combustion"]
    )
    totals = (
        subplants.groupby("plant_id_eia", dropna=False)[
            generation_columns
            + [f"combustion_{column}" for column in heat_input_columns]
        ]
        .sum(min_count=1)
        .reset_index()
        .merge(
            sum_ozone_season_values(
                fill_zeros_when_not_operating(
                    subplants.merge(
                        load_subplant_operating_dates(path_prefix, year),
                        on=["plant_id_eia", "subplant_id"],
                        how="left",
                        validate="m:1",
                    ),
                    seasonal_columns,
                ),
                ["plant_id_eia"],
                seasonal_columns,
            ),
            on="plant_id_eia",
            how="left",
            validate="1:1",
        )
    )
    # the data source of each subplant-month, as eGRID's labels
    subplants = subplants.merge(
        pd.read_csv(
            results_folder(f"{path_prefix}plant_data/plant_metadata.csv"),
            usecols=keys + ["data_source"],
            parse_dates=["report_date"],
        ),
        on=keys,
        how="left",
        validate="1:1",
    )
    subplants["data_source"] = subplants["data_source"].map(DATA_SOURCE_LABELS)
    # flag plants with a pumped storage subplant
    flags = (
        subplants.assign(
            pumped_storage_flag=subplants["subplant_primary_prime_mover_code"]
            == PUMPED_STORAGE_PRIME_MOVER,
        )
        .groupby("plant_id_eia", dropna=False)[["pumped_storage_flag"]]
        .any()
    )
    totals = totals.merge(
        flags.replace({True: "Yes", False: None}).reset_index(),
        on="plant_id_eia",
        how="left",
    )
    for column in ["nox_mass_lb", "fuel_consumed_mmbtu"] + source_columns:
        totals = totals.merge(
            list_data_sources(subplants, column), on="plant_id_eia", how="left"
        )
    ozone_season = subplants["report_date"].dt.month.isin(OZONE_SEASON_MONTHS)
    for column in ["nox_mass_lb", "fuel_consumed_mmbtu"]:
        totals = totals.merge(
            list_data_sources(subplants[ozone_season], column, "_ozone_season"),
            on="plant_id_eia",
            how="left",
        )
    return totals


def load_subplant_operating_dates(path_prefix: str, year: int) -> pd.DataFrame:
    """Loads the dates that each subplant started operating and retired, from the
    operating and retirement dates of its generators.

    A subplant started operating when its first generator did, and retired when its
    last generator did. Subplants with a generator that has not retired have no
    retirement date.

    Args:
        path_prefix (str): prefix of the outputs folder of the year, e.g. "2024/".
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia, subplant_id, generator_operating_date and
            generator_retirement_date of each subplant.
    """
    generators = pd.read_csv(
        outputs_folder(f"{path_prefix}subplant_crosswalk_{year}.csv.zip"),
        usecols=[
            "plant_id_eia",
            "subplant_id",
            "generator_id",
            "generator_operating_date",
            "generator_retirement_date",
        ],
        dtype={"generator_id": "str"},
        parse_dates=["generator_operating_date", "generator_retirement_date"],
    ).drop_duplicates(subset=["plant_id_eia", "generator_id"])
    return (
        generators.groupby(["plant_id_eia", "subplant_id"], dropna=False)
        .agg(
            generator_operating_date=("generator_operating_date", "min"),
            # a generator without a retirement date has not retired, so neither has
            # its subplant
            generator_retirement_date=(
                "generator_retirement_date",
                lambda dates: dates.max() if dates.notna().all() else pd.NaT,
            ),
        )
        .reset_index()
    )


def load_coal_flags(path_prefix: str, year: int) -> pd.DataFrame:
    """Flags the plants that burned coal or generated electricity from coal in OGE's
    EIA-923 data by subplant and energy source, including plants whose coal is not
    their subplants' primary fuel.

    Args:
        path_prefix (str): prefix of the outputs folder of the year, e.g. "2024/".
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia and coal_flag ("Yes") of each plant with coal.
    """
    eia923 = assign_fuel_category_to_esc(
        pd.read_csv(
            outputs_folder(f"{path_prefix}subplant_eia923_{year}.csv.zip"),
            usecols=[
                "plant_id_eia",
                "energy_source_code",
                "fuel_consumed_mmbtu",
                "net_generation_mwh",
            ],
        ),
        fuel_category_names=["fuel_category"],
    )
    coal = eia923[
        eia923["fuel_category"].isin(EGRID_FUEL_GROUPS["coal"])
        & ((eia923["fuel_consumed_mmbtu"] > 0) | (eia923["net_generation_mwh"] > 0))
    ]
    return pd.DataFrame(
        {"plant_id_eia": coal["plant_id_eia"].unique(), "coal_flag": "Yes"}
    )


def list_data_sources(
    subplants: pd.DataFrame, column: str, suffix: str = ""
) -> pd.DataFrame:
    """Lists the data sources of each plant's value of a column, i.e. the sources of
    the subplant-months where the value is positive, e.g. "EPA/CAPD; EIA".

    Args:
        subplants (pd.DataFrame): monthly subplant data with plant_id_eia, the column
            and a "data_source" column of eGRID labels.
        column (str): column whose sources are listed.
        suffix (str, optional): added to the column name of the result, e.g.
            "_ozone_season". Defaults to "".

    Returns:
        pd.DataFrame: plant_id_eia and the sources, in a column named after the column
            with the suffix and "_source". Plants without a positive value are left out.
    """
    order = list(DATA_SOURCE_LABELS.values())
    return (
        subplants.loc[subplants[column] > 0]
        .groupby("plant_id_eia")["data_source"]
        .agg(
            lambda sources: "; ".join(
                sorted(sources.dropna().unique(), key=order.index)
            )
        )
        .rename(f"{column}{suffix}_source")
        .reset_index()
    )


def load_eia860_plant_attributes(year: int) -> pd.DataFrame:
    """Loads the attributes of each plant reported in EIA-860.

    Args:
        year (int): data year.

    Returns:
        pd.DataFrame: plant_id_eia, name, state, county, location, state and county
            FIPS codes, transmission or distribution system owner, utility, eGRID
            sector of each plant. The county FIPS code is the three digits within the
            state, as in eGRID.
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


def add_adjustment_values(plants: pd.DataFrame) -> pd.DataFrame:
    """Adds the biomass and CHP adjustment values and flags of each plant, from OGE's
    unadjusted and adjusted values.

    The biomass adjustment is the unadjusted value minus the biomass-adjusted value,
    and the CHP adjustment is the biomass-adjusted value minus the biomass- and
    CHP-adjusted value, so that each plant's adjusted value is its unadjusted value
    minus both adjustments. A plant is flagged where OGE adjusted it.

    Args:
        plants (pd.DataFrame): annual plant data with OGE's unadjusted and adjusted
            columns, and the ozone season and combustion totals of
            `load_plant_totals_from_subplants()`.

    Returns:
        pd.DataFrame: the plant data with the adjustment and flag columns of
            `PLANT_DATA_SOURCES`.
    """
    plants = plants.copy()
    for pollutant in ADJUSTED_POLLUTANTS:
        mass = f"{pollutant}_mass_lb"
        plants[f"{mass}_biomass_adjustment"] = plants[mass] - plants[f"{mass}_adjusted"]
        plants[f"{mass}_chp_adjustment"] = (
            plants[f"{mass}_adjusted"] - plants[f"{mass}_for_electricity_adjusted"]
        )
    plants["nox_mass_lb_biomass_adjustment_ozone_season"] = (
        plants["nox_mass_lb_ozone_season"] - plants["nox_mass_lb_adjusted_ozone_season"]
    )
    plants["nox_mass_lb_chp_adjustment_ozone_season"] = (
        plants["nox_mass_lb_adjusted_ozone_season"]
        - plants["nox_mass_lb_for_electricity_adjusted_ozone_season"]
    )
    for season in ["", "_ozone_season"]:
        plants[f"combustion_fuel_consumed_mmbtu_chp_adjustment{season}"] = (
            plants[f"combustion_fuel_consumed_mmbtu{season}"]
            - plants[f"combustion_fuel_consumed_for_electricity_mmbtu{season}"]
        )
    biomass_adjusted = (
        plants[[f"{p}_mass_lb_biomass_adjustment" for p in ADJUSTED_POLLUTANTS]] > 0
    ).any(axis=1)
    chp_adjusted = (
        plants["fuel_consumed_mmbtu"] > plants["fuel_consumed_for_electricity_mmbtu"]
    )
    plants["biomass_adjustment_flag"] = biomass_adjusted.map({True: "Yes"})
    plants["chp_adjustment_flag"] = chp_adjusted.map({True: "Yes"})
    return plants


########################################################################################
# State, BA and U.S. sheets
########################################################################################


def build_plant_region_table(year: int, path_prefix: str, sheet: str) -> pd.DataFrame:
    """Builds the state or NERC region sheet of the eGRID data file by summing OGE's
    monthly subplant results to each plant's state or NERC region in OGE's plant static
    attributes, including the emission rates by fuel type.

    Each region's data are shaped like OGE's power sector results, i.e. one row for
    each subplant fuel category and a total row, so that a region's columns are
    calculated the same way as a BA's. Generation and heat input from fuel groups that
    a region does not have are zero. Plants without a region are left out, and if no
    plant has a region, the table has no rows.

    Args:
        year (int): data year.
        path_prefix (str): prefix of the results folder of the year, e.g. "2024/".
        sheet (str): "ST" or "NRL".

    Returns:
        pd.DataFrame: one row for each region with plants in OGE's results, in U.S.
            units, with the year, the state abbreviation and FIPS code or the NERC
            region acronym and name, and the columns of `REGION_DATA_SOURCES`.
    """
    attribute, code = PLANT_REGION_ATTRIBUTES[sheet]
    data_columns = list(
        dict.fromkeys(
            [column for _, column, _ in REGION_DATA_SOURCES.values()]
            + list(FUEL_TYPE_RATE_MASSES.values())
        )
    )
    subplants = (
        pd.read_csv(
            results_folder(
                f"{path_prefix}plant_data/monthly/us_units/subplant_data.csv"
            ),
            usecols=["plant_id_eia", "subplant_id", "report_date"] + data_columns,
            parse_dates=["report_date"],
        )
        .merge(
            pd.read_csv(
                results_folder(f"{path_prefix}plant_data/subplant_attributes.csv"),
                usecols=["plant_id_eia", "subplant_id", "fuel_category"],
            ),
            on=["plant_id_eia", "subplant_id"],
            how="left",
            validate="m:1",
        )
        .merge(
            pd.read_csv(
                results_folder(f"{path_prefix}plant_data/plant_static_attributes.csv"),
                usecols=["plant_id_eia", attribute],
            ).rename(columns={attribute: "region"}),
            on="plant_id_eia",
            how="left",
            validate="m:1",
        )
    )
    # plants without a region are left out of the sheet
    subplants = subplants[subplants["region"].notna()]
    if subplants.empty:
        logger.warning(f"No plant has a {attribute}, so the {sheet} sheet is empty")
        return pd.DataFrame(columns=[c for c, _, _ in get_egrid_columns(sheet, year)])
    # monthly data of each region by fuel category, with a total row for each month
    monthly = pd.concat(
        [
            subplants.groupby(["region", "fuel_category", "report_date"], dropna=False)[
                data_columns
            ]
            .sum(min_count=1)
            .reset_index(),
            subplants.groupby(["region", "report_date"], dropna=False)[data_columns]
            .sum(min_count=1)
            .reset_index()
            .assign(fuel_category="total"),
        ]
    )
    annual = (
        monthly.groupby(["region", "fuel_category"], dropna=False)[data_columns]
        .sum(min_count=1)
        .reset_index()
    )
    prefix = REGION_SHEETS[sheet][0]
    rows = [
        {"YEAR": year, code: region}
        | calculate_region_values(
            prefix,
            annual[annual["region"] == region],
            monthly[monthly["region"] == region],
        )
        for region in sorted(subplants["region"].unique())
    ]
    table = pd.DataFrame(rows)
    if sheet == "ST":
        fips = load_pudl_table(
            "core_pudl__codes_subdivisions",
            columns=["country_code", "subdivision_code", "state_id_fips"],
        )
        table = table.merge(
            fips.loc[
                fips["country_code"] == "USA", ["subdivision_code", "state_id_fips"]
            ].rename(
                columns={"subdivision_code": "PSTATABB", "state_id_fips": "FIPSST"}
            ),
            on="PSTATABB",
            how="left",
            validate="1:1",
        )
    else:
        table["NERCNAME"] = table["NERC"].map(NERC_REGION_NAMES)
        if table["NERCNAME"].isna().any():
            logger.warning(
                "The following NERC regions have no name: "
                f"{table.loc[table['NERCNAME'].isna(), 'NERC'].tolist()}"
            )
    table = convert_mass_to_tons(table, sheet, year)
    return add_emission_rates_and_resource_mix(table, sheet, year)


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
        row = {"YEAR": year} | calculate_region_values(
            prefix,
            pd.read_csv(results_folder(folder.format(resolution="annual"))),
            pd.read_csv(
                results_folder(folder.format(resolution="monthly")),
                parse_dates=["report_date"],
            ),
        )
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


def calculate_region_values(
    prefix: str, annual: pd.DataFrame, monthly: pd.DataFrame
) -> dict[str, float]:
    """Calculates the columns of `REGION_DATA_SOURCES` and the emission rates by fuel
    type of one region from its power sector results, in lb.

    Args:
        prefix (str): prefix of the sheet's codes, e.g. "BA".
        annual (pd.DataFrame): annual power sector results of the region, with one row
            for each fuel category and a "total" row.
        monthly (pd.DataFrame): monthly power sector results of the region, with a
            report_date column.

    Returns:
        dict[str, float]: each value, keyed by its code with the sheet's prefix.
    """
    region_data = {
        "annual": annual.set_index("fuel_category"),
        "ozone": monthly[monthly["report_date"].dt.month.isin(OZONE_SEASON_MONTHS)]
        .groupby("fuel_category")
        .sum(numeric_only=True, min_count=1),
    }
    values = {}
    for code, (fuel_category, column, season) in REGION_DATA_SOURCES.items():
        fuel_categories = EGRID_FUEL_GROUPS.get(fuel_category, [fuel_category])
        data = region_data[season]
        if data.index.isin(fuel_categories).any():
            values[prefix + code] = data.loc[
                data.index.isin(fuel_categories), column
            ].sum(min_count=1)
        elif fuel_category != "total":
            # the region has none of the fuel categories
            values[prefix + code] = 0
    return values | {
        prefix + code: rate
        for code, rate in calculate_fuel_type_rates(
            region_data["annual"], region_data["ozone"]
        ).items()
    }


########################################################################################
# Values from more than one sheet
########################################################################################


def add_regional_capacity(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Fills the nameplate capacity of the state, BA, NERC region and U.S. sheets with
    the sum of the nameplate capacity of their plants on the plant sheet.

    Args:
        tables (dict[str, pd.DataFrame]): tables of the sheets, including PLNT, ST, BA,
            NRL and US.

    Returns:
        dict[str, pd.DataFrame]: the tables, with the regional sheets' capacity.
    """
    tables = tables.copy()
    plants = tables["PLNT"]
    for sheet, key in [("ST", "PSTATABB"), ("BA", "BACODE"), ("NRL", "NERC")]:
        capacity = plants.groupby(key)["NAMEPCAP"].sum(min_count=1)
        tables[sheet] = tables[sheet].assign(
            **{f"{REGION_SHEETS[sheet][0]}NAMEPCAP": tables[sheet][key].map(capacity)}
        )
    tables["US"] = tables["US"].assign(USNAMEPCAP=plants["NAMEPCAP"].sum(min_count=1))
    return tables


def add_regional_nonbaseload_values(
    tables: dict[str, pd.DataFrame], year: int
) -> dict[str, pd.DataFrame]:
    """Fills the nonbaseload generation, output emission rates, generation by fuel type
    and resource mix of the state, BA, NERC region and U.S. sheets from the plants on
    the plant sheet, as in eGRID.

    Each plant's nonbaseload generation and emissions are its annual and ozone season
    values times its nonbaseload factor. A region's nonbaseload rates divide the sum of
    its plants' nonbaseload emissions by the sum of their nonbaseload generation, and
    its nonbaseload generation is the sum of the plant sheet's PLNGENNB, so the two are
    always consistent. Each plant's nonbaseload generation counts toward the fuel type of
    its plant fuel category. Regions without any plant with a nonbaseload factor have
    blank nonbaseload values, and the Hg rate is blank.

    Args:
        tables (dict[str, pd.DataFrame]): tables of the sheets, including PLNT, ST, BA,
            NRL and US, with the plant sheet's masses in the units of its columns.
        year (int): data year.

    Returns:
        dict[str, pd.DataFrame]: the tables, with the regional sheets' nonbaseload
            values.
    """
    tables = tables.copy()
    plants = convert_mass_to_lb(tables["PLNT"], "PLNT", year)
    factor = plants["NBFACTOR"]
    # nonbaseload values of each plant, keyed by their codes without the prefix
    nonbaseload = pd.DataFrame(
        {
            "NGENNB": plants["PLNGENNB"],
            "NGENOZ": plants["PLNGENOZ"] * factor,
            "NOXOZ": plants["PLNOXOZ"] * factor,
        }
        | {
            mass_code: plants[f"PL{mass_code}"] * factor
            for mass_code in EMISSION_RATE_MASSES.values()
        }
        | {
            generation_code: plants["PLNGENNB"].where(plants["PLFUELCT"] == category)
            for category, generation_code in NONBASELOAD_GENERATION_CODES.items()
        }
    )
    for sheet, key in [("ST", "PSTATABB"), ("BA", "BACODE"), ("NRL", "NERC")]:
        # plants without a state, BA or NERC region are only counted in the U.S.
        values = calculate_nonbaseload_values(
            nonbaseload.groupby(plants[key], dropna=True).sum(min_count=1)
        )
        tables[sheet] = tables[sheet].assign(
            **{
                f"{REGION_SHEETS[sheet][0]}{code}": tables[sheet][key].map(values[code])
                for code in values.columns
            }
        )
    values = calculate_nonbaseload_values(nonbaseload.sum(min_count=1).to_frame().T)
    tables["US"] = tables["US"].assign(
        **{f"US{code}": values[code].iloc[0] for code in values.columns}
    )
    return tables


def calculate_nonbaseload_values(sums: pd.DataFrame) -> pd.DataFrame:
    """Calculates the nonbaseload values of regions from the sums of their plants'
    nonbaseload generation and emissions.

    Args:
        sums (pd.DataFrame): sums of each region's nonbaseload generation ("NGENNB"),
            ozone season generation ("NGENOZ"), emissions in lb (the masses of
            `EMISSION_RATE_MASSES` and "NOXOZ") and generation by fuel type (the
            columns of `NONBASELOAD_GENERATION_CODES`), indexed by region.

    Returns:
        pd.DataFrame: each region's nonbaseload generation, output emission rates in
            lb/MWh, generation by fuel type and resource mix, with their codes without
            the sheet's prefix as columns.
    """
    total = sums["NGENNB"]
    values = pd.DataFrame({"NGENNB": total}, index=sums.index)
    for pollutant, mass_code in EMISSION_RATE_MASSES.items():
        values[f"{NONBASELOAD_RATE_PREFIX}{pollutant}"] = sums[mass_code] / total.where(
            total > 0
        )
    values[NONBASELOAD_OZONE_NOX_RATE] = sums["NOXOZ"] / sums["NGENOZ"].where(
        sums["NGENOZ"] > 0
    )
    # regions with nonbaseload generation have zero generation in the fuel types none
    # of their nonbaseload plants are in
    for generation_code in NONBASELOAD_GENERATION_CODES.values():
        values[generation_code] = sums[generation_code].where(
            total.isna(), sums[generation_code].fillna(0)
        )
    for mix_code, generation_code in NONBASELOAD_RESOURCE_MIX_GENERATION.items():
        values[mix_code] = values[generation_code] / total.where(total > 0)
    return values


def add_plant_unit_values(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Fills the values of the plant sheet that come from the plant's rows on the unit
    and generator sheets: the number of units and generators, and the CAPD flag.

    The unit sheet has CEMS units only, so plants without CEMS units have no units. A
    plant is flagged as reporting to EPA's Clean Air Markets Division if it has units
    on the unit sheet, including plants whose CEMS data only covers some of their
    subplants or months, whose plant-level values then come partly from EIA.

    Args:
        tables (dict[str, pd.DataFrame]): tables of the sheets, including UNT, GEN and
            PLNT.

    Returns:
        dict[str, pd.DataFrame]: the tables, with the plant sheet's counts and CAPD
            flag.
    """
    tables = tables.copy()
    plants = tables["PLNT"]
    tables["PLNT"] = plants.assign(
        NUMUNT=plants["ORISPL"].map(tables["UNT"]["ORISPL"].value_counts()).fillna(0),
        NUMGEN=plants["ORISPL"].map(tables["GEN"]["ORISPL"].value_counts()).fillna(0),
        CAPDFLAG=plants["ORISPL"]
        .isin(tables["UNT"]["ORISPL"])
        .map({True: "Yes", False: None}),
    )
    return tables


########################################################################################
# Calculations shared by the sheets
########################################################################################


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


def sum_ozone_season_values(
    df: pd.DataFrame, keys: list[str], columns: list[str]
) -> pd.DataFrame:
    """Sums the ozone season (May to September) values of each group, keeping blank
    the values that are unknown.

    A group's ozone season value is blank when all of its ozone season values are
    blank, so a missing measurement is not published as zero. A group with no ozone
    season records at all did not operate in the ozone season, so its value is zero
    wherever its annual value is known.

    Args:
        df (pd.DataFrame): data with the keys, the columns and report_date.
        keys (list[str]): columns identifying each group, e.g. ["plant_id_eia"].
        columns (list[str]): columns to sum.

    Returns:
        pd.DataFrame: the keys of each group in the data and its ozone season totals,
            with "_ozone_season" added to the column names.
    """
    annual = df.groupby(keys, dropna=False)[columns].sum(min_count=1)
    in_season = df[df["report_date"].dt.month.isin(OZONE_SEASON_MONTHS)]
    seasonal = (
        in_season.groupby(keys, dropna=False)[columns]
        .sum(min_count=1)
        .reindex(annual.index)
    )
    # groups without ozone season records did not operate in the ozone season, so
    # their value is zero, unless their annual value is unknown too
    has_season_records = (
        in_season.groupby(keys, dropna=False)
        .size()
        .reindex(annual.index, fill_value=0)
        .gt(0)
    )
    no_season = annual.loc[~has_season_records]
    seasonal.loc[~has_season_records] = no_season.where(no_season.isna(), 0)
    return seasonal.add_suffix("_ozone_season").reset_index()


def fill_zeros_when_not_operating(
    data: pd.DataFrame, columns: list[str]
) -> pd.DataFrame:
    """Fills blank values with zeros in months before a generator or subplant started
    operating or after it retired, since it had no generation, heat input or
    emissions in those months.

    Args:
        data (pd.DataFrame): monthly data with report_date, the columns,
            generator_operating_date and generator_retirement_date.
        columns (list[str]): columns to fill.

    Returns:
        pd.DataFrame: a copy of the data, with zeros in the blank values of the
            columns in months outside the operating period.
    """
    data = data.copy()
    month = data["report_date"].dt.to_period("M")
    not_operating = (month < data["generator_operating_date"].dt.to_period("M")) | (
        month > data["generator_retirement_date"].dt.to_period("M")
    )
    data.loc[not_operating, columns] = data.loc[not_operating, columns].fillna(0)
    return data


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


def add_emission_rates_and_resource_mix(
    table: pd.DataFrame, sheet: str, year: int
) -> pd.DataFrame:
    """Adds the emission rates and resource mix to the plant sheet or a regional sheet, from
    the sheet's emissions, heat input and generation.

    Rates are in lb per MWh or MMBtu, and are blank where the denominator is zero or
    negative, as in eGRID. Each fuel's share of the resource mix is its positive
    generation divided by the sum of the positive generation of every fuel, so that
    the shares of each row add up to 100%. The share of a group of fuels, such as
    renewables, is the sum of the positive generation of its fuels on the same basis.

    Args:
        table (pd.DataFrame): table of the sheet with masses in the units of eGRID's
            columns.
        sheet (str): "PLNT", "ST", "BA", "NRL" or "US".
        year (int): data year.

    Returns:
        pd.DataFrame: the table with the emission rate and resource mix columns.
    """
    table = table.copy()
    prefix = "PL" if sheet == "PLNT" else REGION_SHEETS[sheet][0]
    masses_lb = convert_mass_to_lb(table, sheet, year)

    def rate(mass: pd.Series, denominator: pd.Series) -> pd.Series:
        return mass / denominator.where(denominator > 0)

    for rate_code, (
        annual,
        ozone_season,
        ozone_rate_code,
    ) in EMISSION_RATE_DENOMINATORS.items():
        for pollutant, mass_code in EMISSION_RATE_MASSES.items():
            table[f"{prefix}{pollutant}{rate_code}"] = rate(
                masses_lb[prefix + mass_code], table[prefix + annual]
            )
        table[f"{prefix}NOX{ozone_rate_code}"] = rate(
            masses_lb[f"{prefix}NOXOZ"], table[prefix + ozone_season]
        )

    total = (
        table[[prefix + code for code in RESOURCE_MIX_FUELS]]
        .clip(lower=0)
        .sum(axis=1, min_count=1)
    )
    for mix_code, generation_code in RESOURCE_MIX_GENERATION.items():
        fuels = RESOURCE_MIX_GROUP_FUELS.get(generation_code, [generation_code])
        positive_generation = (
            table[[prefix + fuel for fuel in fuels]]
            .clip(lower=0)
            .sum(axis=1, min_count=1)
        )
        table[prefix + mix_code] = positive_generation / total.where(total > 0)
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


def convert_mass_to_lb(table: pd.DataFrame, sheet: str, year: int) -> pd.DataFrame:
    """Converts the masses of an eGRID sheet's columns whose headers are in short tons
    back to lb, e.g. to calculate emission rates in lb/MWh.

    Args:
        table (pd.DataFrame): table of an eGRID sheet, with masses in the units of the
            sheet.
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        year (int): data year.

    Returns:
        pd.DataFrame: the table with masses in lb.
    """
    table = table.copy()
    for code, header, _ in get_egrid_columns(sheet, year):
        if code in table and get_unit(header) == "tons":
            table[code] = table[code] * ConversionFactors.short_ton_to_lbs
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
