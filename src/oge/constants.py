# specify the energy_source_codes that are considered clean/carbon-free
CLEAN_FUELS = ["SUN", "MWH", "WND", "WAT", "WH", "PUR", "NUC"]

# specify the energy_source_codes that are considerd to be biomass
BIOMASS_FUELS = [
    "AB",
    "BG",
    "BLQ",
    "DG",
    "LFG",
    "MSB",
    "OBG",
    "OBL",
    "OBS",
    "SLW",
    "WDL",
    "WDS",
]

TIME_RESOLUTIONS = {"hourly": "H", "monthly": "M", "annual": "A"}

# derived from table 2.4-4 of the EPA's AP-42 document
nox_lb_per_mmbtu_flared_landfill_gas = 0.078

# values assumed by eGRID for CHP efficiency
chp_gross_thermal_output_efficiency = 0.8
chp_useful_thermal_output_efficiency = 0.75


class ConversionFactors(float):
    """Defines conversion factors between common units."""

    lb_to_kg = 0.453592
    mmbtu_to_GJ = 1.055056
    mwh_to_mmbtu = 3.412142
    short_ton_to_lbs = 2000
