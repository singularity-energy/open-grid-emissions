# specify the date ranges for OGE data
# earliest_data_year is the earliest year for which cleaned, standardized, and stable
# EIA and CEMS data is available through PUDL. PUDL includes data prior to 2005, but
# for now we do not consider pre-2005 data to be clean enough to work with
earliest_data_year = 2005
# earliest_validated_year is the earliest data for which OGE data is currently published
# and validated. The pipeline could be run for data prior to this year, but it has not
# been validated
earliest_validated_year = 2005
# earliest_hourly_data_year is the most recent year for which OGE can produce hourly
# profiles
earliest_hourly_data_year = 2019
# latest_validated_year is the most recent year for which OGE data has been published
latest_validated_year = 2023
# current_early_release_year is the year for which non-final (early-release) data
# is available from the EIA. This enables running the OGE pipeline for this year
# EIA-860ER data is generally available in June and EIA-923ER data is generally
# available in July of the following year. This should not be updated to the next year
# until ER data is available, so for part of the year, latest_validated_year will equal
# current_early_release_year
# TODO: Change this to 2024 around July 2025 (check PUDL to see when integrated)
current_early_release_year = 2023

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
    kgal_to_barrel = 23.8095
    mmbtu_to_GJ = 1.055056
    mwh_to_mmbtu = 3.412142
    short_ton_to_lbs = 2000
