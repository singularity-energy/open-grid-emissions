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
latest_validated_year = 2025
# current_early_release_year is the year for which non-final (early-release) data
# is available from the EIA. This enables running the OGE pipeline for this year
# EIA-860ER data is generally available in June and EIA-923ER data is generally
# available in July of the following year. This should not be updated to the next year
# until ER data is available, so for part of the year, latest_validated_year will equal
# current_early_release_year
# TODO: Change this to 2026 around July 2027 (check PUDL to see when integrated)
current_early_release_year = 2025

# specify the energy_source_codes that are considered clean/carbon-free
CLEAN_FUELS = ["SUN", "MWH", "WND", "WAT", "WH", "PUR", "NUC"]

# specify the fuel categories included in the fossil-only residual mix, which
# approximates the residual mix by excluding resources that may generate RECs.
# Following the CRS guidance for fossil-only residual mixes, this only includes coal,
# oil, and gas. Waste is excluded even though it includes some fossil-derived fuels
# (e.g. non-biogenic municipal solid waste and tire-derived fuels), since it also
# includes biogenic waste and waste-to-energy can generate RECs in some states, so
# excluding it keeps the residual mix conservative
RESIDUAL_MIX_FUEL_CATEGORIES = ["coal", "natural_gas", "petroleum"]

# specify the energy_source_codes that are considerd to be biomass
BIOMASS_FUELS = [
    "AB",  # Agricultural Byproducts
    "BG",
    "BLQ",  # Black Liquor
    "DG",  # Digester Gas
    "LFG",  # Landfill Gas
    "MSB",  # Municipal Solid Waste (Biogenic Component)
    "OBG",  # Other Biomass Gas
    "OBL",  # Other Biomass Liquid
    "OBS",  # Other Biomass Solid
    "SLW",  # Sludge Waste
    "WDL",  # Wood Derived Liquids
    "WDS",  # Wood Derived Solids
]

ENERGY_STORAGE_PRIME_MOVERS = [
    "BA",  # Battery
    "CE",  # Compressed Air Storage
    "ES",  # Energy Storage (Other)
    "FW",  # Flywheel Storage
    "PS",  # Pumped Storage
]

# storage categories (see energy_storage.STORAGE_CATEGORY_BY_METHOD) of plants whose
# solar and wind generation EIA-930 reports as "with integrated battery storage". The
# EIA-930 instructions don't define "integrated", but describe it as storage whose
# charging is netted from the generator's output. Co-located storage is included,
# since in 2025 most BAs that report solar with integrated battery storage have no
# solar at plants with hybrid storage, while their solar at plants with co-located
# storage matches the reported generation.
EIA930_INTEGRATED_STORAGE_CATEGORIES = ["hybrid", "co_located"]

TIME_RESOLUTIONS = {"hourly": "h", "monthly": "M", "annual": "A"}

# values assumed by eGRID for CHP efficiency
chp_gross_thermal_output_efficiency = 0.8
chp_useful_thermal_output_efficiency = 0.75

# the list of status codes we want to filter out
cancelled_status_codes = [
    "CN",  # Cancelled, but previously reported as planned
    "IP",  # Indefinitely postponed or no longer in resource plan
]


class ConversionFactors(float):
    """Defines conversion factors between common units."""

    lb_to_kg = 0.453592
    kgal_to_barrel = 23.8095
    mmbtu_to_GJ = 1.055056
    mwh_to_mmbtu = 3.412142
    short_ton_to_lbs = 2000
