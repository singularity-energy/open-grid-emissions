"""
Check columns for standard data files output by data_pipeline.

Since file names and column names are hardcoded across several files, calling these checks
during file creation (data_pipeline.py) ensures that changes to file names and
column names are not made accidentally.

To make an intentional change in a file or column name, search the project for all
uses of that column/file, update all of them to the new column name, and then change
the name here.

To add a column, add the name here.

To remove a column, search the project for all uses of that column and remove
those files or uses, then remove it here.

After any change, re-run data_pipeline to regenerate all files and re-run these
checks.
"""

from oge.logging_util import get_logger

import pandas as pd

logger = get_logger(__name__)


COLUMNS = {
    "eia923_allocated": {
        "report_date",
        "plant_id_eia",
        "generator_id",
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "co2e_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "co2e_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "co2_mass_lb_for_electricity_adjusted",
        "ch4_mass_lb_for_electricity_adjusted",
        "n2o_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
        "subplant_id",
        "prime_mover_code",
        "energy_source_code",
        "hourly_data_source",
    },
    "cems_cleaned": {
        "plant_id_eia",
        "emissions_unit_id_epa",
        "datetime_utc",
        "operating_time_hours",
        "gross_generation_mwh",
        "steam_load_1000_lb",
        "fuel_consumed_mmbtu",
        "co2_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "plant_id_epa",
        "co2_mass_measurement_code",
        "nox_mass_measurement_code",
        "so2_mass_measurement_code",
        "report_date",
        "subplant_id",
        "energy_source_code",
        "ch4_mass_lb",
        "n2o_mass_lb",
    },
    "cems_subplant": {
        "plant_id_eia",
        "subplant_id",
        "datetime_utc",
        "report_date",
        "gross_generation_mwh",
        "gtn_method",
        "net_generation_mwh",
        "steam_load_1000_lb",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "co2e_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "co2e_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "co2_mass_lb_for_electricity_adjusted",
        "ch4_mass_lb_for_electricity_adjusted",
        "n2o_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
    },
    "partial_cems_subplant": {
        "report_date",
        "plant_id_eia",
        "subplant_id",
        "datetime_utc",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "net_generation_mwh",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "co2e_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "co2e_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "co2_mass_lb_for_electricity_adjusted",
        "ch4_mass_lb_for_electricity_adjusted",
        "n2o_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
    },
    "partial_cems_plant": {
        "report_date",
        "plant_id_eia",
        "subplant_id",
        "datetime_utc",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "net_generation_mwh",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "co2e_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "co2e_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "co2_mass_lb_for_electricity_adjusted",
        "ch4_mass_lb_for_electricity_adjusted",
        "n2o_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
    },
    "plant_static_attributes": {
        "plant_id_eia",
        "plant_primary_fuel",
        "fuel_category",
        "fuel_category_eia930",
        "ba_code",
        "ba_code_physical",
        "distribution_flag",
        "timezone",
        "data_availability",
        "shaped_plant_id",
        "latitude",
        "longitude",
        "state",
        "county",
        "city",
        "plant_name_eia",
        "capacity_mw",
        "plant_operating_date",
        "plant_retirement_date",
    },
    "plant_metadata": {
        "plant_id_eia",
        "subplant_id",
        "report_date",
        "data_source",
        "hourly_profile_source",
        "net_generation_method",
    },
    "hourly_profiles": {
        "ba_code",
        "fuel_category",
        "datetime_utc",
        "datetime_local",
        "report_date",
        "eia930_profile",
        "cems_profile",
        "residual_profile",
        "scaled_residual_profile",
        "shifted_residual_profile",
        "imputed_profile",
        "profile",
        "flat_profile",
        "profile_method",
    },
    "shaped_eia923_data": {
        "plant_id_eia",
        "datetime_utc",
        "report_date",
        "ba_code",
        "fuel_category",
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "co2e_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "co2e_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "co2_mass_lb_for_electricity_adjusted",
        "ch4_mass_lb_for_electricity_adjusted",
        "n2o_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
        "profile_method",
    },
    "annual_generation_averages_by_fuel": {
        "fuel_category",
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "co2e_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "co2e_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "co2_mass_lb_for_electricity_adjusted",
        "ch4_mass_lb_for_electricity_adjusted",
        "n2o_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
        "generated_co2_rate_lb_per_mwh_for_electricity",
        "generated_ch4_rate_lb_per_mwh_for_electricity",
        "generated_n2o_rate_lb_per_mwh_for_electricity",
        "generated_co2e_rate_lb_per_mwh_for_electricity",
        "generated_nox_rate_lb_per_mwh_for_electricity",
        "generated_so2_rate_lb_per_mwh_for_electricity",
        "generated_co2_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_ch4_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_n2o_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_co2e_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_nox_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_so2_rate_lb_per_mwh_for_electricity_adjusted",
    },
    "gross_to_net_conversions": {
        "plant_id_eia",
        "subplant_id",
        "report_date",
        "gross_generation_mwh",
        "minimum_gross_generation_mwh",
        "maximum_gross_generation_mwh",
        "capacity_mw",
        "net_generation_mwh",
        "data_source",
        "hours_in_month_subplant",
        "hours_in_month_plant",
        "annual_subplant_ratio",
        "annual_plant_ratio",
        "annual_fleet_ratio",
        "default_gtn_ratio",
        "monthly_subplant_ratio",
        "monthly_plant_ratio",
        "annual_plant_shift_mw",
        "monthly_subplant_shift_mw",
        "annual_subplant_shift_mw",
        "subplant_regression_ratio",
        "subplant_regression_shift_mw",
        "subplant_regression_rsq_adj",
        "plant_regression_ratio",
        "plant_regression_shift_mw",
        "plant_regression_rsq_adj",
    },
    "subplant_crosswalk": {
        "plant_id_epa",
        "emissions_unit_id_epa",
        "plant_id_eia",
        "generator_id",
        "subplant_id",
        "unit_id_pudl",
        "generator_operating_date",
        "generator_retirement_date",
        "current_planned_generator_operating_date",
        "prime_mover_code",
    },
    "shaped_aggregated_plants": {
        "plant_id_eia",
        "report_date",
        "fuel_category",
        "ba_code",
        "aggregated_plants",
    },
    "primary_fuel_table": {
        "plant_id_eia",
        "generator_id",
        "subplant_id",
        "energy_source_code",
        "plant_primary_fuel_from_fuel_consumed_for_electricity_mmbtu",
        "plant_primary_fuel_from_net_generation_mwh",
        "plant_primary_fuel_from_capacity_mw",
        "plant_primary_fuel_from_mode",
        "plant_primary_fuel",
        "subplant_primary_fuel_from_fuel_consumed_for_electricity_mmbtu",
        "subplant_primary_fuel_from_net_generation_mwh",
        "subplant_primary_fuel_from_capacity_mw",
        "subplant_primary_fuel_from_mode",
        "subplant_primary_fuel",
    },
}


def check_columns(df: pd.DataFrame, file_name: str):
    """Given a file name and a data frame to export, check that its columns are as
    expected.

    Args:
        df (pd.DataFrame): table to check.
        file_name (str): key name of file in `COLUMNS`.

    Raises:
        ValueError: if `file_name` cannot be found in list of files.
        ValueError: if columns are missing in `df`.
    """

    cols = set(list(df.columns))
    # Get expected columns
    if file_name not in COLUMNS:
        raise ValueError(
            f"Could not find file {file_name} in expected file names {COLUMNS.keys()}"
        )
    expected_cols = COLUMNS[file_name]

    # Check for extra columns. Warning not exception
    extras = cols - expected_cols
    if len(extras) > 0:
        logger.warning(
            f"columns {extras} in {file_name} are not guaranteed by column_checks.py"
        )

    # Raise exception for missing columns
    missing = expected_cols - cols
    if len(missing) > 0:
        raise ValueError(f"Columns {missing} missing from {file_name}")

    return


def get_dtypes() -> dict:
    """Returns a dictionary of dtypes that should be used for each column name.

    We use both the dtypes defined in pudl as well as dtypes for columns specific to OGE

    Returns:
        dict: combined dictionary of pudl and oge dtypes
    """

    oge_dtypes = {
        "acid_gas_removal_efficiency": "float64",
        "annual_nox_emission_rate_lb_per_mmbtu": "float64",
        "ba_code": "string",
        "ba_code_physical": "string",
        "cems_profile": "float64",
        "ch4_mass_lb": "float64",
        "ch4_mass_lb_adjusted": "float64",
        "ch4_mass_lb_for_electricity": "float64",
        "ch4_mass_lb_for_electricity_adjusted": "float64",
        "co2_mass_lb": "float64",
        "co2_mass_lb_adjusted": "float64",
        "co2_mass_lb_for_electricity": "float64",
        "co2_mass_lb_for_electricity_adjusted": "float64",
        "co2_mass_measurement_code": "category",
        "co2e_mass_lb": "float64",
        "co2e_mass_lb_adjusted": "float64",
        "co2e_mass_lb_for_electricity": "float64",
        "co2e_mass_lb_for_electricity_adjusted": "float64",
        "data_availability": "category",
        "distribution_flag": "boolean",
        "eia930_profile": "float64",
        "equipment_tech_description": "string",
        "fgd_electricity_consumption_mwh": "float64",
        "fgd_sorbent_consumption_1000_tons": "float64",
        "flat_profile": "float32",
        "fuel_category": "string",
        "fuel_category_eia930": "string",
        "gross_generation_mwh": "float64",
        "gtn_method": "category",
        "hourly_data_source": "category",
        "hours_in_service": "float64",
        "imputed_profile": "float64",
        "mercury_control_id_eia": "string",
        "mercury_emission_rate_lb_per_trillion_btu": "float64",
        "mercury_removal_efficiency": "float64",
        "n2o_mass_lb": "float64",
        "n2o_mass_lb_adjusted": "float64",
        "n2o_mass_lb_for_electricity": "float64",
        "n2o_mass_lb_for_electricity_adjusted": "float64",
        "nox_mass_lb": "float64",
        "controlled_nox_mass_lb": "float64",
        "nox_mass_lb_adjusted": "float64",
        "nox_mass_lb_for_electricity": "float64",
        "nox_mass_lb_for_electricity_adjusted": "float64",
        "nox_mass_measurement_code": "category",
        "operating_time_hours": "float32",
        "ozone_season_nox_emission_rate_lb_per_mmbtu": "float64",
        "particulate_emission_rate_lb_per_mmbtu": "float64",
        "particulate_removal_efficiency_annual": "float64",
        "particulate_removal_efficiency_at_full_load": "float64",
        "plant_primary_fuel": "string",
        "plant_primary_fuel_from_capacity_mw": "string",
        "plant_primary_fuel_from_fuel_consumed_for_electricity_mmbtu": "string",
        "plant_primary_fuel_from_mode": "string",
        "plant_primary_fuel_from_net_generation_mwh": "string",
        "profile": "float64",
        "profile_method": "string",
        "residual_profile": "float64",
        "scaled_residual_profile": "float64",
        "shifted_residual_profile": "float64",
        "so2_mass_lb": "float64",
        "so2_mass_lb_adjusted": "float64",
        "so2_mass_lb_for_electricity": "float64",
        "so2_mass_lb_for_electricity_adjusted": "float64",
        "so2_mass_measurement_code": "category",
        "so2_removal_efficiency_annual": "float64",
        "so2_removal_efficiency_at_full_load": "float64",
        "steam_load_1000_lb": "float64",
        "subplant_id": "Int64",
        "subplant_primary_fuel": "string",
        "subplant_primary_fuel_from_capacity_mw": "string",
        "subplant_primary_fuel_from_fuel_consumed_for_electricity_mmbtu": "string",
        "subplant_primary_fuel_from_mode": "string",
        "subplant_primary_fuel_from_net_generation_mwh": "string",
    }

    pudl_dtypes = {
        "acid_gas_control": "boolean",
        "air_flow_100pct_load_cubic_feet_per_minute": "float64",
        "ash_content_pct": "float64",
        "ash_impoundment": "boolean",
        "ash_impoundment_lined": "boolean",
        "ash_impoundment_status": "string",
        "associated_combined_heat_power": "boolean",
        "balancing_authority_code_eia": "string",
        "balancing_authority_code_eia_consistent_rate": "float64",
        "balancing_authority_name_eia": "string",
        "bga_source": "string",
        "boiler_fuel_code_1": "string",
        "boiler_fuel_code_2": "string",
        "boiler_fuel_code_3": "string",
        "boiler_fuel_code_4": "string",
        "boiler_generator_assn_type_code": "string",
        "boiler_id": "string",
        "boiler_manufacturer": "string",
        "boiler_manufacturer_code": "string",
        "boiler_status": "string",
        "boiler_type": "string",
        "bypass_heat_recovery": "boolean",
        "capacity_factor": "float64",
        "capacity_mw": "float64",
        "carbon_capture": "boolean",
        "city": "string",
        "co2_mass_measurement_code": "string",
        "co2_mass_tons": "float64",
        "code": "string",
        "cofire_fuels": "boolean",
        "compliance_year_mercury": "Int64",
        "compliance_year_nox": "Int64",
        "compliance_year_particulate": "Int64",
        "compliance_year_so2": "Int64",
        "county": "string",
        "data_maturity": "string",
        "datum": "string",
        "deliver_power_transgrid": "boolean",
        "description": "string",
        "distributed_generation": "boolean",
        "duct_burners": "boolean",
        "efficiency_100pct_load": "float64",
        "efficiency_50pct_load": "float64",
        "emission_control_equipment_type_code": "string",
        "emission_control_id_eia": "string",
        "emission_control_id_pudl": "float64",
        "emission_control_id_type": "string",
        "emission_control_equipment_cost": "float64",
        "emissions_unit_id_epa": "string",
        "energy_source_1_transport_1": "string",
        "energy_source_1_transport_2": "string",
        "energy_source_1_transport_3": "string",
        "energy_source_2_transport_1": "string",
        "energy_source_2_transport_2": "string",
        "energy_source_2_transport_3": "string",
        "energy_source_code": "string",
        "energy_source_code_1": "string",
        "energy_source_code_2": "string",
        "energy_source_code_3": "string",
        "energy_source_code_4": "string",
        "energy_source_code_5": "string",
        "energy_source_code_6": "string",
        "energy_storage": "boolean",
        "energy_storage_capacity_mwh": "float64",
        "ferc_cogen_docket_no": "string",
        "ferc_cogen_status": "boolean",
        "ferc_exempt_wholesale_generator": "boolean",
        "ferc_exempt_wholesale_generator_docket_no": "string",
        "ferc_qualifying_facility": "boolean",
        "ferc_qualifying_facility_docket_no": "string",
        "ferc_small_power_producer": "boolean",
        "ferc_small_power_producer_docket_no": "string",
        "firing_rate_using_coal_tons_per_hour": "float64",
        "firing_rate_using_gas_mcf_per_hour": "float64",
        "firing_rate_using_oil_bbls_per_hour": "float64",
        "firing_rate_using_other_fuels": "float64",
        "firing_type_1": "string",
        "firing_type_2": "string",
        "firing_type_3": "string",
        "fluidized_bed_tech": "boolean",
        "fly_ash_reinjection": "boolean",
        "fuel_consumed_for_electricity_mmbtu": "float64",
        "fuel_consumed_for_electricity_units": "float64",
        "fuel_consumed_mmbtu": "float64",
        "fuel_consumed_units": "float64",
        "fuel_cost_from_eiaapi": "boolean",
        "fuel_cost_per_mmbtu": "float64",
        "fuel_cost_per_mwh": "float64",
        "fuel_mmbtu_per_unit": "float64",
        "fuel_type_code_aer": "string",
        "fuel_type_code_pudl": "string",
        "fuel_type_count": "Int64",
        "generator_id": "string",
        "generator_id_epa": "string",
        "grid_voltage_1_kv": "float64",
        "grid_voltage_2_kv": "float64",
        "grid_voltage_3_kv": "float64",
        "gross_load_mw": "float64",
        "heat_content_mmbtu": "float64",
        "hrsg": "boolean",
        "iso_rto_code": "string",
        "label": "string",
        "latitude": "float64",
        "liquefied_natural_gas_storage": "boolean",
        "longitude": "float64",
        "max_steam_flow_1000_lbs_per_hour": "float64",
        "mercury_control_existing_strategy_1": "string",
        "mercury_control_existing_strategy_2": "string",
        "mercury_control_existing_strategy_3": "string",
        "mercury_control_existing_strategy_4": "string",
        "mercury_control_existing_strategy_5": "string",
        "mercury_control_existing_strategy_6": "string",
        "mercury_control_id_eia": "string",
        "mercury_control_proposed_strategy_1": "string",
        "mercury_control_proposed_strategy_2": "string",
        "mercury_control_proposed_strategy_3": "string",
        "minimum_load_mw": "float64",
        "multiple_fuels": "boolean",
        "nameplate_power_factor": "float64",
        "natural_gas_local_distribution_company": "string",
        "natural_gas_pipeline_name_1": "string",
        "natural_gas_pipeline_name_2": "string",
        "natural_gas_pipeline_name_3": "string",
        "natural_gas_storage": "boolean",
        "nerc_region": "string",
        "net_capacity_mwdc": "float64",
        "net_generation_mwh": "float64",
        "net_metering": "boolean",
        "new_source_review": "boolean",
        "new_source_review_permit": "string",
        "nox_control_existing_caaa_compliance_strategy_1": "string",
        "nox_control_existing_caaa_compliance_strategy_2": "string",
        "nox_control_existing_caaa_compliance_strategy_3": "string",
        "nox_control_existing_strategy_1": "string",
        "nox_control_existing_strategy_2": "string",
        "nox_control_existing_strategy_3": "string",
        "nox_control_id_eia": "string",
        "nox_control_manufacturer": "string",
        "nox_control_manufacturer_code": "string",
        "nox_control_out_of_compliance_strategy_1": "string",
        "nox_control_out_of_compliance_strategy_2": "string",
        "nox_control_out_of_compliance_strategy_3": "string",
        "nox_control_planned_caaa_compliance_strategy_1": "string",
        "nox_control_planned_caaa_compliance_strategy_2": "string",
        "nox_control_planned_caaa_compliance_strategy_3": "string",
        "nox_control_proposed_strategy_1": "string",
        "nox_control_proposed_strategy_2": "string",
        "nox_control_proposed_strategy_3": "string",
        "nox_control_status_code": "string",
        "nox_mass_lbs": "float64",
        "nox_mass_measurement_code": "string",
        "operating_switch": "string",
        "operating_time_hours": "float64",
        "operational_status": "string",
        "operational_status_code": "string",
        "other_combustion_tech": "boolean",
        "other_planned_modifications": "boolean",
        "owned_by_non_utility": "boolean",
        "ownership_code": "string",
        "particulate_control_id_eia": "string",
        "particulate_control_out_of_compliance_strategy_1": "string",
        "particulate_control_out_of_compliance_strategy_2": "string",
        "particulate_control_out_of_compliance_strategy_3": "string",
        "pipeline_notes": "string",
        "planned_energy_source_code_1": "string",
        "planned_modifications": "boolean",
        "planned_net_summer_capacity_derate_mw": "float64",
        "planned_net_summer_capacity_uprate_mw": "float64",
        "planned_net_winter_capacity_derate_mw": "float64",
        "planned_net_winter_capacity_uprate_mw": "float64",
        "planned_new_capacity_mw": "float64",
        "planned_new_prime_mover_code": "string",
        "plant_id_eia": "Int64",
        "plant_id_epa": "Int64",
        "plant_id_pudl": "Int64",
        "plant_name_eia": "string",
        "previously_canceled": "boolean",
        "primary_purpose_id_naics": "Int64",
        "prime_mover_code": "string",
        "pulverized_coal_tech": "boolean",
        "reactive_power_output_mvar": "float64",
        "regulation_mercury": "string",
        "regulation_nox": "string",
        "regulation_particulate": "string",
        "regulation_so2": "string",
        "regulatory_status_code": "string",
        "report_year": "Int64",
        "reporting_frequency_code": "string",
        "rto_iso_lmp_node_id": "string",
        "rto_iso_location_wholesale_reporting_id": "string",
        "sector_id_eia": "Int64",
        "sector_name_eia": "string",
        "service_area": "string",
        "so2_control_existing_caaa_compliance_strategy_1": "string",
        "so2_control_existing_caaa_compliance_strategy_2": "string",
        "so2_control_existing_caaa_compliance_strategy_3": "string",
        "so2_control_existing_strategy_1": "string",
        "so2_control_existing_strategy_2": "string",
        "so2_control_existing_strategy_3": "string",
        "so2_control_id_eia": "string",
        "so2_control_out_of_compliance_strategy_1": "string",
        "so2_control_out_of_compliance_strategy_2": "string",
        "so2_control_out_of_compliance_strategy_3": "string",
        "so2_control_planned_caaa_compliance_strategy_1": "string",
        "so2_control_planned_caaa_compliance_strategy_2": "string",
        "so2_control_planned_caaa_compliance_strategy_3": "string",
        "so2_control_proposed_strategy_1": "string",
        "so2_control_proposed_strategy_2": "string",
        "so2_control_proposed_strategy_3": "string",
        "so2_mass_lbs": "float64",
        "so2_mass_measurement_code": "string",
        "solid_fuel_gasification": "boolean",
        "standard_nox_rate": "float64",
        "standard_particulate_rate": "float64",
        "standard_so2_percent_scrubbed": "float64",
        "standard_so2_rate": "float64",
        "startup_source_code_1": "string",
        "startup_source_code_2": "string",
        "startup_source_code_3": "string",
        "startup_source_code_4": "string",
        "state": "string",
        "steam_load_1000_lbs": "float64",
        "steam_plant_type_code": "Int64",
        "stoker_tech": "boolean",
        "street_address": "string",
        "subcritical_tech": "boolean",
        "sulfur_content_pct": "float64",
        "summer_capacity_estimate": "boolean",
        "summer_capacity_mw": "float64",
        "summer_estimated_capability_mw": "float64",
        "supercritical_tech": "boolean",
        "switch_oil_gas": "boolean",
        "syncronized_transmission_grid": "boolean",
        "technology_description": "string",
        "time_cold_shutdown_full_load_code": "string",
        "timezone": "string",
        "topping_bottoming_code": "string",
        "total_fuel_cost": "float64",
        "total_mmbtu": "float64",
        "transmission_distribution_owner_id": "Int64",
        "transmission_distribution_owner_name": "string",
        "transmission_distribution_owner_state": "string",
        "turbines_inverters_hydrokinetics": "Int64",
        "turbines_num": "Int64",
        "turndown_ratio": "float64",
        "ultrasupercritical_tech": "boolean",
        "unit_heat_rate_mmbtu_per_mwh": "float64",
        "unit_id_eia": "string",
        "unit_id_pudl": "Int64",
        "unit_nox": "string",
        "unit_particulate": "string",
        "unit_so2": "string",
        "uprate_derate_during_year": "boolean",
        "utility_id_eia": "Int64",
        "utility_id_pudl": "Int64",
        "utility_name_eia": "string",
        "waste_heat_input_mmbtu_per_hour": "float64",
        "water_source": "string",
        "wet_dry_bottom": "string",
        "winter_capacity_estimate": "boolean",
        "winter_capacity_mw": "float64",
        "winter_estimated_capability_mw": "float64",
        "year": "Int64",
        "zip_code": "string",
    }

    pudl_dtypes.update(oge_dtypes)

    return pudl_dtypes


def apply_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Applies specified types to a data frame and identifies if a type is not
    specified for a column.

    Args:
        df (pd.DataFrame): table whose columns will be converted.

    Returns:
        pd.DataFrame: original data frame with type converted columns.
    """
    dtypes = get_dtypes()
    datetime_cols = [
        "boiler_operating_date",
        "boiler_retirement_date",
        "current_planned_generator_operating_date",
        "datetime_utc",
        "datetime_local",
        "emission_control_operating_date",
        "emission_control_retirement_date",
        "generator_operating_date",
        "generator_retirement_date",
        "new_source_review_date",
        "operating_datetime_utc",
        "original_planned_generator_operating_date",
        "other_modifications_date",
        "planned_derate_date",
        "planned_generator_retirement_date",
        "planned_repower_date",
        "planned_uprate_date",
        "plant_operating_date",
        "plant_retirement_date",
        "report_date",
        "uprate_derate_completed_date",
    ]

    cols_missing_dtypes = [
        col
        for col in df.columns
        if ((col not in dtypes) and (col not in datetime_cols))
    ]
    if len(cols_missing_dtypes) > 0:
        logger.warning(
            "The following columns do not have dtypes assigned in "
            "`column_checks.get_dtypes()`"
        )
        logger.warning(cols_missing_dtypes)

    # apply the dtype if the column is in the dtype dict
    df = df.astype({col: dtypes[col] for col in df.columns if (col in dtypes)})
    # apply datetime columns
    for col in datetime_cols:
        if col in df.columns:
            try:
                df[col] = df[col].astype("datetime64[s]")
            except TypeError:
                df[col] = df[col].dt.tz_localize(None).astype("datetime64[s]")
            if "_utc" in col:
                df[col] = df[col].dt.tz_localize("UTC")

    return df
