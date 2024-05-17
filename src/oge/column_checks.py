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
        "annual_subplant_shift_mw",
        "annual_subplant_ratio",
        "annual_plant_shift_mw",
        "annual_plant_ratio",
        "plant_primary_fuel",
        "annual_fuel_ratio",
        "default_gtn_ratio",
        "monthly_subplant_shift_mw",
        "monthly_subplant_ratio",
        "monthly_plant_ratio",
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


def get_dtypes():
    """Returns a dictionary of dtypes that should be used for each column name."""
    dtypes_to_use = {
        "acid_gas_removal_efficiency": "float64",
        "annual_nox_emission_rate_lb_per_mmbtu": "float64",
        "ba_code": "str",
        "ba_code_physical": "str",
        "boiler_id": "str",
        "capacity_mw": "float64",
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
        "distribution_flag": "bool",
        "eia930_profile": "float64",
        "emissions_unit_id_epa": "str",
        "energy_source_code": "str",
        "energy_source_code_1": "str",
        "equipment_tech_description": "str",
        "fgd_electricity_consumption_mwh": "float64",
        "fgd_sorbent_consumption_1000_tons": "float64",
        "firing_type_1": "str",
        "firing_type_2": "str",
        "firing_type_3": "str",
        "flat_profile": "float32",
        "fuel_category": "str",
        "fuel_category_eia930": "str",
        "fuel_consumed_for_electricity_mmbtu": "float64",
        "fuel_consumed_mmbtu": "float64",
        "fuel_mmbtu_per_unit": "float64",
        "generator_id": "str",
        "gross_generation_mwh": "float64",
        "gtn_method": "category",
        "hourly_data_source": "category",
        "hours_in_service": "float64",
        "imputed_profile": "float64",
        "mercury_control_id_eia": "str",
        "mercury_emission_rate_lb_per_trillion_btu": "float64",
        "mercury_removal_efficiency": "float64",
        "n2o_mass_lb": "float64",
        "n2o_mass_lb_adjusted": "float64",
        "n2o_mass_lb_for_electricity": "float64",
        "n2o_mass_lb_for_electricity_adjusted": "float64",
        "net_generation_mwh": "float64",
        "nox_control_id_eia": "str",
        "nox_mass_lb": "float64",
        "nox_mass_lb_adjusted": "float64",
        "nox_mass_lb_for_electricity": "float64",
        "nox_mass_lb_for_electricity_adjusted": "float64",
        "nox_mass_measurement_code": "category",
        "operating_time_hours": "float16",
        "operational_status": "str",
        "ozone_season_nox_emission_rate_lb_per_mmbtu": "float64",
        "particulate_control_id_eia": "str",
        "particulate_emission_rate_lb_per_mmbtu": "float64",
        "particulate_removal_efficiency_annual": "float64",
        "particulate_removal_efficiency_at_full_load": "float64",
        "plant_id_eia": "Int32",
        "plant_id_epa": "Int32",
        "plant_primary_fuel": "str",
        "plant_primary_fuel_from_capacity_mw": "str",
        "plant_primary_fuel_from_fuel_consumed_for_electricity_mmbtu": "str",
        "plant_primary_fuel_from_mode": "str",
        "plant_primary_fuel_from_net_generation_mwh": "str",
        "prime_mover_code": "str",
        "profile": "float64",
        "profile_method": "str",
        "residual_profile": "float64",
        "scaled_residual_profile": "float64",
        "shifted_residual_profile": "float64",
        "so2_control_id_eia": "str",
        "so2_mass_lb": "float64",
        "so2_mass_lb_adjusted": "float64",
        "so2_mass_lb_for_electricity": "float64",
        "so2_mass_lb_for_electricity_adjusted": "float64",
        "so2_mass_measurement_code": "category",
        "so2_removal_efficiency_annual": "float64",
        "so2_removal_efficiency_at_full_load": "float64",
        "state": "str",
        "steam_load_1000_lb": "float64",
        "subplant_id": "Int16",
        "subplant_primary_fuel": "str",
        "subplant_primary_fuel_from_capacity_mw": "str",
        "subplant_primary_fuel_from_fuel_consumed_for_electricity_mmbtu": "str",
        "subplant_primary_fuel_from_mode": "str",
        "subplant_primary_fuel_from_net_generation_mwh": "str",
        "timezone": "str",
        "wet_dry_bottom": "str",
        "latitude": "float64",
        "longitude": "float64",
        "county": "str",
        "city": "str",
        "plant_name_eia": "str",
    }

    return dtypes_to_use


def apply_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Applies specified types to a data frame and identifies if a type is not
    specified for a column.

    Args:
        df (pd.DataFrame): table whose columns will be converted.

    Returns:
        pd.DataFrame: original data frame with type converted columns.
    """
    dtypes = get_dtypes()
    datetime_columns = ["datetime_utc", "datetime_local", "report_date"]
    cols_missing_dtypes = [
        col
        for col in df.columns
        if (col not in dtypes) and (col not in datetime_columns)
    ]
    if len(cols_missing_dtypes) > 0:
        logger.warning(
            "The following columns do not have dtypes assigned in "
            "`column_checks.get_dtypes()`"
        )
        logger.warning(cols_missing_dtypes)
    return df.astype({col: dtypes[col] for col in df.columns if col in dtypes})
