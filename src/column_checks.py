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

COLUMNS = {
    "eia923_allocated_": {
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
    "cems_": {
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
        # columns dropped when aggregating to subplant
        # "plant_id_epa",
        # "unitid",
        # "energy_source_code",
        # "operating_time_hours",
        # "co2_mass_measurement_code",
        # "nox_mass_measurement_code",
        # "so2_mass_measurement_code",
    },
    "partial_cems_": {
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
    "plant_static_attributes_": {
        "plant_id_eia",
        "plant_primary_fuel",
        "fuel_category",
        "fuel_category_eia930",
        "ba_code",
        "ba_code_physical",
        "state",
        "distribution_flag",
        "timezone",
        "data_availability",
    },
    "plant_attributes_with_synthetic_": {
        "plant_id_eia",
        "plant_primary_fuel",
        "fuel_category",
        "fuel_category_eia930",
        "ba_code",
        "ba_code_physical",
        "state",
        "distribution_flag",
        "timezone",
        "data_availability",
    },
    "hourly_profiles_": {
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
    "shaped_eia923_data_": {
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
    "annual_generation_averages_by_fuel_": {
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
        "generated_co2_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_nox_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_so2_rate_lb_per_mwh_for_electricity_adjusted",
    },
    "gross_to_net_conversions_": {
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
        "monthly_subplant_shift_mw",
        "monthly_subplant_ratio",
        "monthly_plant_ratio",
        "subplant_regression_ratio",
        "subplant_regression_shift_mw",
        "plant_regression_ratio",
        "plant_regression_shift_mw",
    },
}


def check_columns(file_path):
    """
    Given a file name or path to file, check that its columns are as expected.
    """
    file = file_path.split("/")[-1]
    file = file.replace(".csv", "")

    # If file is appended by year, remove it because column names are standard across years
    maybe_year = file[-4:]
    if maybe_year.isnumeric():
        year = int(maybe_year)
        if not 2000 < year < 2050:
            print(f"Got unexpected year {maybe_year}")
        file = file.replace(maybe_year, "")

    # Get actual columns
    with open(file_path) as f:
        firstline = f.readline().rstrip()
    cols = set(firstline.split(","))

    # Get expected columns
    if file not in COLUMNS:
        raise ValueError(
            f"Could not find file prefix {file} from {file_path} in expected file names {COLUMNS.keys()}"
        )
    expected_cols = COLUMNS[file]

    # Check for extra columns. Warning not exception
    extras = cols - expected_cols
    if len(extras) > 0:
        print(
            f"Warning: columns {extras} in {file_path} are not guaranteed by column_checks.py"
        )

    # Raise exception for missing columns
    missing = expected_cols - cols
    if len(missing) > 0:
        raise ValueError(f"Columns {missing} missing from {file_path}")

    return


def get_dtypes():
    dtypes_to_use = {
        "plant_id_eia": "Int32",
        "plant_id_epa": "Int32",
        "subplant_id": "Int16",
        "generator_id": "str",
        "unitid": "str",
        "operating_time_hours": "float16",
        "gross_generation_mwh": "float64",
        "steam_load_1000_lb": "float64",
        "fuel_consumed_mmbtu": "float64",
        "co2_mass_lb": "float64",
        "co2_mass_measurement_code": "category",
        "nox_mass_lb": "float64",
        "nox_mass_measurement_code": "category",
        "so2_mass_lb": "float64",
        "so2_mass_measurement_code": "category",
        "energy_source_code": "str",
        "ch4_mass_lb": "float64",
        "n2o_mass_lb": "float64",
        "fuel_consumed_for_electricity_mmbtu": "float64",
        "co2_mass_lb_for_electricity": "float64",
        "ch4_mass_lb_for_electricity": "float64",
        "n2o_mass_lb_for_electricity": "float64",
        "nox_mass_lb_for_electricity": "float64",
        "so2_mass_lb_for_electricity": "float64",
        "co2_mass_lb_adjusted": "float64",
        "ch4_mass_lb_adjusted": "float64",
        "n2o_mass_lb_adjusted": "float64",
        "nox_mass_lb_adjusted": "float64",
        "so2_mass_lb_adjusted": "float64",
        "co2_mass_lb_for_electricity_adjusted": "float64",
        "ch4_mass_lb_for_electricity_adjusted": "float64",
        "n2o_mass_lb_for_electricity_adjusted": "float64",
        "nox_mass_lb_for_electricity_adjusted": "float64",
        "so2_mass_lb_for_electricity_adjusted": "float64",
        "co2e_mass_lb": "float64",
        "co2e_mass_lb_for_electricity": "float64",
        "co2e_mass_lb_adjusted": "float64",
        "co2e_mass_lb_for_electricity_adjusted": "float64",
        "gtn_method": "category",
        "net_generation_mwh": "float64",
        "prime_mover_code": "str",
        "hourly_data_source": "category",
        "fuel_category": "str",
        "ba_code": "str",
        "ba_code_physical": "str",
        "plant_primary_fuel": "str",
        "fuel_category_eia930": "str",
        "state": "str",
        "distribution_flag": "bool",
        "timezone": "str",
        "eia930_profile": "float64",
        "cems_profile": "float64",
        "residual_profile": "float64",
        "scaled_residual_profile": "float64",
        "shifted_residual_profile": "float64",
        "imputed_profile": "float64",
        "profile": "float64",
        "flat_profile": "float32",
        "profile_method": "str",
        "data_availability": "category",
    }

    return dtypes_to_use


def apply_dtypes(df):
    dtypes = get_dtypes()
    datetime_columns = ["datetime_utc", "datetime_local", "report_date"]
    cols_missing_dtypes = [
        col
        for col in df.columns
        if (col not in dtypes) and (col not in datetime_columns)
    ]
    if len(cols_missing_dtypes) > 0:
        print(
            "Warning: The following columns do not have dtypes assigned in `column_checks.get_dtypes()`"
        )
        print(cols_missing_dtypes)
    return df.astype({col: dtypes[col] for col in df.columns if col in dtypes})
