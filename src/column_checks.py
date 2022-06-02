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
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "subplant_id",
        "prime_mover_code",
        "energy_source_code",
        "plant_primary_fuel",
        "hourly_data_source",
    },
    "cems_": {
        "plant_id_eia",
        "unitid",
        "cems_id",
        "operating_datetime_utc",
        "operating_time_hours",
        "gross_generation_mwh",
        "steam_load_1000_lb",
        "fuel_consumed_mmbtu",
        "co2_mass_lb",
        "co2_mass_measurement_code",
        "nox_mass_lb",
        "nox_mass_measurement_code",
        "so2_mass_lb",
        "so2_mass_measurement_code",
        "plant_id_epa",
        "unit_id_epa",
        "report_date",
        "energy_source_code",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "fuel_consumed_for_electricity_mmbtu",
        "electric_allocation_factor",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
        "subplant_id",
        "gtn_method",
        "net_generation_mwh",
    },
    "partial_cems_scaled_": {
        "report_date",
        "plant_id_eia",
        "subplant_id",
        "operating_datetime_utc",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "net_generation_mwh",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
    },
    "plant_static_attributes": {
        "plant_id_eia",
        "plant_primary_fuel",
        "fuel_category",
        "fuel_category_eia930",
        "ba_code",
        "ba_code_physical",
        "state",
        "distribution_flag",
    },
    "residual_profiles_": {
        "ba_code",
        "fuel_category",
        "datetime_utc",
        "datetime_local",
        "report_date",
        "residual_scaled",
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
