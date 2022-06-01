"""
Entry point for creating final dataset and intermediate cleaned data products.

Run from `src` as `python data_pipeline.py` after installing conda environment

Optional arguments are --year (default 2020), --gtn_years (default 5)

# Overview of Data Pipeline 

## Data Used
EPA Continuous Emissions Monitoring System (CEMS) data
 - What is it: Measured hourly gross generation, fuel consumption, and emissions data for emitting power generation units > 25MW
 - How we use it: Primary source for hourly emissions and generation data

EIA Form 923
 - What is it: Reported monthly net generation and fuel consumption data for power generators > 1 MW
 - How we use it: To convert gross generation data from CEMS to net generation, and to calculate emissions that are not reported to CEMS

EIA Form 860
 - What is it: Inventory of all generators and plants and their static characteristics
 - How we use it: to transform and aggregate the data reported in CEMS and EIA-923 based on plant and generator characteristics

EPA-EIA Power Sector Data Crosswalk
 - What is it: Maps EPA plant IDs and unit IDs to EIA plant IDs and generator IDs
 - How we use it: To match data between CEMS and EIA-923

EIA Form 930 / Hourly Electric Grid Monitor
 - What is it: Reported hourly net generation by fuel category, demand, and interchange for each Balancing Area in the U.S. 
 - How we use it: To assign an hourly profile to the monthly generation and fuel data reported in EIA-923

EPA eGRID database
 - What is it: Reports annual-level generation and emissions statistics at the plant and BA level 
 - How we use it: to validate our outputs

## Process
1. Download data, including CEMS (via PUDL), EIA Forms 860 and 923 (via PUDL), EPA-EIA Power Sector Data Crosswalk, EIA-930 data
    - Downloads are cached on first run so do not need to be redownloaded  
2. Identify subplants and gross-to-net generation factors using multiple years of historical data.
    - Using Power Sector Data Crosswalk, identify distinct subplant clusters of EPA units and EIA generators in each plant
    - Using multiple years of generation data from CEMS and EIA-923, run linear regressions of net generation on gross generation at teh subplant and plant level
    - Calculate simple monthly ratios between gross and net generation at teh subplant and plant level.
3. Clean monthly generation and fuel data from EIA-923
    - allocate monthly net generation and fuel consumption data reported for each plant prime mover to each plant generator
    - Calculate monthly emissions for each generator based on its fuel consumption and fuel source
    - Remove data for non grid-connected plants and plants in Puerto Rico
    - Assign a primary fuel type and balancing authority location to each generator
4. Clean hourly generation, fuel, and emissions data from CEMS
    - Remove data for non grid-connected plants, plants in Puerto Rico, and certain steam-only units
    - Assign a monthly "report_date" to each hourly observation based on the date of the local timestamp (this allows us to match the data to EIA-923 report dates)
    - Assign a fuel type to each unit
    - Fill in missing hourly emissions data using the assigned fuel type and reported hourly fuel consumption data
    - Remove all observations for each unit-month when no operation is reported for that unit in that month
    - Allocate hourly data for combined heat and power plants between electricity generation and steam production
    - Remove data for units for which we are unable to fill missing emissions data
5. Convert hourly gross generation in CEMS to hourly net generation
    - aggregate CEMS gross generation to monthly level to match with monthly-reported net generation
    - Apply several methodologies to calculate gross-to-net generation conversion factors
    - apply GTN factors to convert hourly gross generation to hourly net generation using the following hierarchy:
        - Use regression value if regression has good r2
        - If there is not a good regrssion, use monthly ratio unless ratio is outside of normal bounds (negative, >>1, missing)
        - Where there are outliers (eg gross generation is very different from net generation):
            - if EIA reported monthly and not distributed, maybe trust EIA (monthly ratio).
            - Otherwise, trust general regression
6. Crosswalk the CEMS data to the EIA-923 data to identify for which generator-months there is no hourly data reported in CEMS
    - Use the EPA-EIA Power Sector Data Crosswalk
    - Assign subplant groupings to data
7. Assign static plant characteristics to CEMS and EIA data to allow for data aggregation and matching with EIA-930
    - assign generator and plant-level primary fuel
    - assign Balancing Authority and State to each plant
    - assign fuel categories to each plant that match EIA-930 categories
8. Clean and reconcile EIA-930 data
    - Fix timezone/timestamp issues with raw 930 data
    - Perform physics-based reconciliation so that data satisfies conservation of energy equations
9. Calculate residual net generation profiles for each BA-fuel category by comparing EIA-930 and CEMS hourly net generation data
10. Assign monthly EIA-923 data an hourly profile based on the residual net generation profile
11. Concatenate the shaped hourly EIA-923 data to the hourly CEMS data
12. Run validation checks on processed data
13. Aggregate the hourly data to the BA level and output


## Outputs
 - Processed hourly subplant-level data
 - Aggregated hourly data for each BA (total emissions, total generation, generated carbon intensity)

## Output Validation Checks
 - Aggregate data to annual level and compare with published eGRID results
 - Check that aggregated heat rates and emissions rates by fuel type are within reasonable ranges for each BA
 - Plant-level checks for anomolous data

"""


# import packages
import numpy as np
import pandas as pd
import argparse
import os

# # Tell python where to look for modules.
import sys

sys.path.append("../../hourly-egrid/")

# import local modules
import src.data_cleaning as data_cleaning
import src.gross_to_net_generation as gross_to_net_generation
import src.load_data as load_data
import src.residual as residual
import src.column_checks as column_checks


def get_args():
    """
    Specify arguments here. 
    Returns dictionary of {arg_name: arg_value}
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", help="Year for analysis", default=2020, type=int)
    parser.add_argument(
        "--gtn_years",
        help="Number of years to use to calculate GTN ratio regressions, ending at `year`",
        default=5,
        type=int,
    )
    parser.add_argument(
        "--small",
        help="Run on subset of data for quicker testing, outputs to outputs/small and results to results/small.",
        type=bool,
        default=False,
    )

    args = parser.parse_args()
    return args


def write_final_result(combined_data, path_prefix):
    """
    Helper function to write combined data by BA
    """

    # only keep relevant columns
    combined_data = combined_data[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
            "co2_mass_lb",
            "co2_mass_lb_adjusted",
            "data_source",
        ]
    ]

    for ba in list(combined_data.ba_code.unique()):

        # filter the data for a single BA
        ba_table = combined_data[combined_data["ba_code"] == ba].drop(columns="ba_code")

        # convert the datetime_utc column back to a datetime
        ba_table["datetime_utc"] = pd.to_datetime(ba_table["datetime_utc"], utc=True)

        # combine the data from CEMS and EIA for each fuel-hour
        ba_table = (
            ba_table.groupby(["fuel_category", "datetime_utc"]).sum().reset_index()
        )

        # calculate a total for the BA
        ba_total = (
            ba_table.groupby(["datetime_utc"])
            .sum()[
                [
                    "net_generation_mwh",
                    "fuel_consumed_mmbtu",
                    "fuel_consumed_for_electricity_mmbtu",
                    "co2_mass_lb",
                    "co2_mass_lb_adjusted",
                ]
            ]
            .reset_index()
        )
        ba_total["fuel_category"] = "total"

        # concat the totals to the fuel-specific totals
        ba_table = pd.concat([ba_table, ba_total], axis=0, ignore_index=True)

        # calculate a generated emission rate
        ba_table["generated_co2_rate_lb_per_mwh"] = (
            (ba_table["co2_mass_lb"] * 2000 / ba_table["net_generation_mwh"])
            .fillna(0)
            .replace(np.inf, np.NaN)
        )
        ba_table["adjusted_generated_co2_rate_lb_per_mwh"] = (
            (ba_table["co2_mass_lb_adjusted"] * 2000 / ba_table["net_generation_mwh"])
            .fillna(0)
            .replace(np.inf, np.NaN)
        )

        ba_table = ba_table.pivot(index="datetime_utc", columns="fuel_category")

        # round all values to one decimal place
        ba_table = ba_table.round(1)

        # flatten the multilevel column into a single column name like data_fuelname
        ba_table.columns = ["_".join(col) for col in ba_table.columns.values]

        # export to a csv
        ba_table.to_csv(f"../data/results/{path_prefix}carbon_accounting/{ba}.csv")


def main():
    args = get_args()
    year = args.year

    # 0. Set up directory structure
    path_prefix = "" if not args.small else "small/"
    os.makedirs("../data/downloads", exist_ok=True)
    os.makedirs(f"../data/outputs/{path_prefix}", exist_ok=True)
    os.makedirs(f"../data/results/{path_prefix}", exist_ok=True)
    os.makedirs(f"../data/results/{path_prefix}carbon_accounting", exist_ok=True)
    os.makedirs(f"../data/results/{path_prefix}plant_data", exist_ok=True)

    # 1. Download data
    # PUDL
    load_data.download_pudl_data(
        zenodo_url="https://zenodo.org/record/6349861/files/pudl-v0.6.0-2022-03-12.tgz"
    )
    load_data.download_updated_pudl_database(download=False)
    # eGRID
    # the 2019 and 2020 data appear to be hosted on different urls
    egrid_files_to_download = [
        "https://www.epa.gov/sites/default/files/2021-02/egrid2019_data.xlsx",
        "https://www.epa.gov/system/files/documents/2022-01/egrid2020_data.xlsx",
    ]
    load_data.download_egrid_files(egrid_files_to_download)
    # EIA-930
    load_data.download_eia930_data(years_to_download=[year])
    load_data.download_chalendar_files()
    # Power Sector Data Crosswalk
    # NOTE: Check for new releases at https://github.com/USEPA/camd-eia-crosswalk
    load_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.2.1/epa_eia_crosswalk.csv"
    )

    # 2. Identify subplants and gross-to net ratios
    # GTN ratios are saved for reloading, as this is computationally intensive
    if not os.path.isdir("../data/outputs/gross_to_net/"):
        print('Generating subplant IDs and gross to net calcuations')
        number_of_years = args.gtn_years
        gross_to_net_generation.identify_subplants_and_gtn_conversions(
            year, number_of_years
        )

    # 3. Clean EIA-923 Generation and Fuel Data at the Monthly Level
    print('Cleaning EIA-923 data')
    eia923_allocated, primary_fuel_table = data_cleaning.clean_eia923(year, args.small)

    # Add primary fuel data to each generator
    eia923_allocated = eia923_allocated.merge(
        primary_fuel_table,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )

    # 4. Clean Hourly Data from CEMS
    print('Cleaning CEMS data')
    cems = data_cleaning.clean_cems(year, args.small)

    # 5. Convert CEMS Hourly Gross Generation to Hourly Net Generation
    print('Converting CEMS gross generation to net generation')
    cems = data_cleaning.convert_gross_to_net_generation(cems)

    # 6. Crosswalk CEMS and EIA data
    print('Identifying source for hourly data')
    eia923_allocated = data_cleaning.identify_hourly_data_source(
        eia923_allocated, cems, year
    )

    # 7. Calculate hourly data for partial_cems plants
    print('Scaling partial CEMS data')
    partial_cems_scaled, eia923_allocated = data_cleaning.scale_partial_cems_data(
        cems, eia923_allocated
    )

    # Export data cleaned by above for later validation, visualization, analysis
    print('Exporting intermediate output files')
    cems.to_csv(f"../data/outputs/{path_prefix}cems_{year}.csv", index=False)
    column_checks.check_columns(f"../data/outputs/{path_prefix}cems_{year}.csv")
    eia923_allocated.to_csv(
        f"../data/outputs/{path_prefix}eia923_allocated_{year}.csv", index=False
    )
    column_checks.check_columns(
        f"../data/outputs/{path_prefix}eia923_allocated_{year}.csv"
    )
    partial_cems_scaled.to_csv(
        f"../data/outputs/{path_prefix}partial_cems_scaled_{year}.csv", index=False
    )
    column_checks.check_columns(
        f"../data/outputs/{path_prefix}partial_cems_scaled_{year}.csv"
    )

    # 8. Assign static characteristics to CEMS and EIA data to aid in aggregation
    # assign a BA code and state code to each plant
    eia923_allocated = data_cleaning.assign_ba_code_to_plant(eia923_allocated, year)
    # assign a fuel category to each plant based on what is most likely to match with the category used in EIA-930
    # TODO: Add two different fuel categories (one for 930, one that is more specific)
    eia923_allocated = data_cleaning.assign_fuel_category_to_ESC(
        df=eia923_allocated, esc_column="plant_primary_fuel",
    )
    # add a flag about whether the plant is distribution connected
    eia923_allocated = data_cleaning.identify_distribution_connected_plants(
        eia923_allocated, year, voltage_threshold_kv=60
    )
    # Repeat for CEMS
    cems = data_cleaning.assign_ba_code_to_plant(cems, year)
    cems = data_cleaning.identify_distribution_connected_plants(
        cems, year, voltage_threshold_kv=60
    )
    # add a plant primary fuel and a fuel category for eia930
    cems = cems.merge(
        primary_fuel_table.drop_duplicates(subset="plant_id_eia")[
            ["plant_id_eia", "plant_primary_fuel"]
        ],
        how="left",
        on="plant_id_eia",
    )
    cems = data_cleaning.assign_fuel_category_to_ESC(
        df=cems, esc_column="plant_primary_fuel"
    )

    partial_cems_scaled = data_cleaning.assign_ba_code_to_plant(
        partial_cems_scaled, year
    )
    # add a plant primary fuel and a fuel category for eia930
    partial_cems_scaled = partial_cems_scaled.merge(
        primary_fuel_table.drop_duplicates(subset="plant_id_eia")[
            ["plant_id_eia", "plant_primary_fuel"]
        ],
        how="left",
        on="plant_id_eia",
    )
    partial_cems_scaled = data_cleaning.assign_fuel_category_to_ESC(
        df=partial_cems_scaled, esc_column="plant_primary_fuel"
    )

    # export plant frame
    plant_static_columns = [
        "plant_id_eia",
        "plant_primary_fuel",
        "fuel_category",
        "fuel_category_eia930",
        "ba_code",
        "ba_code_physical",
        "state",
        "distribution_flag",
    ]
    plant_frame = eia923_allocated[plant_static_columns].drop_duplicates(
        subset="plant_id_eia"
    )
    plant_frame.to_csv(
        f"../data/outputs/{path_prefix}plant_static_attributes.csv", index=False
    )
    column_checks.check_columns(
        f"../data/outputs/{path_prefix}plant_static_attributes.csv"
    )
    plant_frame.to_csv(
        f"../data/results/{path_prefix}plant_data/plant_static_attributes.csv",
        index=False,
    )

    # 9. Clean and Reconcile EIA-930 data
    print('Cleaning EIA-930 data')
    # TODO
    # Load raw EIA-930 data, fix timestamp issues, perform physics-based reconciliation
    # Currently implemented in `notebooks/930_lag` and the `gridemissions` repository
    # Output: `data/outputs/EBA_adjusted_elec.csv`

    # 10. Calculate Residual Net Generation Profile
    print('Calculating residual net generation profiles from EIA-930')
    # TODO
    # Currently implemented in `notebooks/calculate_residual_net_generation`

    # 11. Assign hourly profile to monthly data
    print('Assigning hourly profile to monthly EIA-923 data')
    # create a separate dataframe containing only the generators for which we do not have CEMS data
    monthly_eia_data_to_distribute = eia923_allocated[
        (eia923_allocated["hourly_data_source"] == "eia")
        & ~(eia923_allocated["fuel_consumed_mmbtu"].isna())
    ]
    # load profile data and format for use in the pipeline
    # TODO: once this is in the pipeline (step 10), may not need to read file
    hourly_profiles = pd.read_csv(
        "../data/outputs/residual_profiles.csv", parse_dates=["report_date"]
    )
    hourly_profiles = residual.assign_flat_profiles(
        monthly_eia_data_to_distribute, hourly_profiles, year
    )
    hourly_eia_data = data_cleaning.distribute_monthly_eia_data_to_hourly(
        monthly_eia_data_to_distribute, hourly_profiles, "residual_scaled"
    )
    # Export data
    columns_for_output = [
        "ba_code",
        "fuel_category",
        "datetime_utc",
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "co2_mass_lb_adjusted",
    ]
    hourly_eia_data[columns_for_output].to_csv(
        f"../data/outputs/{path_prefix}hourly_data_distributed_from_eia_{year}.csv",
        index=False,
    )
    # column_checks.check_columns(
    #     f"../data/outputs/{path_prefix}hourly_data_distributed_from_eia_{year}.csv"
    # )

    # 12. Aggregate CEMS data to BA-fuel and combine with hourly shaped EIA data
    print('Outputting final results')
    cems_ba_fuel = (
        cems.groupby(["ba_code", "fuel_category_eia930", "operating_datetime_utc"])
        .sum()[
            [
                "gross_generation_mwh",
                "net_generation_mwh",
                "fuel_consumed_mmbtu",
                "fuel_consumed_for_electricity_mmbtu",
                "co2_mass_lb",
                "co2_mass_lb_adjusted",
            ]
        ]
        .reset_index()
    )
    cems_ba_fuel["data_source"] = "CEMS"
    # rename the datetime_utc column
    cems_ba_fuel = cems_ba_fuel.rename(
        columns={
            "operating_datetime_utc": "datetime_utc",
            "fuel_category_eia930": "fuel_category",
        }
    )
    combined_data = pd.concat(
        [cems_ba_fuel, hourly_eia_data.drop(columns=["datetime_local", "report_date"])],
        axis=0,
    )

    write_final_result(combined_data, path_prefix)


if __name__ == "__main__":
    main()
