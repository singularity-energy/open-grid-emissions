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
import argparse
import os

# import local modules
# # Tell python where to look for modules.
import sys

sys.path.append("../../hourly-egrid/")
import src.download_data as download_data
import src.data_cleaning as data_cleaning
import src.gross_to_net_generation as gross_to_net_generation
import src.impute_hourly_profiles as impute_hourly_profiles
import src.eia930 as eia930
import src.validation as validation
import src.output_data as output_data
import src.consumed as consumed


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
    parser.add_argument(
        "--skip_outputs",
        help="Skip outputting data to csv files for quicker testing.",
        type=bool,
        default=False,
    )

    args = parser.parse_args()
    return args


def main():
    args = get_args()
    year = args.year

    # 0. Set up directory structure
    path_prefix = "" if not args.small else "small/"
    path_prefix += f"{year}/"
    os.makedirs("../data/downloads", exist_ok=True)
    os.makedirs(f"../data/outputs/{path_prefix}", exist_ok=True)
    os.makedirs(f"../data/outputs/{path_prefix}/eia930", exist_ok=True)
    os.makedirs(f"../data/results/{path_prefix}", exist_ok=True)
    os.makedirs(f"../data/results/{path_prefix}data_quality_metrics", exist_ok=True)
    for unit in ["us_units", "metric_units"]:
        for time_resolution in output_data.TIME_RESOLUTIONS.keys():
            for subfolder in ["plant_data", "carbon_accounting", "power_sector_data"]:
                os.makedirs(
                    f"../data/results/{path_prefix}/{subfolder}/{time_resolution}/{unit}",
                    exist_ok=True,
                )

    # 1. Download data
    ####################################################################################
    print("1. Downloading data")
    # PUDL
    download_data.download_pudl_data(
        zenodo_url="https://zenodo.org/record/6349861/files/pudl-v0.6.0-2022-03-12.tgz"
    )
    # eGRID
    # the 2019 and 2020 data appear to be hosted on different urls
    egrid_files_to_download = [
        "https://www.epa.gov/sites/default/files/2020-03/egrid2018_data_v2.xlsx",
        "https://www.epa.gov/sites/default/files/2021-02/egrid2019_data.xlsx",
        "https://www.epa.gov/system/files/documents/2022-01/egrid2020_data.xlsx",
    ]
    download_data.download_egrid_files(egrid_files_to_download)
    # EIA-930
    # for `small` run, we'll only clean 1 week, so need chalander file for making profiles
    download_data.download_chalendar_files()
    # We use balance files for imputing missing hourly profiles. TODO use cleaned instead?
    # need last year for rolling data cleaning
    download_data.download_eia930_data(years_to_download=[year, year - 1])
    # Power Sector Data Crosswalk
    # NOTE: Check for new releases at https://github.com/USEPA/camd-eia-crosswalk
    download_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.2.1/epa_eia_crosswalk.csv"
    )
    # download the raw EIA-923 and EIA-860 files for use in NOx/SO2 calculations until integrated into pudl
    download_data.download_raw_eia860(year)
    download_data.download_raw_eia923(year)

    # 2. Identify subplants
    ####################################################################################
    print("2. Identifying subplant IDs")
    # GTN ratios are saved for reloading, as this is computationally intensive
    if not os.path.exists(f"../data/outputs/{year}/subplant_crosswalk.csv"):
        print("    Generating subplant IDs")
        number_of_years = args.gtn_years
        data_cleaning.identify_subplants(year, number_of_years)
    else:
        print("    Subplant IDs already created")

    # 3. Clean EIA-923 Generation and Fuel Data at the Monthly Level
    ####################################################################################
    print("3. Cleaning EIA-923 data")
    eia923_allocated, primary_fuel_table = data_cleaning.clean_eia923(year, args.small)
    # Add primary fuel data to each generator
    eia923_allocated = eia923_allocated.merge(
        primary_fuel_table,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )

    # 4. Clean Hourly Data from CEMS
    ####################################################################################
    print("4. Cleaning CEMS data")
    cems = data_cleaning.clean_cems(year, args.small)

    # calculate biomass-adjusted emissions while cems data is at the unit level
    cems = data_cleaning.adjust_emissions_for_biomass(cems)

    # 5. Assign static characteristics to CEMS and EIA data to aid in aggregation
    ####################################################################################
    print("5. Loading plant static attributes")
    plant_attributes = data_cleaning.create_plant_attributes_table(
        cems, eia923_allocated, year, primary_fuel_table
    )
    output_data.output_intermediate_data(
        plant_attributes,
        "plant_static_attributes",
        path_prefix,
        year,
        args.skip_outputs,
    )

    # 6. Crosswalk CEMS and EIA data
    ####################################################################################
    print("6. Identifying source for hourly data")
    eia923_allocated = data_cleaning.identify_hourly_data_source(
        eia923_allocated, cems, year
    )

    # 7. Aggregating CEMS data to subplant
    ####################################################################################
    print("7. Aggregating CEMS data from unit to subplant")
    # aggregate cems data to subplant level
    cems = data_cleaning.aggregate_cems_to_subplant(cems)

    # 8. Calculate hourly data for partial_cems plants
    ####################################################################################
    print("8. Shaping partial CEMS data")
    (
        cems,
        partial_cems,
    ) = impute_hourly_profiles.shape_partial_cems_data(cems, eia923_allocated)
    # Export data cleaned by above for later validation, visualization, analysis
    output_data.output_intermediate_data(
        eia923_allocated.drop(columns="plant_primary_fuel"),
        "eia923_allocated",
        path_prefix,
        year,
        args.skip_outputs,
    )
    output_data.output_intermediate_data(
        partial_cems, "partial_cems", path_prefix, year, args.skip_outputs
    )

    # 9. Convert CEMS Hourly Gross Generation to Hourly Net Generation
    ####################################################################################
    print("9. Converting CEMS gross generation to net generation")
    cems, gtn_conversions = gross_to_net_generation.convert_gross_to_net_generation(
        cems, eia923_allocated, plant_attributes, year
    )
    # calculate the percent of gross generation converted using each method
    output_data.output_data_quality_metrics(
        validation.identify_cems_gtn_method(cems),
        "cems_gross_to_net_methods",
        path_prefix,
        args.skip_outputs,
    )
    # export the gtn conversion data
    output_data.output_intermediate_data(
        gtn_conversions,
        "gross_to_net_conversions",
        path_prefix,
        year,
        args.skip_outputs,
    )

    # 10. Adjust CEMS emission data for CHP
    ####################################################################################
    print("10. Adjusting CEMS emissions for CHP")
    cems = data_cleaning.adjust_cems_for_chp(cems, eia923_allocated)
    cems = data_cleaning.calculate_co2e_mass(
        cems, year, gwp_horizon=100, ar5_climate_carbon_feedback=True
    )
    validation.test_emissions_adjustments(cems)
    output_data.output_intermediate_data(
        cems, "cems", path_prefix, year, args.skip_outputs
    )

    # 11. Export monthly and annual plant-level results
    ####################################################################################
    print("11. Exporting monthly and annual plant-level results")
    # create a separate dataframe containing only the EIA data that is missing from cems
    monthly_eia_data_to_shape = eia923_allocated[
        (eia923_allocated["hourly_data_source"] == "eia")
        & ~(eia923_allocated["fuel_consumed_mmbtu"].isna())
    ]
    del eia923_allocated
    output_data.output_data_quality_metrics(
        validation.identify_percent_of_data_by_input_source(
            cems, partial_cems, monthly_eia_data_to_shape, year
        ),
        "input_data_source",
        path_prefix,
        args.skip_outputs,
    )
    # combine and export plant data at monthly and annual level
    monthly_plant_data = data_cleaning.combine_plant_data(
        cems, partial_cems, monthly_eia_data_to_shape, "monthly", True
    )
    output_data.output_plant_data(
        monthly_plant_data, path_prefix, "monthly", args.skip_outputs
    )
    output_data.output_plant_data(
        monthly_plant_data, path_prefix, "annual", args.skip_outputs
    )
    del monthly_plant_data

    # 12. Clean and Reconcile EIA-930 data
    ####################################################################################
    print("12. Cleaning EIA-930 data")
    # Scrapes and cleans data in data/downloads, outputs cleaned file at EBA_elec.csv
    if args.small or not (
        os.path.exists(f"../data/outputs/{path_prefix}/eia930/eia930_elec.csv")
    ):
        eia930.clean_930(year, small=args.small, path_prefix=path_prefix)
    else:
        print(
            f"    Not re-running 930 cleaning. If you'd like to re-run, please delete data/outputs/{path_prefix}/eia930/"
        )
    # If running small, we didn't clean the whole year, so need to use the Chalender file to build residual profiles.
    clean_930_file = (
        "../data/downloads/eia930/chalendar/EBA_elec.csv"
        if args.small
        else f"../data/outputs/{path_prefix}/eia930/eia930_elec.csv"
    )
    eia930_data = eia930.load_chalendar_for_pipeline(clean_930_file, year=year)
    # until we can fix the physics reconciliation, we need to apply some post-processing steps
    eia930_data = eia930.remove_imputed_ones(eia930_data)
    eia930_data = eia930.remove_months_with_zero_data(eia930_data)

    # 13. Calculate hourly profiles for monthly EIA data
    ####################################################################################
    print("13. Estimating hourly profiles for EIA data")
    hourly_profiles = impute_hourly_profiles.calculate_hourly_profiles(
        cems,
        partial_cems,
        eia930_data,
        plant_attributes,
        monthly_eia_data_to_shape,
        year,
        transmission_only=False,
        ba_column_name="ba_code",
    )
    del eia930_data
    # validate how well the wind and solar imputation methods work
    output_data.output_data_quality_metrics(
        validation.validate_diba_imputation_method(hourly_profiles, year),
        "diba_imputation_performance",
        path_prefix,
        args.skip_outputs,
    )
    output_data.output_data_quality_metrics(
        validation.validate_national_imputation_method(hourly_profiles),
        "national_imputation_performance",
        path_prefix,
        args.skip_outputs,
    )
    output_data.output_intermediate_data(
        hourly_profiles, "hourly_profiles", path_prefix, year, args.skip_outputs
    )

    # 14. Assign hourly profile to monthly data
    ####################################################################################
    print("14. Assigning hourly profiles to monthly EIA-923 data")
    hourly_profiles = impute_hourly_profiles.convert_profile_to_percent(hourly_profiles)
    # Aggregate EIA data to BA/fuel/month, then assign hourly profile per BA/fuel
    (
        monthly_eia_data_to_shape,
        plant_attributes,
    ) = impute_hourly_profiles.aggregate_eia_data_to_ba_fuel(
        monthly_eia_data_to_shape, plant_attributes, path_prefix
    )
    shaped_eia_data = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
        monthly_eia_data_to_shape, hourly_profiles
    )
    output_data.output_data_quality_metrics(
        validation.hourly_profile_source_metric(cems, partial_cems, shaped_eia_data),
        "hourly_profile_method",
        path_prefix,
        args.skip_outputs,
    )
    # Export data
    output_data.output_intermediate_data(
        shaped_eia_data, "shaped_eia923_data", path_prefix, year, args.skip_outputs
    )
    if not args.skip_outputs:
        plant_attributes.to_csv(
            f"../data/results/{path_prefix}plant_data/plant_static_attributes.csv"
        )
    # validate that the shaping did not alter data at the monthly level
    validation.validate_shaped_totals(shaped_eia_data, monthly_eia_data_to_shape)

    # 15. Combine plant-level data from all sources
    ####################################################################################
    print("15. Combining and exporting plant-level hourly results")
    # write metadata and remove metadata columns
    cems, partial_cems, shaped_eia_data = output_data.write_plant_metadata(
        cems, partial_cems, shaped_eia_data, path_prefix, args.skip_outputs
    )
    combined_plant_data = data_cleaning.combine_plant_data(
        cems, partial_cems, shaped_eia_data, "hourly"
    )
    del shaped_eia_data, cems, partial_cems  # free memory back to python
    # export to a csv.
    output_data.output_plant_data(
        combined_plant_data, path_prefix, "hourly", args.skip_outputs
    )

    # 16. Aggregate CEMS data to BA-fuel and write power sector results
    ####################################################################################
    print("16. Creating and exporting BA-level power sector results")
    ba_fuel_data = data_cleaning.aggregate_plant_data_to_ba_fuel(
        combined_plant_data, plant_attributes
    )
    del combined_plant_data
    # Output intermediate data: produced per-fuel annual averages
    output_data.write_generated_averages(
        ba_fuel_data, year, path_prefix, args.skip_outputs
    )
    # Output final data: per-ba hourly generation and rate
    output_data.write_power_sector_results(ba_fuel_data, path_prefix, args.skip_outputs)

    # 17. Calculate consumption-based emissions and write carbon accounting results
    ####################################################################################
    print("17. Calculating and exporting consumption-based results")
    hourly_consumed_calc = consumed.HourlyBaDataEmissionsCalc(
        clean_930_file,
        year=year,
        small=args.small,
        path_prefix=path_prefix,
    )
    hourly_consumed_calc.process()
    hourly_consumed_calc.output_data(path_prefix=path_prefix)


if __name__ == "__main__":
    main()
