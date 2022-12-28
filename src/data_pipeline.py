"""
Entry point for creating final dataset and intermediate cleaned data products.

Run from `src` as `python data_pipeline.py` after installing conda environment

Optional arguments are --year (default 2021), --shape_individual_plants (default True)
Optional arguments for development are --small, --flat, and --skip_outputs
"""


# import packages
import argparse
import os
import shutil

# import local modules
# import local modules
# # # Tell python where to look for modules.
import download_data
import data_cleaning
import emissions
import gross_to_net_generation
import impute_hourly_profiles
import eia930
import validation
import output_data
import consumed
from filepaths import downloads_folder, outputs_folder, results_folder


def get_args():
    """
    Specify arguments here.
    Returns dictionary of {arg_name: arg_value}
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", help="Year for analysis", default=2021, type=int)
    parser.add_argument(
        "--shape_individual_plants",
        help="Assign an hourly profile to each individual plant with EIA-only data, instead of aggregating to the fleet level before shaping.",
        type=bool,
        default=True,
    )
    parser.add_argument(
        "--small",
        help="Run on subset of data for quicker testing, outputs to outputs/small and results to results/small.",
        type=bool,
        default=False,
    )
    parser.add_argument(
        "--flat",
        help="Use flat hourly profiles?",
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

    validation.validate_year(year)

    # 0. Set up directory structure
    path_prefix = "" if not args.small else "small/"
    path_prefix += "flat/" if args.flat else ""
    path_prefix += f"{year}/"
    os.makedirs(downloads_folder(), exist_ok=True)
    os.makedirs(outputs_folder(f"{path_prefix}"), exist_ok=True)
    os.makedirs(outputs_folder(f"{path_prefix}/eia930"), exist_ok=True)
    if not args.skip_outputs:
        # If we are outputing, wipe results dir so we can be confident there are no old result files (eg because of a file name change)
        if os.path.exists(results_folder(f"{path_prefix}")):
            shutil.rmtree(results_folder(f"{path_prefix}"))
        os.makedirs(results_folder(f"{path_prefix}"), exist_ok=False)
    else:  # still make sure results dir exists, but exist is ok and we won't be writing to it
        os.makedirs(results_folder(f"{path_prefix}"), exist_ok=True)
    os.makedirs(
        results_folder(f"{path_prefix}data_quality_metrics"),
        exist_ok=True,
    )
    # Make results subfolders
    for unit in ["us_units", "metric_units"]:
        for time_resolution in output_data.TIME_RESOLUTIONS.keys():
            for subfolder in ["plant_data", "carbon_accounting", "power_sector_data"]:
                os.makedirs(
                    results_folder(
                        f"{path_prefix}/{subfolder}/{time_resolution}/{unit}"
                    ),
                    exist_ok=True,
                )

    # 1. Download data
    ####################################################################################
    print("1. Downloading data")
    # PUDL
    download_data.download_pudl_data(
        zenodo_url="https://zenodo.org/record/7472137/files/pudl-v2022.11.30.tgz"
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
    if args.small or args.flat:
        download_data.download_chalendar_files()
    # We use balance files for imputing missing hourly profiles. TODO use cleaned instead?
    # need last year for rolling data cleaning
    download_data.download_eia930_data(years_to_download=[year, year - 1])
    # Power Sector Data Crosswalk
    # NOTE: Check for new releases at https://github.com/USEPA/camd-eia-crosswalk
    download_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.3/epa_eia_crosswalk.csv"
    )
    # download the raw EIA-923 and EIA-860 files for use in NOx/SO2 calculations until integrated into pudl
    download_data.download_raw_eia860(year)
    download_data.download_raw_eia923(year)

    # 2. Identify subplants
    ####################################################################################
    print("2. Identifying subplant IDs")
    data_cleaning.identify_subplants(year)

    # 3. Clean EIA-923 Generation and Fuel Data at the Monthly Level
    ####################################################################################
    print("3. Cleaning EIA-923 data")
    (
        eia923_allocated,
        primary_fuel_table,
        subplant_emission_factors,
    ) = data_cleaning.clean_eia923(year, args.small)
    # Add primary fuel data to each generator
    eia923_allocated = eia923_allocated.merge(
        primary_fuel_table,
        how="left",
        on=["plant_id_eia", "subplant_id", "generator_id"],
        validate="m:1",
    )

    # 4. Clean Hourly Data from CEMS
    ####################################################################################
    print("4. Cleaning CEMS data")
    cems = data_cleaning.clean_cems(
        year, args.small, primary_fuel_table, subplant_emission_factors
    )
    # output data quality metrics about measured vs imputed CEMS data
    output_data.output_data_quality_metrics(
        validation.summarize_cems_measurement_quality(cems),
        "cems_pollutant_measurement_quality",
        path_prefix,
        args.skip_outputs,
    )

    # output cleaned cems data
    output_data.output_intermediate_data(
        cems,
        "cems_cleaned",
        path_prefix,
        year,
        args.skip_outputs,
    )

    # calculate biomass-adjusted emissions while cems data is at the unit level
    cems = emissions.adjust_emissions_for_biomass(cems)

    # 5. Assign static characteristics to CEMS and EIA data to aid in aggregation
    ####################################################################################
    print("5. Loading plant static attributes")
    plant_attributes = data_cleaning.create_plant_attributes_table(
        cems, eia923_allocated, year, primary_fuel_table
    )

    # 6. Crosswalk CEMS and EIA data
    ####################################################################################
    print("6. Identifying source for hourly data")
    eia923_allocated = data_cleaning.identify_hourly_data_source(
        eia923_allocated, cems, year
    )
    # Export data cleaned by above for later validation, visualization, analysis
    output_data.output_intermediate_data(
        eia923_allocated.drop(columns=["plant_primary_fuel", "subplant_primary_fuel"]),
        "eia923_allocated",
        path_prefix,
        year,
        args.skip_outputs,
    )
    # output data quality metrics about annually-reported EIA-923 data
    output_data.output_data_quality_metrics(
        validation.summarize_annually_reported_eia_data(eia923_allocated, year),
        "annually_reported_eia_data",
        path_prefix,
        args.skip_outputs,
    )

    # 7. Aggregating CEMS data to subplant
    ####################################################################################
    print("7. Aggregating CEMS data from unit to subplant")
    # aggregate cems data to subplant level
    cems = data_cleaning.aggregate_cems_to_subplant(cems)

    # 8. Calculate hourly data for partial_cems plants
    ####################################################################################
    print("8. Shaping partial CEMS data")
    # shape partial CEMS plant data
    partial_cems_plant = impute_hourly_profiles.shape_partial_cems_plants(
        cems, eia923_allocated
    )
    validation.validate_unique_datetimes(
        df=partial_cems_plant,
        df_name="partial_cems_plant",
        keys=["plant_id_eia", "subplant_id"],
    )
    output_data.output_intermediate_data(
        partial_cems_plant,
        "partial_cems_plant",
        path_prefix,
        year,
        args.skip_outputs,
    )
    # shape partial CEMS subplant data
    (
        cems,
        partial_cems_subplant,
    ) = impute_hourly_profiles.shape_partial_cems_subplants(cems, eia923_allocated)

    validation.validate_unique_datetimes(
        df=partial_cems_subplant,
        df_name="partial_cems_subplant",
        keys=["plant_id_eia", "subplant_id"],
    )
    output_data.output_intermediate_data(
        partial_cems_subplant,
        "partial_cems_subplant",
        path_prefix,
        year,
        args.skip_outputs,
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
    cems = emissions.calculate_co2e_mass(
        cems, year, gwp_horizon=100, ar5_climate_carbon_feedback=True
    )
    validation.test_emissions_adjustments(cems)
    validation.validate_unique_datetimes(
        df=cems,
        df_name="cems_subplant",
        keys=["plant_id_eia", "subplant_id"],
    )
    output_data.output_intermediate_data(
        cems, "cems_subplant", path_prefix, year, args.skip_outputs
    )

    # 11. Export monthly and annual plant-level results
    ####################################################################################
    print("11. Exporting monthly and annual plant-level results")
    # create a separate dataframe containing only the EIA data that is missing from cems
    monthly_eia_data_to_shape = eia923_allocated[
        (eia923_allocated["hourly_data_source"] == "eia")
        & ~(eia923_allocated["fuel_consumed_mmbtu"].isna())
    ]
    output_data.output_data_quality_metrics(
        validation.identify_percent_of_data_by_input_source(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            monthly_eia_data_to_shape,
            year,
            plant_attributes,
        ),
        "input_data_source",
        path_prefix,
        args.skip_outputs,
    )
    # combine and export plant data at monthly and annual level
    monthly_plant_data = data_cleaning.combine_plant_data(
        cems,
        partial_cems_subplant,
        partial_cems_plant,
        monthly_eia_data_to_shape,
        "monthly",
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
    if args.flat:
        print("    Not running 930 cleaning because we'll be using a flat profile.")
    elif not (os.path.exists(outputs_folder(f"{path_prefix}/eia930/eia930_elec.csv"))):
        eia930.clean_930(year, small=args.small, path_prefix=path_prefix)
    else:
        print(
            f"    Not re-running 930 cleaning. If you'd like to re-run, please delete data/outputs/{path_prefix}/eia930/"
        )

    # If running small, we didn't clean the whole year, so need to use the Chalender file to build residual profiles.
    clean_930_file = (
        downloads_folder("eia930/chalendar/EBA_elec.csv")
        if (args.small or args.flat)
        else outputs_folder(f"{path_prefix}/eia930/eia930_elec.csv")
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
        partial_cems_subplant,
        partial_cems_plant,
        eia930_data,
        plant_attributes,
        monthly_eia_data_to_shape,
        year,
        transmission_only=False,
        ba_column_name="ba_code",
        use_flat=args.flat,
    )
    del eia930_data
    # validate how well the wind and solar imputation methods work
    output_data.output_data_quality_metrics(
        validation.validate_wind_solar_imputation(hourly_profiles, year),
        "wind_solar_profile_imputation_performance",
        path_prefix,
        args.skip_outputs,
    )
    output_data.output_intermediate_data(
        hourly_profiles, "hourly_profiles", path_prefix, year, args.skip_outputs
    )

    hourly_profiles = impute_hourly_profiles.convert_profile_to_percent(
        hourly_profiles,
        group_keys=["ba_code", "fuel_category", "profile_method"],
        columns_to_convert=["profile", "flat_profile"],
    )

    # 14. Export hourly plant-level data
    ####################################################################################
    print("14. Exporting Hourly Plant-level data for each BA")
    if args.shape_individual_plants and not args.small:
        impute_hourly_profiles.combine_and_export_hourly_plant_data(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            monthly_eia_data_to_shape,
            plant_attributes,
            hourly_profiles,
            path_prefix,
            args.skip_outputs,
            region_to_group="ba_code",
        )
    else:
        print(
            "    Not shaping and exporting individual plant data since `shape_individual_plants` is False."
        )
        print(
            "    Plants that only report to EIA will be aggregated to the fleet level before shaping."
        )

    # 15. Shape fleet-level data
    ####################################################################################
    print("15. Assigning hourly profiles to monthly EIA-923 data")
    hourly_profiles = impute_hourly_profiles.convert_profile_to_percent(
        hourly_profiles,
        group_keys=["ba_code", "fuel_category", "profile_method"],
        columns_to_convert=["profile", "flat_profile"],
    )
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
        validation.hourly_profile_source_metric(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            shaped_eia_data,
            plant_attributes,
        ),
        "hourly_profile_method",
        path_prefix,
        args.skip_outputs,
    )
    # Export data
    validation.validate_unique_datetimes(
        df=shaped_eia_data,
        df_name="shaped_eia_data",
        keys=["plant_id_eia"],
    )
    output_data.output_intermediate_data(
        shaped_eia_data, "shaped_eia923_data", path_prefix, year, args.skip_outputs
    )
    output_data.output_intermediate_data(
        plant_attributes,
        "plant_static_attributes",
        path_prefix,
        year,
        args.skip_outputs,
    )
    if not args.skip_outputs:
        plant_attributes.to_csv(
            results_folder(f"{path_prefix}plant_data/plant_static_attributes.csv"),
            index=False,
        )
    # validate that the shaping did not alter data at the monthly level
    validation.validate_shaped_totals(
        shaped_eia_data,
        monthly_eia_data_to_shape,
        group_keys=["ba_code", "fuel_category"],
    )

    # 16. Combine plant-level data from all sources
    ####################################################################################
    print("16. Combining plant-level hourly data")
    # write metadata outputs
    output_data.write_plant_metadata(
        plant_attributes,
        eia923_allocated,
        cems,
        partial_cems_subplant,
        partial_cems_plant,
        shaped_eia_data,
        path_prefix,
        args.skip_outputs,
    )
    # set validate parameter to False since validating non-overlapping data requires subplant-level data
    # since the shaped eia data is at the fleet level, this check will not work.
    # However, we already checked for non-overlapping data in step 11 when combining monthly data
    combined_plant_data = data_cleaning.combine_plant_data(
        cems,
        partial_cems_subplant,
        partial_cems_plant,
        shaped_eia_data,
        "hourly",
        False,
    )
    del (
        shaped_eia_data,
        cems,
        partial_cems_subplant,
        partial_cems_plant,
    )  # free memory back to python
    # export to a csv.
    validation.validate_unique_datetimes(
        df=combined_plant_data,
        df_name="combined_plant_data",
        keys=["plant_id_eia"],
    )
    if not args.shape_individual_plants:
        output_data.output_plant_data(
            combined_plant_data,
            path_prefix,
            "hourly",
            args.skip_outputs,
        )

    # 17. Aggregate CEMS data to BA-fuel and write power sector results
    ####################################################################################
    print("17. Creating and exporting BA-level power sector results")
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

    # 18. Calculate consumption-based emissions and write carbon accounting results
    ####################################################################################
    print("18. Calculating and exporting consumption-based results")
    hourly_consumed_calc = consumed.HourlyConsumed(
        clean_930_file,
        path_prefix,
        year,
        small=args.small,
        skip_outputs=args.skip_outputs,
    )
    hourly_consumed_calc.run()
    hourly_consumed_calc.output_results()


if __name__ == "__main__":
    main()
