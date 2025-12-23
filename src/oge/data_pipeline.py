"""
Entry point for creating final dataset and intermediate cleaned data products.

Run from `src` as `python data_pipeline.py` after installing conda environment

Optional arguments are --year (default 2024)
Optional arguments for development are --flat and --skip_outputs
"""

import argparse
import gc
import os
import pandas as pd
import shutil

# import local modules
import oge.data_cleaning as data_cleaning
import oge.consumed as consumed
import oge.download_data as download_data
import oge.eia930 as eia930
import oge.emissions as emissions
import oge.gross_to_net_generation as gross_to_net_generation
import oge.helpers as helpers
import oge.impute_hourly_profiles as impute_hourly_profiles
import oge.output_data as output_data
import oge.subplant_identification as subplant_identification
import oge.validation as validation

from oge import PUDL_ENGINE
from oge.constants import (
    TIME_RESOLUTIONS,
    latest_validated_year,
    current_early_release_year,
    earliest_hourly_data_year,
)
from oge.column_checks import DATA_COLUMNS
from oge.filepaths import downloads_folder, outputs_folder, results_folder
from oge.logging_util import get_logger, configure_root_logger


def get_args() -> argparse.Namespace:
    """Specify arguments here.

    Returns dictionary of {arg_name: arg_value}
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", help="Year for analysis", default=2024, type=int)
    parser.add_argument(
        "--flat",
        help="Use flat hourly profiles?",
        default=False,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--skip_outputs",
        help="Skip outputting data to csv files for quicker testing.",
        default=False,
        action=argparse.BooleanOptionalAction,
    )

    args = parser.parse_args()

    return args


def print_args(args: argparse.Namespace, logger):
    """Print out the command line arguments."""
    argstring = "\n".join([f"  * {k} = {v}" for k, v in vars(args).items()])
    logger.info(f"\n\nRunning with the following options:\n{argstring}\n")


def main(args):
    """Runs the OGE data pipeline."""
    args = get_args()
    year = args.year

    validation.validate_year(year)

    if os.getenv("OGE_DATA_STORE") in ["s3", "2"]:
        raise OSError(
            "Invalid OGE_DATA_STORE environment variable. Should be 'local' or '1'"
        )
    # 0. Set up directory structure
    path_prefix = "flat/" if args.flat else ""
    path_prefix += f"{year}/"
    os.makedirs(downloads_folder(), exist_ok=True)
    os.makedirs(outputs_folder(f"{path_prefix}"), exist_ok=True)
    os.makedirs(outputs_folder(f"{path_prefix}/eia930"), exist_ok=True)
    if not args.skip_outputs:
        # If we are outputing, wipe results dir so we can be confident there are no old
        # result files (eg because of a file name change)
        if os.path.exists(results_folder(f"{path_prefix}")):
            shutil.rmtree(results_folder(f"{path_prefix}"))
        os.makedirs(results_folder(f"{path_prefix}"), exist_ok=False)
    else:
        # still make sure results dir exists, but exist is ok and we won't be writing
        # to it
        os.makedirs(results_folder(f"{path_prefix}"), exist_ok=True)
    os.makedirs(
        results_folder(f"{path_prefix}data_quality_metrics"),
        exist_ok=True,
    )
    # Make results subfolders
    for unit in ["us_units", "metric_units"]:
        for time_resolution in TIME_RESOLUTIONS.keys():
            for subfolder in ["plant_data", "carbon_accounting", "power_sector_data"]:
                os.makedirs(
                    results_folder(
                        f"{path_prefix}/{subfolder}/{time_resolution}/{unit}"
                    ),
                    exist_ok=True,
                )

    # configure the logger
    # Log the print statements to a file for debugging.
    configure_root_logger(
        logfile=results_folder(f"{year}/data_quality_metrics/data_pipeline.log")
    )
    logger = get_logger("data_pipeline")
    print_args(args, logger)

    logger.info(f"Running data pipeline for year {year}")

    # 1. Download data
    ####################################################################################
    logger.info("1. Downloading data")
    logger.info(f"Using {os.getenv('PUDL_BUILD', default='stable')} PUDL build")
    # PUDL
    if PUDL_ENGINE:
        download_data.download_pudl_data(source="aws")
    # eGRID
    download_data.download_egrid_files()
    # EIA-930
    # for `flat` run, we need chalander file for making profiles
    if args.flat:
        download_data.download_chalendar_files()
    # Power Sector Data Crosswalk
    # NOTE: Check for new releases at https://github.com/USEPA/camd-eia-crosswalk
    download_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.3/epa_eia_crosswalk.csv"
    )
    # download the raw EIA-923 and EIA-860 files for use in NOx/SO2 calculations until
    # integrated into pudl
    download_data.download_raw_eia860(year)
    # download raw EIA-860 from the latest validated year for use in subplant
    # identification
    download_data.download_raw_eia860(
        max(latest_validated_year, current_early_release_year)
    )
    download_data.download_raw_eia923(year)

    # 2. Identify subplants
    ####################################################################################
    logger.info("2. Identifying subplant IDs")
    subplant_crosswalk = subplant_identification.generate_subplant_ids()
    output_data.output_intermediate_data(
        subplant_crosswalk,
        "subplant_crosswalk",
        path_prefix,
        year,
        skip_outputs=False,  # always output the crosswalk because it is loaded later
    )
    del subplant_crosswalk

    # 3. Clean EIA-923 Generation and Fuel Data at the Monthly Level
    ####################################################################################
    logger.info("3. Cleaning EIA-923 data")
    (
        eia923_allocated,
        primary_fuel_table,
        subplant_emission_factors,
        subplant_eia923,
    ) = data_cleaning.clean_eia923(year)
    # output primary fuel table
    output_data.output_intermediate_data(
        primary_fuel_table,
        "primary_fuel_table",
        path_prefix,
        year,
        skip_outputs=False,
    )
    # output subplant-level EIA-923 data
    output_data.output_intermediate_data(
        subplant_eia923,
        "subplant_eia923",
        path_prefix,
        year,
        skip_outputs=args.skip_outputs,
    )
    del subplant_eia923
    # Add primary fuel data to each generator
    eia923_allocated = eia923_allocated.merge(
        primary_fuel_table[
            [
                "plant_id_eia",
                "subplant_id",
                "generator_id",
                "energy_source_code",
                "plant_primary_fuel",
                "subplant_primary_fuel",
            ]
        ],
        how="left",
        on=["plant_id_eia", "subplant_id", "generator_id"],
        validate="m:1",
    )

    # 4. Clean Hourly Data from CEMS
    ####################################################################################
    logger.info("4. Cleaning CEMS data")
    cems = data_cleaning.clean_cems(year, primary_fuel_table, subplant_emission_factors)
    del subplant_emission_factors
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

    # remove emissions measurement quality columns and steam load columns to reduce memory
    # NOTE: steam load columns may be used in the future
    cems.drop(
        columns=[
            "steam_load_1000_lb",
            "co2_mass_measurement_code",
            "nox_mass_measurement_code",
            "so2_mass_measurement_code",
        ],
        inplace=True,
    )

    # calculate biomass-adjusted emissions while cems data is at the unit level
    cems = emissions.adjust_emissions_for_biomass(cems)

    # 5. Assign static characteristics to CEMS and EIA data to aid in aggregation
    ####################################################################################
    logger.info("5. Loading plant static attributes")
    plant_attributes = helpers.create_plant_attributes_table(
        cems, eia923_allocated, year, primary_fuel_table
    )
    validation.test_for_missing_values(
        plant_attributes, skip_cols=["plant_retirement_date"]
    )

    # 6. Crosswalk CEMS and EIA data
    ####################################################################################
    logger.info("6. Identifying source for hourly data")
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
        output_data.summarize_annually_reported_eia_data(eia923_allocated, year),
        "annually_reported_eia_data",
        path_prefix,
        args.skip_outputs,
    )

    # 7. Aggregate CEMS data to subplant
    ####################################################################################
    logger.info("7. Aggregating CEMS data from unit to subplant")
    # aggregate cems data to subplant level
    cems = data_cleaning.aggregate_cems_to_subplant(cems)

    # 8. Calculate hourly data for partial_cems plants
    ####################################################################################
    logger.info("8. Shaping partial CEMS data")
    # shape partial CEMS plant data
    partial_cems_plant = impute_hourly_profiles.shape_partial_cems_plants(
        cems, eia923_allocated, year
    )
    validation.validate_unique_datetimes(
        year,
        df=partial_cems_plant,
        df_name="partial_cems_plant",
        keys=["plant_id_eia", "subplant_id"],
    )
    validation.check_for_complete_hourly_timeseries(
        df=partial_cems_plant,
        df_name="partial_cems_plant",
        keys=["plant_id_eia", "subplant_id"],
        period="month",
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
    ) = impute_hourly_profiles.shape_partial_cems_subplants(
        cems, eia923_allocated, year
    )

    validation.validate_unique_datetimes(
        year,
        df=partial_cems_subplant,
        df_name="partial_cems_subplant",
        keys=["plant_id_eia", "subplant_id"],
    )
    validation.check_for_complete_hourly_timeseries(
        df=partial_cems_subplant,
        df_name="partial_cems_subplant",
        keys=["plant_id_eia", "subplant_id"],
        period="month",
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
    logger.info("9. Converting CEMS gross generation to net generation")
    cems, gtn_conversions = gross_to_net_generation.convert_gross_to_net_generation(
        cems, eia923_allocated, primary_fuel_table, year
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
    del gtn_conversions

    # 10. Adjust CEMS emission data for CHP
    ####################################################################################
    logger.info("10. Adjusting CEMS emissions for CHP")
    cems = data_cleaning.adjust_cems_for_chp(cems, eia923_allocated)
    cems = emissions.calculate_co2e_mass(
        cems, year, gwp_horizon=100, ar5_climate_carbon_feedback=True
    )
    validation.test_emissions_adjustments(cems, year)
    validation.validate_unique_datetimes(
        year,
        df=cems,
        df_name="cems_subplant",
        keys=["plant_id_eia", "subplant_id"],
    )
    validation.check_for_complete_hourly_timeseries(
        df=cems,
        df_name="cems_subplant",
        keys=["plant_id_eia", "subplant_id"],
        period="month",
    )
    # now that GTN calculations are complete, we can remove zero observations from the
    # cems data to free up memory
    # dropping here means cems_subplant will contain missing timestamps
    cems = data_cleaning.remove_cems_with_zero_monthly_data(
        cems,
        ["gross_generation_mwh", "net_generation_mwh", "fuel_consumed_mmbtu"],
        remove_all_zeros=True,
    )
    partial_cems_plant = data_cleaning.remove_cems_with_zero_monthly_data(
        partial_cems_plant,
        ["net_generation_mwh", "fuel_consumed_mmbtu"],
        remove_all_zeros=True,
    )
    partial_cems_subplant = data_cleaning.remove_cems_with_zero_monthly_data(
        partial_cems_subplant,
        ["net_generation_mwh", "fuel_consumed_mmbtu"],
        remove_all_zeros=True,
    )
    output_data.output_intermediate_data(
        cems, "cems_subplant", path_prefix, year, args.skip_outputs
    )

    # 11. Export monthly and annual plant-level results
    ####################################################################################
    logger.info("11. Exporting monthly and annual plant-level results")
    # create a separate dataframe containing only the EIA data that is missing from cems
    monthly_eia_data_to_shape = eia923_allocated[
        (eia923_allocated["hourly_data_source"] == "eia")
    ]
    output_data.output_data_quality_metrics(
        output_data.identify_percent_of_data_by_input_source(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            monthly_eia_data_to_shape,
            year,
            plant_attributes,
            primary_fuel_table,
        ),
        "input_data_source",
        path_prefix,
        args.skip_outputs,
    )
    # group EIA data to subplant level
    monthly_eia_data_to_shape = (
        monthly_eia_data_to_shape.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )[DATA_COLUMNS]
        .sum()
        .reset_index()
    )
    # combine and export plant data at monthly and annual level
    validation.ensure_non_overlapping_data_from_all_sources(
        cems, partial_cems_subplant, partial_cems_plant, monthly_eia_data_to_shape
    )
    monthly_subplant_data = helpers.combine_subplant_data(
        cems,
        partial_cems_subplant,
        partial_cems_plant,
        monthly_eia_data_to_shape,
        resolution="monthly",
    )

    # export subplant attributes table
    helpers.create_subplant_attributes_table(
        monthly_subplant_data, plant_attributes, primary_fuel_table, year, path_prefix
    )
    # For years before 2019, export the plant attributes table now
    if year < earliest_hourly_data_year:
        # export plant static attributes to csv
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

    validation.check_for_complete_monthly_timeseries(
        df=monthly_subplant_data,
        df_name="monthly_plant_data",
        keys=["plant_id_eia", "subplant_id"],
        columns_to_check=["net_generation_mwh", "fuel_consumed_for_electricity_mmbtu"],
        year=year,
    )
    # output plant and subplant data
    for plant_part in ["subplant", "plant"]:
        for resolution in ["monthly", "annual"]:
            output_data.write_plant_data_to_results(
                monthly_subplant_data,
                year,
                path_prefix=path_prefix,
                plant_part=plant_part,
                resolution=resolution,
                skip_outputs=args.skip_outputs,
                plant_attributes=plant_attributes,
            )
    # 12. Export monthly and annual power sector data
    ####################################################################################
    logger.info("12. Exporting monthly and annual fleet-level results")
    # output monthly/annual power sector results
    fleet_data = data_cleaning.aggregate_subplant_data_to_fleet(
        monthly_subplant_data, plant_attributes, primary_fuel_table, year
    )
    output_data.write_power_sector_results(
        fleet_data,
        year,
        path_prefix,
        args.skip_outputs,
        include_hourly=False,
        include_monthly=True,
        include_annual=True,
    )
    # free up memory
    del monthly_subplant_data
    del fleet_data
    gc.collect()

    # calculate hourly outputs for years after 2019
    if year >= earliest_hourly_data_year:
        # 13. Clean and Reconcile EIA-930 data
        ################################################################################
        logger.info("13. Cleaning EIA-930 data")
        # Scrapes and cleans data in data/downloads, outputs cleaned file at
        # EBA_elec.csv
        if args.flat:
            logger.info(
                "Not running 930 cleaning because we'll be using a flat profile."
            )
        elif not (
            os.path.exists(outputs_folder(f"{path_prefix}/eia930/eia930_elec.csv"))
        ):
            eia930.clean_930(year, path_prefix=path_prefix)
        else:
            logger.info(
                "Not re-running 930 cleaning. If you'd like to re-run, "
                f"please delete data/outputs/{path_prefix}/eia930/"
            )

        # If running flat, we need to use the Chalender file to build residual profiles.
        clean_930_file = (
            downloads_folder("eia930/chalendar/EBA_elec.csv")
            if args.flat
            else outputs_folder(f"{path_prefix}/eia930/eia930_elec.csv")
        )
        eia930_data = eia930.load_chalendar_for_pipeline(clean_930_file, year=year)
        # until we can fix the physics reconciliation, we need to apply some
        # post-processing steps
        eia930_data = eia930.remove_imputed_ones(eia930_data)

        # 14. Calculate hourly profiles for monthly EIA data
        ################################################################################
        logger.info("14. Estimating hourly profiles for EIA data")
        hourly_profiles = impute_hourly_profiles.calculate_hourly_profiles(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            eia930_data,
            plant_attributes,
            primary_fuel_table,
            monthly_eia_data_to_shape,
            year,
            transmission_only=False,
            ba_column_name="ba_code",
            use_flat=args.flat,
        )
        del eia930_data
        gc.collect()
        # validate how well the wind and solar imputation methods work
        output_data.output_data_quality_metrics(
            impute_hourly_profiles.validate_wind_solar_imputation(
                hourly_profiles, year
            ),
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

        # 15. Export hourly plant-level data
        ################################################################################
        # This provides a complete hourly timeseries for each plant, but is more
        # memory-intensive and results in larger output files.
        # This data is immediately exported and not held in memory for the rest of the
        # pipeline. Instead, this data is re-combined again using the aggregated shaped
        # plants in step 17.
        logger.info("15. Exporting Hourly Plant-level data for each BA")
        impute_hourly_profiles.combine_and_export_hourly_plant_data(
            year,
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            monthly_eia_data_to_shape,
            plant_attributes,
            primary_fuel_table,
            hourly_profiles,
            path_prefix,
            args.skip_outputs,
            region_to_group="ba_code",
        )

        # 16. Shape fleet-level data
        ################################################################################
        logger.info("16. Assigning hourly profiles to monthly EIA-923 data")
        # Aggregate EIA data to BA/fuel/month, then assign hourly profile per BA/fuel
        monthly_eia_fleet_data = impute_hourly_profiles.aggregate_eia_data_to_fleet(
            monthly_eia_data_to_shape, plant_attributes, primary_fuel_table, year
        )
        shaped_eia_fleet_data = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
            monthly_eia_fleet_data,
            hourly_profiles,
            year,
            fuel_category_col_for_shaping="fuel_category_for_shaping",
        )
        del hourly_profiles
        gc.collect()

        output_data.output_data_quality_metrics(
            validation.hourly_profile_source_metric(
                cems,
                partial_cems_subplant,
                partial_cems_plant,
                shaped_eia_fleet_data,
                plant_attributes,
            ),
            "hourly_profile_method",
            path_prefix,
            args.skip_outputs,
        )
        # Export data
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
            shaped_eia_fleet_data,
            monthly_eia_fleet_data,
            year,
            group_keys=["ba_code", "fuel_category"],
        )
        del monthly_eia_fleet_data

        # write metadata outputs
        output_data.write_plant_metadata(
            plant_attributes,
            eia923_allocated,
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            shaped_eia_fleet_data,
            path_prefix,
            args.skip_outputs,
            year,
        )
        del eia923_allocated
        gc.collect()
        # group shaped data by fleet, since the fuel category used for shaping might
        # not match the fuel category for fleet aggregation
        shaped_eia_fleet_data = (
            shaped_eia_fleet_data.groupby(
                [
                    "ba_code",
                    "fuel_category",
                    "datetime_utc",
                    "report_date",
                ],
                dropna=False,
                sort=False,
            )[DATA_COLUMNS]
            .sum(numeric_only=True)
            .reset_index()
        )

        # 17. Combine plant-level data from all sources
        ################################################################################
        logger.info("17. Combining hourly CEMS data")
        # Because EIA data is already aggregated to the fleet level, at this point, we
        # only want to combine CEMS data together. We pass a blank dataframe as the
        # EIA data
        validation.ensure_non_overlapping_data_from_all_sources(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            eia_data=pd.DataFrame(
                columns=["plant_id_eia", "subplant_id", "report_date", "datetime_utc"]
            ),
        )
        combined_cems_subplant_data = helpers.combine_subplant_data(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            eia_data=pd.DataFrame(
                columns=["plant_id_eia", "subplant_id", "report_date", "datetime_utc"]
            ),
            resolution="hourly",
        )
        # free memory back to python
        del (
            cems,
            partial_cems_subplant,
            partial_cems_plant,
        )
        gc.collect()
        # ensure complete timeseries
        combined_cems_subplant_data = data_cleaning.complete_hourly_timeseries(
            combined_cems_subplant_data,
            year,
            group_cols=["plant_id_eia", "subplant_id"],
            columns_to_fill_with_zero=DATA_COLUMNS,
        )
        validation.validate_unique_datetimes(
            year,
            df=combined_cems_subplant_data,
            df_name="combined_plant_data",
            keys=["plant_id_eia", "subplant_id"],
        )

        # 18. Aggregate CEMS data to fleet and write power sector results
        ################################################################################
        logger.info(
            "18. Creating and exporting hourly, fleet-level power sector results"
        )
        # aggregate CEMS data to the fleet level
        cems_fleet_data = data_cleaning.aggregate_subplant_data_to_fleet(
            combined_cems_subplant_data, plant_attributes, primary_fuel_table, year
        )
        del combined_cems_subplant_data
        del plant_attributes
        del primary_fuel_table
        gc.collect()

        # combine fleet-level CEMS data and EIA data into a single df and group
        # fleets together
        combined_fleet_data = pd.concat(
            [cems_fleet_data, shaped_eia_fleet_data], axis=0, copy=False
        )
        del cems_fleet_data
        del shaped_eia_fleet_data
        gc.collect()
        combined_fleet_data = (
            combined_fleet_data.groupby(
                ["ba_code", "fuel_category", "datetime_utc", "report_date"],
                dropna=False,
            )[data_cleaning.DATA_COLUMNS]
            .sum()
            .reset_index()
        )

        # Output final data: per-ba hourly generation and rate
        output_data.write_power_sector_results(
            combined_fleet_data,
            year,
            path_prefix,
            args.skip_outputs,
            include_hourly=True,
            include_monthly=False,
            include_annual=False,
        )
        # Write US-average fleet data
        output_data.write_national_fleet_averages(
            combined_fleet_data, year, path_prefix, skip_outputs=False
        )

        # 19. Calculate consumption-based emissions and write carbon accounting results
        ################################################################################
        logger.info("19. Calculating and exporting consumption-based results")
        hourly_consumed_calc = consumed.HourlyConsumed(
            clean_930_file,
            path_prefix,
            year,
            skip_outputs=args.skip_outputs,
        )
        hourly_consumed_calc.run()
        hourly_consumed_calc.output_results()


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
