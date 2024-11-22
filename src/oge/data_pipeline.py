"""
Entry point for creating final dataset and intermediate cleaned data products.

Run from `src` as `python data_pipeline.py` after installing conda environment

Optional arguments are --year (default 2022)
Optional arguments for development are --small, --flat, and --skip_outputs
"""

import argparse
import os
import shutil

import pandas as pd

# import local modules
import oge.download_data as download_data
import oge.data_cleaning as data_cleaning
import oge.subplant_identification as subplant_identification
import oge.emissions as emissions
import oge.gross_to_net_generation as gross_to_net_generation
import oge.helpers as helpers
import oge.impute_hourly_profiles as impute_hourly_profiles
import oge.eia930 as eia930
import oge.validation as validation
import oge.output_data as output_data
import oge.consumed as consumed
from oge.filepaths import downloads_folder, outputs_folder, results_folder
from oge.logging_util import get_logger, configure_root_logger
from oge.constants import (
    TIME_RESOLUTIONS,
    latest_validated_year,
    current_early_release_year,
    earliest_hourly_data_year,
)


def get_args() -> argparse.Namespace:
    """Specify arguments here.

    Returns dictionary of {arg_name: arg_value}
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", help="Year for analysis", default=2022, type=int)
    parser.add_argument(
        "--small",
        help=(
            "Run on subset of data for quicker testing, outputs to outputs/small and "
            "results to results/small."
        ),
        default=False,
        action=argparse.BooleanOptionalAction,
    )
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
    path_prefix = "" if not args.small else "small/"
    path_prefix += "flat/" if args.flat else ""
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
    # PUDL
    download_data.download_pudl_data(source="aws")
    logger.info(f"Using {os.getenv('PUDL_BUILD', default='stable')} PUDL build")
    # eGRID
    download_data.download_egrid_files()
    # EIA-930
    # for `small` run, we'll only clean 1 week, so need chalander file for making
    # profiles
    if args.small or args.flat:
        download_data.download_chalendar_files()
    # We use balance files for imputing missing hourly profiles.
    # need last year for rolling data cleaning
    download_data.download_eia930_data(years_to_download=[year, year - 1])
    # Power Sector Data Crosswalk
    # NOTE: Check for new releases at https://github.com/USEPA/camd-eia-crosswalk
    download_data.download_epa_psdc(
        psdc_url="https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.3/epa_eia_crosswalk.csv"
    )
    # download the raw EIA-923 and EIA-860 files for use in NOx/SO2 calculations until
    # integrated into pudl
    download_data.download_raw_eia860(year)
    # download eia860 from the latest validated year for use in subplant identification
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
        args.skip_outputs,
    )

    # 3. Clean EIA-923 Generation and Fuel Data at the Monthly Level
    ####################################################################################
    logger.info("3. Cleaning EIA-923 data")
    (
        eia923_allocated,
        primary_fuel_table,
        subplant_emission_factors,
    ) = data_cleaning.clean_eia923(year, args.small)
    # output primary fuel table
    output_data.output_intermediate_data(
        primary_fuel_table,
        "primary_fuel_table",
        path_prefix,
        year,
        args.skip_outputs,
    )
    # remove intermediate columns from primary fuel table
    primary_fuel_table = primary_fuel_table[
        [
            "plant_id_eia",
            "subplant_id",
            "generator_id",
            "energy_source_code",
            "plant_primary_fuel",
            "subplant_primary_fuel",
        ]
    ]
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
    logger.info("5. Loading plant static attributes")
    plant_attributes = helpers.create_plant_attributes_table(
        cems, eia923_allocated, year, primary_fuel_table
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
        validation.summarize_annually_reported_eia_data(eia923_allocated, year),
        "annually_reported_eia_data",
        path_prefix,
        args.skip_outputs,
    )

    # 7. Aggregating CEMS data to subplant
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
        validation.identify_percent_of_data_by_input_source(
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
    # combine and export plant data at monthly and annual level
    monthly_subplant_data = data_cleaning.combine_subplant_data(
        cems,
        partial_cems_subplant,
        partial_cems_plant,
        monthly_eia_data_to_shape,
        resolution="monthly",
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
                path_prefix,
                resolution,
                plant_part,
                args.skip_outputs,
                plant_attributes,
            )
    # 12. Export monthly and annual power sector data
    ####################################################################################
    logger.info("12. Exporting monthly and annual fleet-level results")
    # output monthly/annual power sector results
    fleet_data = data_cleaning.aggregate_subplant_data_to_fleet(
        monthly_subplant_data, plant_attributes, primary_fuel_table
    )
    output_data.write_power_sector_results(
        fleet_data, year, path_prefix, args.skip_outputs, include_hourly=False
    )

    # TODO: does fleet_data need to be kept around?

    # For 2019 onward, calculate hourly data, otherwise skip these steps
    if year >= earliest_hourly_data_year:
        del monthly_subplant_data
        # 12. Clean and Reconcile EIA-930 data
        ################################################################################
        logger.info("12. Cleaning EIA-930 data")
        # Scrapes and cleans data in data/downloads, outputs cleaned file at
        # EBA_elec.csv
        if args.flat:
            logger.info(
                "Not running 930 cleaning because we'll be using a flat profile."
            )
        elif not (
            os.path.exists(outputs_folder(f"{path_prefix}/eia930/eia930_elec.csv"))
        ):
            eia930.clean_930(year, small=args.small, path_prefix=path_prefix)
        else:
            logger.info(
                "Not re-running 930 cleaning. If you'd like to re-run, "
                f"please delete data/outputs/{path_prefix}/eia930/"
            )

        # If running small, we didn't clean the whole year, so need to use the
        # Chalender file to build residual profiles.
        clean_930_file = (
            downloads_folder("eia930/chalendar/EBA_elec.csv")
            if (args.small or args.flat)
            else outputs_folder(f"{path_prefix}/eia930/eia930_elec.csv")
        )
        eia930_data = eia930.load_chalendar_for_pipeline(clean_930_file, year=year)
        # until we can fix the physics reconciliation, we need to apply some
        # post-processing steps
        eia930_data = eia930.remove_imputed_ones(eia930_data)
        eia930_data = eia930.remove_months_with_zero_data(eia930_data)

        # 13. Calculate hourly profiles for monthly EIA data
        ################################################################################
        logger.info("13. Estimating hourly profiles for EIA data")
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
        ################################################################################
        # The data pipeline offers two options for exporting plant-level data. The
        # first option, executed here in step 14 is to shape hourly data for each
        # individual plant.
        # This option provides a complete hourly timeseries for each plant, but is more
        # memory-intensive and results in larger output files.
        # This data is immediately exported and not held in memory for the rest of the
        # pipeline. Instead, this data is re-combined again using the aggregated shaped
        # plants in step 16.
        # The other option happens in step 16 if step 14 is not run: instead of shaping
        # each plant, only an aggregate fleet-level synthetic plant is exported, which
        # represents the characteristics of all plants that do not report hourly data
        # elsewhere in each month.
        logger.info("14. Exporting Hourly Plant-level data for each BA")
        if not args.small:
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

        # 15. Shape fleet-level data
        ################################################################################
        logger.info("15. Assigning hourly profiles to monthly EIA-923 data")
        hourly_profiles = impute_hourly_profiles.convert_profile_to_percent(
            hourly_profiles,
            group_keys=["ba_code", "fuel_category", "profile_method"],
            columns_to_convert=["profile", "flat_profile"],
        )
        # Aggregate EIA data to BA/fuel/month, then assign hourly profile per BA/fuel
        (
            monthly_eia_fleet_data,
            plant_attributes,
        ) = impute_hourly_profiles.aggregate_eia_data_to_fleet(
            monthly_eia_data_to_shape, plant_attributes, primary_fuel_table
        )
        shaped_eia_fleet_data = impute_hourly_profiles.shape_monthly_eia_data_as_hourly(
            monthly_eia_fleet_data, hourly_profiles, year
        )
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

        # 16. Combine plant-level data from all sources
        ################################################################################
        logger.info("16. Combining plant-level hourly data")
        # write metadata outputs
        # TODO: remove shaped plant id from this and use ba-fuel instead
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

        # Because EIA data is already aggregated to the fleet level, at this point, we
        # only want to combine CEMS data together. We pass a blank dataframe as the
        # EIA data
        combined_cems_subplant_data = data_cleaning.combine_subplant_data(
            cems,
            partial_cems_subplant,
            partial_cems_plant,
            eia_data=pd.DataFrame(
                columns=["plant_id_eia", "subplant_id", "report_date", "datetime_utc"]
            ),
            resolution="hourly",
            validate=True,
        )
        # free memory back to python
        del (
            cems,
            partial_cems_subplant,
            partial_cems_plant,
        )
        # export to a csv.
        validation.validate_unique_datetimes(
            year,
            df=combined_cems_subplant_data,
            df_name="combined_plant_data",
            keys=["plant_id_eia", "subplant_id"],
        )

        # 17. Aggregate CEMS data to BA-fuel and write power sector results
        ################################################################################
        logger.info("17. Creating and exporting BA-level power sector results")
        # aggregate CEMS data to the fleet level
        cems_fleet_data = data_cleaning.aggregate_subplant_data_to_fleet(
            combined_cems_subplant_data, plant_attributes, primary_fuel_table
        )
        del combined_cems_subplant_data

        # combine fleet-level CEMS data and EIA data into a single df and group
        # fleets together
        combined_fleet_data = pd.concat(
            [cems_fleet_data, shaped_eia_fleet_data], axis=0
        )
        combined_fleet_data = (
            combined_fleet_data.groupby(
                ["ba_code", "fuel_category", "datetime_utc", "report_date"],
                dropna=False,
            )[data_cleaning.DATA_COLUMNS]
            .sum()
            .reset_index()
        )

        # Output intermediate data: produced per-fuel annual averages
        output_data.write_generated_averages(
            combined_fleet_data, year, path_prefix, args.skip_outputs
        )
        # Output final data: per-ba hourly generation and rate
        output_data.write_power_sector_results(
            combined_fleet_data,
            year,
            path_prefix,
            args.skip_outputs,
            include_hourly=True,
        )

        # 18. Calculate consumption-based emissions and write carbon accounting results
        ################################################################################
        logger.info("18. Calculating and exporting consumption-based results")
        hourly_consumed_calc = consumed.HourlyConsumed(
            clean_930_file,
            path_prefix,
            year,
            small=args.small,
            skip_outputs=args.skip_outputs,
        )
        hourly_consumed_calc.run()
        hourly_consumed_calc.output_results()

    # for years prior to 2019, do not export carbon accounting data
    elif year < earliest_hourly_data_year:
        # export plant static attributes to csv
        output_data.output_intermediate_data(
            plant_attributes.assign(shaped_plant_id=pd.NA),
            "plant_static_attributes",
            path_prefix,
            year,
            args.skip_outputs,
        )
        if not args.skip_outputs:
            plant_attributes.assign(shaped_plant_id=pd.NA).to_csv(
                results_folder(f"{path_prefix}plant_data/plant_static_attributes.csv"),
                index=False,
            )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
