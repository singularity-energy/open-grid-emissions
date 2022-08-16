import pandas as pd
import numpy as np

import load_data
import impute_hourly_profiles
from column_checks import get_dtypes
from filepaths import *


# DATA PIPELINE VALIDATION FUNCTIONS
########################################################################################


def validate_year(year):
    """Returns a warning if the year specified is not known to work with the pipeline."""

    earliest_validated_year = 2019
    latest_validated_year = 2020

    if year < earliest_validated_year:
        year_warning = f"""
        ################################################################################
        Warning: The data pipeline has only been validated to work for years {earliest_validated_year}-{latest_validated_year}.
        Running the pipeline for {year} may cause it to fail or may lead to poor-quality
        or anomalous results. To check on the progress of validating additional years of
        data, see: https://github.com/singularity-energy/open-grid-emissions/issues/117
        ################################################################################
        """
        print(year_warning)
    elif year > latest_validated_year:
        year_warning = f"""
        ################################################################################
        Warning: The most recent available year of input data is currently {latest_validated_year}.
        Input data for {year} should be available from the EIA in Fall {year+1} and we will
        work to validate that the pipeline works with {year} data as soon as possible
        after the data is released.
        ################################################################################
        """
        raise UserWarning(year_warning)


def test_for_negative_values(df):
    print("    Checking that fuel and emissions values are positive...  ", end="")
    columns_that_should_be_positive = [
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
    ]
    columns_to_test = [
        col for col in columns_that_should_be_positive if col in df.columns
    ]
    negative_warnings = 0
    for column in columns_to_test:
        negative_test = df[df[column] < 0]
        if not negative_test.empty:
            print(" ")
            print(
                f"Warning: There are {len(negative_test)} records where {column} is negative. Check `negative_test` for complete list"
            )
            print(" ")
            negative_warnings += 1
    if negative_warnings > 0:
        raise UserWarning("The above negative values are errors and must be fixed")
    else:
        print("OK")
    return negative_test


def test_chp_allocation(df):
    print(
        "    Checking that total fuel consumed >= fuel consumed for electricity...  ",
        end="",
    )
    chp_allocation_test = df[
        df["fuel_consumed_for_electricity_mmbtu"] > df["fuel_consumed_mmbtu"]
    ]
    if not chp_allocation_test.empty:
        raise UserWarning(
            f"Warning: There are {len(chp_allocation_test)} records where fuel consumed for electricity is greater than total fuel consumption. Check `chp_allocation_test` for complete list"
        )
    else:
        print("OK")

    return chp_allocation_test


def test_for_missing_energy_source_code(df):
    print(
        "    Checking that there are no missing energy source codes associated with non-zero fuel consumption...  ",
        end="",
    )
    missing_esc_test = df[
        (df["energy_source_code"].isna()) & (df["fuel_consumed_mmbtu"] > 0)
    ]
    if not missing_esc_test.empty:
        print(" ")
        print(
            f"Warning: There are {len(missing_esc_test)} records where there is a missing energy source code associated with non-zero fuel consumption. Check `missing_esc_test` for complete list"
        )
    else:
        print("OK")

    return missing_esc_test


def test_for_missing_subplant_id(df):
    print("    Checking that all data has an associated `subplant_id`...  ", end="")
    missing_subplant_test = df[df["subplant_id"].isna()]
    if not missing_subplant_test.empty:
        print(" ")
        print(
            f"Warning: There are {len(missing_subplant_test)} records for {len(missing_subplant_test[['plant_id_eia']].drop_duplicates())} plants without a subplant ID. See `missing_subplant_test` for details"
        )
    else:
        print("OK")
    return missing_subplant_test


def validate_gross_to_net_conversion(cems, eia923_allocated):
    "checks whether the calculated net generation matches the reported net generation from EIA-923 at the annual plant level."
    print(
        "    Checking that calculated net generation matches reported net generation in EIA-923...  ",
        end="",
    )
    # merge together monthly subplant totals from EIA and calculated from CEMS
    eia_netgen = (
        eia923_allocated.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .sum(min_count=1)["net_generation_mwh"]
        .reset_index()
        .dropna(subset="net_generation_mwh")
    )
    calculated_netgen = (
        cems.groupby(["plant_id_eia", "subplant_id", "report_date"], dropna=False)
        .sum()["net_generation_mwh"]
        .reset_index()
    )
    validated_ng = eia_netgen.merge(
        calculated_netgen,
        how="inner",
        on=["plant_id_eia", "subplant_id", "report_date"],
        suffixes=("_eia", "_calc"),
        validate="1:1",
    )

    validated_ng = validated_ng.groupby("plant_id_eia").sum()[
        ["net_generation_mwh_eia", "net_generation_mwh_calc"]
    ]

    validated_ng = validated_ng.round(3)
    validated_ng = validated_ng[
        validated_ng[["net_generation_mwh_eia", "net_generation_mwh_calc"]].sum(axis=1)
        != 0
    ]

    validated_ng["pct_error"] = (
        validated_ng["net_generation_mwh_calc"] - validated_ng["net_generation_mwh_eia"]
    ) / validated_ng["net_generation_mwh_eia"]

    cems_net_not_equal_to_eia = validated_ng[validated_ng["pct_error"] != 0]

    if len(cems_net_not_equal_to_eia) > 0:
        print(" ")
        print(
            f"Warning: There are {len(cems_net_not_equal_to_eia)} plants where calculated annual net generation does not match EIA annual net generation."
        )
        print(cems_net_not_equal_to_eia)
    else:
        print("OK")


def test_emissions_adjustments(df):
    """For each emission, tests that mass_lb >= mass_lb_for_electricity >= mass_lb_for_electricity_adjusted."""

    print(
        "    Checking that adjusted emission values are less than total emissions...  ",
        end="",
    )

    pollutants = ["co2", "ch4", "n2o", "co2e", "nox", "so2"]

    bad_adjustments = 0

    for pollutant in pollutants:
        # test that mass_lb >= mass_lb_for_electricity
        bad_adjustment = df[
            (df[f"{pollutant}_mass_lb"] < df[f"{pollutant}_mass_lb_for_electricity"])
        ]
        if len(bad_adjustment) > 0:
            print(
                f"Warning: There are {len(bad_adjustment)} records where {pollutant}_mass_lb_for_electricity > {pollutant}_mass_lb"
            )
            bad_adjustment += 1

        # test that mass_lb >= mass_lb_adjusted
        bad_adjustment = df[
            (df[f"{pollutant}_mass_lb"] < df[f"{pollutant}_mass_lb_adjusted"])
        ]
        if len(bad_adjustment) > 0:
            print(
                f"Warning: There are {len(bad_adjustment)} records where {pollutant}_mass_lb_adjusted > {pollutant}_mass_lb"
            )
            bad_adjustment += 1

        # test that mass_lb_for_electricity >= mass_lb_for_electricity_adjusted
        bad_adjustment = df[
            (
                df[f"{pollutant}_mass_lb_for_electricity"]
                < df[f"{pollutant}_mass_lb_for_electricity_adjusted"]
            )
        ]
        if len(bad_adjustment) > 0:
            print(" ")
            print(
                f"Warning: There are {len(bad_adjustment)} records where {pollutant}_mass_lb_for_electricity_adjusted > {pollutant}_mass_lb_for_electricity"
            )
            bad_adjustment += 1

    # if there were any bad adjustments, raise a userwarning.
    if bad_adjustments > 0:
        raise UserWarning("The above issues with emissions adjustments must be fixed.")
    else:
        print("OK")


def ensure_non_overlapping_data_from_all_sources(cems, partial_cems, eia_data):

    print("    Checking that all data to be combined is unique...  ", end="")

    if "hourly_data_source" in eia_data.columns:
        eia_only_data = eia_data.loc[
            eia_data["hourly_data_source"] == "eia",
            ["plant_id_eia", "subplant_id", "report_date"],
        ].drop_duplicates()
    else:
        eia_only_data = eia_data[
            ["plant_id_eia", "subplant_id", "report_date"]
        ].drop_duplicates()
    eia_only_data["in_eia"] = 1

    cems_data = cems[["plant_id_eia", "subplant_id", "report_date"]].drop_duplicates()
    cems_data["in_cems"] = 1

    partial_cems_data = partial_cems[
        ["plant_id_eia", "subplant_id", "report_date"]
    ].drop_duplicates()
    partial_cems_data["in_partial_cems"] = 1

    data_overlap = eia_only_data.merge(
        cems_data,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="1:1",
    )
    data_overlap = data_overlap.merge(
        partial_cems_data,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="1:1",
    )
    data_overlap[["in_eia", "in_cems", "in_partial_cems"]] = data_overlap[
        ["in_eia", "in_cems", "in_partial_cems"]
    ].fillna(0)
    data_overlap["number_of_locations"] = (
        data_overlap["in_eia"]
        + data_overlap["in_cems"]
        + data_overlap["in_partial_cems"]
    )

    if len(data_overlap[data_overlap["number_of_locations"] > 1]) > 0:
        eia_cems_overlap = data_overlap[
            (data_overlap["in_eia"] == 1) & (data_overlap["in_cems"] == 1)
        ]
        if len(eia_cems_overlap) > 0:
            print(" ")
            print(
                f"Warning: There are {len(eia_cems_overlap)} subplant-months that exist in both shaped EIA data and CEMS"
            )
        eia_pc_overlap = data_overlap[
            (data_overlap["in_eia"] == 1) & (data_overlap["in_partial_cems"] == 1)
        ]
        if len(eia_pc_overlap) > 0:
            print(" ")
            print(
                f"Warning: There are {len(eia_pc_overlap)} subplant-months that exist in both shaped EIA data and partial CEMS data"
            )
        cems_pc_overlap = data_overlap[
            (data_overlap["in_cems"] == 1) & (data_overlap["in_partial_cems"] == 1)
        ]
        if len(cems_pc_overlap) > 0:
            print(" ")
            print(
                f"Warning: There are {len(cems_pc_overlap)} subplant-months that exist in both CEMS data and partial CEMS data"
            )
        all_overlap = data_overlap[data_overlap["number_of_locations"] == 3]
        if len(all_overlap) > 0:
            print(" ")
            print(
                f"Warning: There are {len(all_overlap)} subplant-months that exist in shaped EIA data, CEMS data, and partial CEMS data."
            )
        raise UserWarning("The above overlaps must be fixed before proceeding.")
    else:
        print("OK")


def validate_shaped_totals(shaped_eia_data, monthly_eia_data_to_shape):

    print("    Checking that shaped hourly data matches monthly totals...  ", end="")

    # aggregate data to ba fuel month
    shaped_data_agg = shaped_eia_data.groupby(
        ["ba_code", "fuel_category", "report_date"], dropna=False
    ).sum()[["net_generation_mwh", "fuel_consumed_mmbtu"]]
    eia_data_agg = monthly_eia_data_to_shape.groupby(
        ["ba_code", "fuel_category", "report_date"], dropna=False
    ).sum()[["net_generation_mwh", "fuel_consumed_mmbtu"]]

    # calculate the difference between the two datasets
    compare = (shaped_data_agg - eia_data_agg).round(0)

    if compare.sum().sum() > 0:
        print(" ")
        print(
            compare[
                (compare["net_generation_mwh"] != 0)
                | (compare["fuel_consumed_mmbtu"] != 0)
            ]
        )
        raise UserWarning(
            "The EIA process is changing the monthly total values compared to reported EIA values. This process should only shape the data, not alter it."
        )
    else:
        print("OK")


# DATA QUALITY METRIC FUNCTIONS
########################################################################################


def hourly_profile_source_metric(cems, partial_cems, shaped_eia_data):
    """Calculates the percentage of data whose hourly profile was determined by method"""
    data_metrics = ["net_generation_mwh", "co2_mass_lb"]

    # determine the source of the hourly profile
    profile_from_cems = pd.DataFrame(cems[data_metrics].sum(axis=0)).T
    profile_from_cems["profile_method"] = "cems_reported"

    profile_from_partial_cems = pd.DataFrame(partial_cems[data_metrics].sum(axis=0)).T
    profile_from_partial_cems["profile_method"] = "eia_shaped_partial_cems"

    profile_from_eia = (
        shaped_eia_data.groupby("profile_method", dropna=False)[data_metrics]
        .sum()
        .reset_index()
    )
    profile_from_eia["profile_method"] = (
        "eia_shaped_" + profile_from_eia["profile_method"]
    )

    profile_source = pd.concat(
        [profile_from_cems, profile_from_partial_cems, profile_from_eia]
    )
    profile_source = profile_source.set_index("profile_method")
    profile_source = profile_source / profile_source.sum(axis=0)
    profile_source = profile_source.reset_index()

    return profile_source


def identify_percent_of_data_by_input_source(cems, partial_cems, eia_only_data, year):
    data_sources = {
        "cems": cems,
        "partial_cems": partial_cems,
        "eia": eia_only_data,
    }
    if year % 4 == 0:
        hours_in_year = 8784
    else:
        hours_in_year = 8760
    source_of_input_data = []
    for name, df in data_sources.items():
        if len(df) == 0:  # Empty df. May occur when running `small`
            print(f"WARNING: data source {name} has zero entries")
            continue
        if name == "eia":
            subplant_data = df.groupby(
                ["plant_id_eia", "subplant_id"], dropna=False
            ).sum()[
                ["net_generation_mwh", "co2_mass_lb", "co2_mass_lb_for_electricity"]
            ]
            subplant_hours = len(subplant_data) * hours_in_year
        else:
            subplant_data = df.groupby(
                ["plant_id_eia", "subplant_id", "datetime_utc"], dropna=False
            ).sum()[
                ["net_generation_mwh", "co2_mass_lb", "co2_mass_lb_for_electricity"]
            ]
            subplant_hours = len(subplant_data)
        summary = pd.DataFrame.from_dict(
            {
                "source": [name],
                "subplant_hours": [subplant_hours],
                "net_generation_mwh": [subplant_data["net_generation_mwh"].sum()],
                "co2_mass_lb": [subplant_data["co2_mass_lb"].sum()],
                "co2_mass_lb_for_electricity": [
                    subplant_data["co2_mass_lb_for_electricity"].sum()
                ],
            }
        )
        source_of_input_data.append(summary)
    source_of_input_data = pd.concat(source_of_input_data)

    source_of_input_data["source"] = source_of_input_data["source"].replace(
        "partial_cems", "eia"
    )
    source_of_input_data = source_of_input_data.groupby("source").sum()
    source_of_input_data = source_of_input_data / source_of_input_data.sum(axis=0)

    source_of_input_data = source_of_input_data.reset_index()

    return source_of_input_data


def identify_cems_gtn_method(cems):
    method_summary = cems.groupby("gtn_method", dropna=False).sum()[
        "gross_generation_mwh"
    ]
    method_summary = method_summary / method_summary.sum(axis=0)
    method_summary = method_summary.reset_index()
    method_summary["gtn_method"] = method_summary["gtn_method"].astype(str)
    method_summary = method_summary.sort_values(by="gtn_method", axis=0)
    return method_summary


def validate_diba_imputation_method(hourly_profiles, year):

    # only keep wind and solar data
    data_to_validate = hourly_profiles[
        (hourly_profiles["fuel_category"].isin(["wind", "solar"]))
        & (~hourly_profiles["eia930_profile"].isna())
    ]
    data_to_validate = data_to_validate[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "eia930_profile",
        ]
    ]

    profiles_to_impute = data_to_validate[
        ["ba_code", "fuel_category", "report_date"]
    ].drop_duplicates()

    profiles_to_impute = profiles_to_impute[
        profiles_to_impute["report_date"].dt.year == 2020
    ]

    dibas = load_data.load_diba_data(year)

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    for index, row in profiles_to_impute.iterrows():
        ba = row["ba_code"]
        fuel = row["fuel_category"]
        report_date = row["report_date"]

        # for wind and solar, average the wind and solar generation profiles from
        # nearby interconnected BAs
        if fuel in ["wind", "solar"]:
            # get a list of diba located in the same region and located in the same time zone
            ba_dibas = list(
                dibas.loc[
                    (dibas.ba_code == ba)
                    & (dibas.ba_region == dibas.diba_region)
                    & (dibas.timezone_local == dibas.timezone_local_diba),
                    "diba_code",
                ].unique()
            )
            if len(ba_dibas) > 0:
                df_temporary = impute_hourly_profiles.average_diba_wind_solar_profiles(
                    data_to_validate, ba, fuel, report_date, ba_dibas, True
                )

                hourly_profiles_to_add.append(df_temporary)
            # if there are no neighboring DIBAs, calculate a national average profile
            else:
                pass

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )

    # calculate the correlations
    compare_method = data_to_validate.merge(
        hourly_profiles_to_add,
        how="left",
        on=[
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "ba_code",
        ],
        validate="1:1",
    )

    compare_method = (
        compare_method.groupby(["fuel_category", "report_date", "ba_code"])
        .corr()
        .reset_index()
    )
    compare_method = compare_method[compare_method["level_3"] == "eia930_profile"]

    compare_method = (
        compare_method.groupby(["fuel_category", "ba_code"])
        .mean()["imputed_profile"]
        .reset_index()
    )

    compare_method = compare_method.rename(
        columns={"imputed_profile": "average_correlation_coefficient"}
    )

    return compare_method


def validate_national_imputation_method(hourly_profiles):

    # only keep wind and solar data
    data_to_validate = hourly_profiles[
        (hourly_profiles["fuel_category"].isin(["wind", "solar"]))
        & (~hourly_profiles["eia930_profile"].isna())
    ]
    data_to_validate = data_to_validate[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "eia930_profile",
        ]
    ]

    profiles_to_impute = data_to_validate[
        ["ba_code", "fuel_category", "report_date"]
    ].drop_duplicates()

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    for index, row in profiles_to_impute.iterrows():
        ba = row["ba_code"]
        fuel = row["fuel_category"]
        report_date = row["report_date"]

        # for wind and solar, average the wind and solar generation profiles from
        # nearby interconnected BAs
        if fuel in ["wind", "solar"]:
            # get a list of diba located in the same region and located in the same time zone
            df_temporary = impute_hourly_profiles.average_national_wind_solar_profiles(
                data_to_validate, ba, fuel, report_date
            )

        hourly_profiles_to_add.append(df_temporary)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )

    # calculate the correlations
    compare_method = data_to_validate.merge(
        hourly_profiles_to_add,
        how="left",
        on=["fuel_category", "datetime_utc", "report_date", "ba_code"],
        validate="1:1",
    )

    compare_method = (
        compare_method.groupby(["fuel_category", "report_date", "ba_code"])
        .corr()
        .reset_index()
    )
    compare_method = compare_method[compare_method["level_3"] == "eia930_profile"]

    compare_method = (
        compare_method.groupby(["fuel_category", "ba_code"])
        .mean()["imputed_profile"]
        .reset_index()
    )

    compare_method = compare_method.rename(
        columns={"imputed_profile": "average_correlation_coefficient"}
    )

    return compare_method


# VALIDATION NOTEBOOK FUNCTIONS
########################################################################################


def test_for_missing_fuel(df, generation_column):
    missing_fuel_test = df[
        (df[generation_column] > 0)
        & (
            (df["fuel_consumed_for_electricity_mmbtu"].isnull())
            | (df["fuel_consumed_for_electricity_mmbtu"] == 0)
        )
    ]
    if not missing_fuel_test.empty:
        print(
            f"Warning: There are {len(missing_fuel_test)} records where {generation_column} is positive but no fuel consumption is reported. Check `missing_fuel_test` for complete list"
        )

    return missing_fuel_test


def test_for_missing_co2(df):
    missing_co2_test = df[df["co2_mass_lb"].isna() & ~df["fuel_consumed_mmbtu"].isna()]
    if not missing_co2_test.empty:
        print(
            f"Warning: There are {len(missing_co2_test)} records where co2 data is missing. Check `missing_co2_test` for complete list"
        )
    return missing_co2_test


def test_for_missing_data(df, columns_to_test):
    missing_data_test = df[df[columns_to_test].isnull().all(axis=1)]
    if not missing_data_test.empty:
        print(
            f"Warning: There are {len(missing_data_test)} records for which no data was reported. Check `missing_data_test` for complete list"
        )
    return missing_data_test


def test_for_missing_incorrect_prime_movers(df, year):

    # cehck for incorrect PM by comparing to EIA-860 data
    pudl_out = load_data.initialize_pudl_out(year)
    pms_in_eia860 = pudl_out.gens_eia860()[
        ["plant_id_eia", "generator_id", "prime_mover_code"]
    ]
    incorrect_pm_test = df.copy()[["plant_id_eia", "generator_id", "prime_mover_code"]]
    incorrect_pm_test = incorrect_pm_test.merge(
        pms_in_eia860,
        how="left",
        on=["plant_id_eia", "generator_id"],
        suffixes=("_allocated", "_eia860"),
        validate="m:1",
    )
    incorrect_pm_test = incorrect_pm_test[
        incorrect_pm_test["prime_mover_code_allocated"]
        != incorrect_pm_test["prime_mover_code_eia860"]
    ]
    if not incorrect_pm_test.empty:
        print(
            f"Warning: There are {len(incorrect_pm_test)} records for which the allocated prime mover does not match the reported prime mover. Check `incorrect_pm_test` for complete list"
        )

    # check for missing PM code
    missing_pm_test = df[df["prime_mover_code"].isna()]
    if not missing_pm_test.empty:
        print(
            f"Warning: There are {len(missing_pm_test)} records for which no prime mover was assigned. Check `missing_pm_test` for complete list"
        )

    return incorrect_pm_test, missing_pm_test


def test_for_outlier_heat_rates(df):
    # check heat rates
    print("Heat Rate Test")
    # remove non-fossil fuel types
    thermal_generators = df[
        ~df["energy_source_code"].isin(["SUN", "MWH", "WND", "WAT", "WH", "PUR"])
    ]
    heat_rate_test_all = []
    for fuel_type in sorted(
        list(thermal_generators.energy_source_code.dropna().unique())
    ):
        # identify all generators with a given fuel type
        generators = thermal_generators[
            thermal_generators["energy_source_code"] == fuel_type
        ]
        # identify all unique prime mover codes for generators of that fuel type
        for pm in sorted(list(generators.prime_mover_code.dropna().unique())):
            generators_with_pm = generators[generators["prime_mover_code"] == pm]
            # calculate a heat rate for each generator of the given fuel type
            heat_rate = (
                generators_with_pm["fuel_consumed_for_electricity_mmbtu"]
                / generators_with_pm["net_generation_mwh"]
            )
            # calculate descriptive statistics for all nonnegative heat rates
            heat_rate_stats = heat_rate[
                (heat_rate >= 0) & (heat_rate != np.inf)
            ].describe()
            # set the outlier threhshold to 1.5 x IQR
            outlier_threshold = (
                (heat_rate_stats["75%"] - heat_rate_stats["25%"]) * 1.5
            ) + heat_rate_stats["75%"]
            # identify all generators whose heatrate is outside this threshold
            heat_rate_test = generators_with_pm[
                (
                    (heat_rate > outlier_threshold)
                    & (heat_rate != np.inf)
                    & (generators_with_pm["net_generation_mwh"] >= 1)
                    | (heat_rate.round(2) == 0)
                )
            ]
            if not heat_rate_test.empty:
                print(
                    f"    Warning: {len(heat_rate_test)} of {len(generators_with_pm)} records for {fuel_type} generators with {pm} prime mover have heat rate of zero or > {outlier_threshold.round(2)} mmbtu/MWh"
                )
                print(
                    f'             median = {heat_rate_stats["50%"].round(2)}, max = {heat_rate_stats["max"].round(2)}, min = {heat_rate_stats["min"].round(2)}'
                )
                heat_rate_test_all.append(heat_rate_test)

    heat_rate_test_all = pd.concat(heat_rate_test_all, axis=0)[
        [
            "report_date",
            "plant_id_eia",
            "generator_id",
            "energy_source_code",
            "prime_mover_code",
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ]
    ]
    heat_rate_test_all["heat_rate"] = (
        heat_rate_test_all["fuel_consumed_for_electricity_mmbtu"]
        / heat_rate_test_all["net_generation_mwh"]
    )
    return heat_rate_test_all


def test_for_zero_data(df, columns_to_test):
    zero_data_test = df[
        (~df[columns_to_test].isnull().all(axis=1))
        & (df[columns_to_test].sum(axis=1) == 0)
    ]
    if not zero_data_test.empty:
        print(
            f"Warning: There are {len(zero_data_test)} records where all operating data are zero. Check `zero_data_test` for complete list"
        )
    return zero_data_test


def test_gtn_results(df):
    gtn_test = df[df["net_generation_mwh"] > df["gross_generation_mwh"]]
    if not gtn_test.empty:
        print(
            f"Warning: There are {round(len(gtn_test)/len(df)*100, 1)}% of records where net generation > gross generation. See `gtn_test` for details"
        )
    return gtn_test


# EGRID VALIDATION METRIC FUNCTIONS
########################################################################################


def load_egrid_plant_file(year):
    # load plant level data from egrid
    egrid_plant = pd.read_excel(
        f"{downloads_folder()}egrid/egrid{year}_data.xlsx",
        sheet_name=f"PLNT{str(year)[-2:]}",
        header=1,
        usecols=[
            "BACODE",
            "PSTATABB",
            "PLPRMFL",
            "ORISPL",
            "PNAME",
            "PLGENATN",
            "PLGENATR",
            "PLHTIANT",
            "UNCO2",
            "UNHTIT",
            "UNHTIOZT",
            "UNHTISRC",
            "UNHOZSRC",
            "PLCO2AN",
            "CHPFLAG",
        ],
    )
    # calculate total net generation from reported renewable and nonrenewable generation
    egrid_plant["net_generation_mwh"] = (
        egrid_plant["PLGENATN"] + egrid_plant["PLGENATR"]
    )
    egrid_plant = egrid_plant.drop(columns=["PLGENATN", "PLGENATR"])
    # rename the columns
    egrid_plant = egrid_plant.rename(
        columns={
            "BACODE": "ba_code",
            "PSTATABB": "state",
            "PLPRMFL": "plant_primary_fuel",
            "ORISPL": "plant_id_egrid",
            "PNAME": "plant_name",
            "UNHTIT": "fuel_consumed_mmbtu",
            "PLHTIANT": "fuel_consumed_for_electricity_mmbtu",
            "UNCO2": "co2_mass_lb",  # this is actually in tons, but we are converting in the next step
            "PLCO2AN": "co2_mass_lb_adjusted",  # this is actually in tons, but we are converting in the next step
            "CHPFLAG": "chp_flag",
            "UNHTIOZT": "fuel_consumed_mmbtu_ozone_season",
            "UNHTISRC": "fuel_data_source_annual",
            "UNHOZSRC": "fuel_data_source_ozone",
        }
    )

    # convert co2 mass tons to lb
    egrid_plant["co2_mass_lb"] = egrid_plant["co2_mass_lb"] * 2000
    egrid_plant["co2_mass_lb_adjusted"] = egrid_plant["co2_mass_lb_adjusted"] * 2000

    # if egrid has a missing value for co2 for a clean plant, replace with zero
    clean_fuels = ["SUN", "MWH", "WND", "WAT", "WH", "PUR", "NUC"]
    egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(clean_fuels), "co2_mass_lb_adjusted"
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(clean_fuels), "co2_mass_lb_adjusted"
    ].fillna(
        0
    )
    egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(clean_fuels), "co2_mass_lb"
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(clean_fuels), "co2_mass_lb"
    ].fillna(
        0
    )

    # reorder the columns
    egrid_plant = egrid_plant[
        [
            "ba_code",
            "state",
            "plant_id_egrid",
            "plant_name",
            "plant_primary_fuel",
            "chp_flag",
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
            "co2_mass_lb",
            "co2_mass_lb_adjusted",
            "fuel_consumed_mmbtu_ozone_season",
            "fuel_data_source_annual",
            "fuel_data_source_ozone",
        ]
    ]

    # remove any plants that have no reported data
    # NOTE: it seems that egrid includes a lot of proposed projects that are not yet operating, but just has missing data for them
    plants_with_no_data_in_egrid = list(
        egrid_plant[
            egrid_plant[
                [
                    "net_generation_mwh",
                    "fuel_consumed_mmbtu",
                    "fuel_consumed_for_electricity_mmbtu",
                    "co2_mass_lb",
                    "co2_mass_lb_adjusted",
                ]
            ].sum(axis=1)
            == 0
        ]["plant_id_egrid"]
    )
    egrid_plant = egrid_plant[
        ~egrid_plant["plant_id_egrid"].isin(plants_with_no_data_in_egrid)
    ]

    # We also want to remove any plants that are located in Puerto Rico
    egrid_plant = egrid_plant[(egrid_plant["state"] != "PR")]

    # create a column for eia id
    egrid_plant = add_egrid_plant_id(egrid_plant, from_id="egrid", to_id="eia")

    return egrid_plant


def load_egrid_ba_file(year):
    # load egrid BA totals
    egrid_ba = pd.read_excel(
        f"{downloads_folder()}egrid/egrid{year}_data.xlsx",
        sheet_name=f"BA{str(year)[-2:]}",
        header=1,
        usecols=["BANAME", "BACODE", "BAHTIANT", "BANGENAN", "BACO2AN"],
    )
    # rename the columns
    egrid_ba = egrid_ba.rename(
        columns={
            "BANAME": "ba_name",
            "BACODE": "ba_code",
            "BAHTIANT": "fuel_consumed_for_electricity_mmbtu",
            "BANGENAN": "net_generation_mwh",
            "BACO2AN": "co2_mass_lb_adjusted",
        }
    )
    egrid_ba = egrid_ba.sort_values(by="ba_code", ascending=True)
    egrid_ba["co2_mass_lb_adjusted"] = egrid_ba["co2_mass_lb_adjusted"] * 2000

    return egrid_ba


def add_egrid_plant_id(df, from_id, to_id):
    # For plants that have different EPA and EIA plant IDs, the plant ID in eGRID is usually the EPA ID, but sometimes the EIA ID
    # however, there are sometime 2 EIA IDs for a single eGRID ID, so we need to group the data in the EIA table by the egrid id
    # We need to update all of the egrid plant IDs to the EIA plant IDs
    egrid_crosswalk = pd.read_csv(
        f"{manual_folder()}eGRID2020_crosswalk_of_EIA_ID_to_EPA_ID.csv",
        dtype=get_dtypes(),
    )
    id_map = dict(
        zip(
            list(egrid_crosswalk[f"plant_id_{from_id}"]),
            list(egrid_crosswalk[f"plant_id_{to_id}"]),
        )
    )

    df[f"plant_id_{to_id}"] = df[f"plant_id_{from_id}"]
    df[f"plant_id_{to_id}"].update(df[f"plant_id_{to_id}"].map(id_map))

    return df


def compare_plant_level_results_to_egrid(
    plant_data, egrid_plant, PLANTS_MISSING_FROM_EGRID
):
    # standardize column names and index so that the two dfs can be divided
    calculated_to_compare = (
        plant_data.groupby("plant_id_egrid", dropna=False)
        .sum()
        .drop(columns=["plant_id_eia"])
    )

    # drop the plants that have no data in eGRID
    plants_with_no_data_in_egrid = list(
        egrid_plant[
            egrid_plant[
                [
                    "net_generation_mwh",
                    "fuel_consumed_mmbtu",
                    "fuel_consumed_for_electricity_mmbtu",
                    "co2_mass_lb",
                    "co2_mass_lb_adjusted",
                ]
            ].sum(axis=1)
            == 0
        ]["plant_id_egrid"]
    )
    egrid_plant = egrid_plant[
        ~egrid_plant["plant_id_eia"].isin(plants_with_no_data_in_egrid)
    ]

    egrid_to_compare = egrid_plant.set_index(["plant_id_egrid"]).drop(
        columns=["ba_code", "state", "plant_name", "plant_id_eia"]
    )
    # only keep plants that are in the comparison data
    egrid_to_compare = egrid_to_compare[
        egrid_to_compare.index.isin(list(calculated_to_compare.index.unique()))
    ]

    # divide calculated value by egrid value
    compared = (
        calculated_to_compare.div(egrid_to_compare)
        .merge(
            egrid_plant[["plant_id_egrid", "plant_name", "ba_code", "state"]],
            how="left",
            left_index=True,
            right_on="plant_id_egrid",
            validate="1:1",
        )
        .set_index("plant_id_egrid")
    )
    compared["plant_name"] = compared["plant_name"].fillna("unknown")

    # create a dataframe that merges the two sources of data together
    compared_merged = calculated_to_compare.merge(
        egrid_to_compare,
        how="left",
        on="plant_id_egrid",
        suffixes=("_calc", "_egrid"),
        validate="1:1",
    )

    # for each column, change missing values to zero if both values are zero (only nan b/c divide by zero)
    for col in [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb_adjusted",
        "co2_mass_lb",
    ]:
        # identify plants with zero values for both
        plant_ids = list(
            compared_merged[
                (compared_merged[f"{col}_calc"] == 0)
                & (compared_merged[f"{col}_egrid"] == 0)
            ].index
        )
        compared.loc[compared.index.isin(plant_ids), col] = 1

    # for each column, categorize the data based on how far it is off from egrid
    for col in [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb_adjusted",
        "co2_mass_lb",
    ]:
        # add a new column
        compared[f"{col}_status"] = pd.cut(
            x=compared[col],
            bins=[
                -999999999,
                -0.0001,
                0.5,
                0.9,
                0.99,
                0.9999,
                1,
                1.0001,
                1.01,
                1.1,
                1.5,
                999999999,
            ],
            labels=[
                "negative",
                "<50%",
                "-50% to -10%",
                "-10% to -1%",
                "+/-1%",
                "!exact",
                "!exact",
                "+/-1%",
                "+1% to 10%",
                "+10% to 50%",
                ">50%",
            ],
            ordered=False,
        )
        # replace any missing values with missing
        compared[f"{col}_status"] = compared[f"{col}_status"].astype(str)
        compared[f"{col}_status"] = compared[f"{col}_status"].fillna("missing")
        compared[f"{col}_status"] = compared[f"{col}_status"].replace("nan", "missing")
        compared.loc[
            (compared.index.isin(PLANTS_MISSING_FROM_EGRID)), f"{col}_status"
        ] = "not_in_egrid"

        # identify which plants are missing from egrid vs calculated values
    for col in [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb_adjusted",
        "co2_mass_lb",
    ]:
        # identify plants that are missing in egrid
        plants_missing_egrid = list(
            compared_merged[
                (compared_merged[f"{col}_calc"] > 0)
                & (compared_merged[f"{col}_egrid"].isna())
            ].index
        )
        compared.loc[
            compared.index.isin(plants_missing_egrid), f"{col}_status"
        ] = "missing_in_egrid"
        # identify plants that are missing from our calculations
        plants_missing_calc = list(
            compared_merged[
                (compared_merged[f"{col}_calc"].isna())
                & (compared_merged[f"{col}_egrid"] > 0)
            ].index
        )
        compared.loc[
            compared.index.isin(plants_missing_calc), f"{col}_status"
        ] = "missing_in_calc"
        # identify where our calculations are missing a zero value
        plants_missing_zero_calc = list(
            compared_merged[
                (compared_merged[f"{col}_calc"].isna())
                & (compared_merged[f"{col}_egrid"] == 0)
            ].index
        )
        compared.loc[
            compared.index.isin(plants_missing_zero_calc), f"{col}_status"
        ] = "calc_missing_zero_value_from_egrid"
        # identify where egrid has a missing value instead of a zero
        plants_missing_zero_egrid = list(
            compared_merged[
                (compared_merged[f"{col}_calc"] == 0)
                & (compared_merged[f"{col}_egrid"].isna())
            ].index
        )
        compared.loc[
            compared.index.isin(plants_missing_zero_egrid), f"{col}_status"
        ] = "egrid_missing_zero_value_from_calc"
        # identify where egrid has a zero value where we have a positive value
        plants_incorrect_zero_egrid = list(
            compared_merged[
                (compared_merged[f"{col}_calc"] > 0)
                & (compared_merged[f"{col}_egrid"] == 0)
            ].index
        )
        compared.loc[
            compared.index.isin(plants_incorrect_zero_egrid), f"{col}_status"
        ] = "calc_positive_but_egrid_zero"

    # create a dataframe that counts how many plants are in each category
    comparison_count = []
    for col in [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb_adjusted",
        "co2_mass_lb",
    ]:
        count = (
            compared.groupby(f"{col}_status", dropna=False)
            .count()["plant_name"]
            .rename(col)
        )
        count.index = count.index.rename("status")
        comparison_count.append(count)

    comparison_count = pd.concat(comparison_count, axis=1).fillna(0).astype(int)
    comparison_count = pd.concat(
        [comparison_count, pd.DataFrame(comparison_count.sum().rename("Total")).T],
        axis=0,
    )
    return comparison_count, compared


def identify_plants_missing_from_our_calculations(
    egrid_plant, annual_plant_results, year
):
    # identify any plants that are in egrid but not our totals, and any plants that are in our totals, but not egrid
    PLANTS_MISSING_FROM_CALCULATION = list(
        set(egrid_plant["plant_id_eia"].unique())
        - set(annual_plant_results["plant_id_eia"].unique())
    )

    # Which plants are included in eGRID but are missing from our calculations?
    missing_from_calc = egrid_plant[
        egrid_plant["plant_id_egrid"].isin(PLANTS_MISSING_FROM_CALCULATION)
    ]

    # see if any of these plants are retired
    generators_eia860 = load_data.load_pudl_table("generators_eia860", year=year)
    missing_from_calc.merge(
        generators_eia860[
            [
                "plant_id_eia",
                "operational_status",
                "current_planned_operating_date",
                "retirement_date",
            ]
        ].drop_duplicates(),
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )

    return missing_from_calc, PLANTS_MISSING_FROM_CALCULATION


def identify_plants_missing_from_egrid(egrid_plant, annual_plant_results):
    # Which plants are in our calculations, but are missing from eGRID?
    PLANTS_MISSING_FROM_EGRID = list(
        set(annual_plant_results["plant_id_egrid"].unique())
        - set(egrid_plant["plant_id_egrid"].unique())
    )

    plant_names = load_data.load_pudl_table("plants_entity_eia")[
        ["plant_id_eia", "plant_name_eia", "sector_name_eia"]
    ]
    missing_from_egrid = annual_plant_results[
        annual_plant_results["plant_id_egrid"].isin(PLANTS_MISSING_FROM_EGRID)
    ].merge(plant_names, how="left", on="plant_id_eia", validate="m:1")

    return missing_from_egrid, PLANTS_MISSING_FROM_EGRID


def segment_plants_by_known_issues(
    annual_plant_results,
    egrid_plant,
    eia923_allocated,
    pudl_out,
    PLANTS_MISSING_FROM_EGRID,
):
    annual_plant_results_segmented = annual_plant_results.copy()
    # missing plants
    annual_plant_results_segmented["flag_missing_egrid"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_id_eia"].isin(PLANTS_MISSING_FROM_EGRID),
        "flag_missing_egrid",
    ] = 1

    # geothermal
    annual_plant_results_segmented["flag_geothermal"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_primary_fuel"] == "GEO", "flag_geothermal"
    ] = 1

    # nuclear
    annual_plant_results_segmented["flag_nuclear"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_primary_fuel"] == "NUC", "flag_nuclear"
    ] = 1

    # fuel cells
    gens_eia860 = pudl_out.gens_eia860()
    PLANTS_WITH_FUEL_CELLS = list(
        gens_eia860.loc[
            gens_eia860["prime_mover_code"] == "FC", "plant_id_eia"
        ].unique()
    )
    annual_plant_results_segmented["flag_fuel_cell"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_id_eia"].isin(PLANTS_WITH_FUEL_CELLS),
        "flag_fuel_cell",
    ] = 1

    # partial
    # identify all of the plants with generators only report part year to CEMS
    partial_year_reporters = eia923_allocated[
        ["plant_id_eia", "generator_id", "hourly_data_source"]
    ].drop_duplicates()
    PARTIAL_YEAR_PLANTS = list(
        partial_year_reporters.loc[
            partial_year_reporters.duplicated(
                subset=["plant_id_eia", "generator_id"], keep=False
            ),
            "plant_id_eia",
        ].unique()
    )
    annual_plant_results_segmented["flag_partial_year"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_id_eia"].isin(PARTIAL_YEAR_PLANTS),
        "flag_partial_year",
    ] = 1

    # CHP plants
    PLANTS_WITH_CHP = list(
        egrid_plant.loc[egrid_plant["chp_flag"] == "Yes", "plant_id_eia"].unique()
    )
    annual_plant_results_segmented["flag_chp"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_id_eia"].isin(PLANTS_WITH_CHP), "flag_chp"
    ] = 1

    # identify plants that report data to the bf or gen table
    bf_reporter = list(pudl_out.bf_eia923()["plant_id_eia"].unique())
    gen_reporter = list(pudl_out.gen_original_eia923()["plant_id_eia"].unique())
    annual_plant_results_segmented["flag_bf_gen_reporter"] = 0
    annual_plant_results_segmented.loc[
        (
            annual_plant_results_segmented["plant_id_eia"].isin(bf_reporter)
            | annual_plant_results_segmented["plant_id_eia"].isin(gen_reporter)
        ),
        "flag_bf_gen_reporter",
    ] = 1

    return annual_plant_results_segmented


def compare_egrid_fuel_total(plant_data, egrid_plant_df):
    """Calculates the difference in fuel and emissions in our calculation and egrid"""

    # standardize column names and index so that the two dfs can be divided
    calculated_to_compare = (
        plant_data.groupby("plant_id_egrid", dropna=False)
        .sum()
        .drop(columns=["plant_id_eia"])
    )

    # drop the plants that have no data in eGRID
    plants_with_no_data_in_egrid = list(
        egrid_plant_df[
            egrid_plant_df[
                [
                    "net_generation_mwh",
                    "fuel_consumed_mmbtu",
                    "fuel_consumed_for_electricity_mmbtu",
                    "co2_mass_lb",
                    "co2_mass_lb_adjusted",
                ]
            ].sum(axis=1)
            == 0
        ]["plant_id_egrid"]
    )
    egrid_plant_df = egrid_plant_df[
        ~egrid_plant_df["plant_id_eia"].isin(plants_with_no_data_in_egrid)
    ]

    egrid_to_compare = egrid_plant_df.set_index(["plant_id_egrid"]).drop(
        columns=["ba_code", "state", "plant_name", "plant_id_eia"]
    )
    # only keep plants that are in the comparison data
    egrid_to_compare = egrid_to_compare[
        egrid_to_compare.index.isin(list(calculated_to_compare.index.unique()))
    ]

    compare_fuel = calculated_to_compare[["fuel_consumed_mmbtu", "co2_mass_lb"]].merge(
        egrid_to_compare[["fuel_consumed_mmbtu", "co2_mass_lb"]],
        how="left",
        left_index=True,
        right_index=True,
        suffixes=("_calc", "_egrid"),
        validate="1:1",
    )
    compare_fuel["difference_fuel"] = (
        compare_fuel["fuel_consumed_mmbtu_egrid"]
        - compare_fuel["fuel_consumed_mmbtu_calc"]
    )
    compare_fuel["difference_co2"] = (
        compare_fuel["co2_mass_lb_egrid"] - compare_fuel["co2_mass_lb_calc"]
    )

    return compare_fuel


def identify_potential_missing_fuel_in_egrid(pudl_out, year, egrid_plant, cems):
    # load the EIA generator fuel data
    IDX_PM_ESC = [
        "report_date",
        "plant_id_eia",
        "energy_source_code",
        "prime_mover_code",
    ]
    gf = pudl_out.gf_eia923().loc[
        :,
        IDX_PM_ESC
        + [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ],
    ]

    # add egrid plant ids
    egrid_crosswalk = pd.read_csv(
        f"{manual_folder()}eGRID2020_crosswalk_of_EIA_ID_to_EPA_ID.csv"
    )
    eia_to_egrid_id = dict(
        zip(
            list(egrid_crosswalk["plant_id_eia"]),
            list(egrid_crosswalk["plant_id_egrid"]),
        )
    )
    gf["plant_id_egrid"] = gf["plant_id_eia"]
    gf["plant_id_egrid"].update(gf["plant_id_egrid"].map(eia_to_egrid_id))

    # calculate an annual total for each plant
    gf_total = gf.groupby(["plant_id_egrid"]).sum().reset_index()

    # choose a metric to compare
    metric = "fuel_consumed_mmbtu"

    # merge the annual EIA-923 data into the egrid data
    egrid_eia_comparison = (
        egrid_plant[
            ["plant_id_egrid", "plant_name", "ba_code", "plant_primary_fuel", metric]
        ]
        .merge(
            gf_total[["plant_id_egrid", metric]],
            how="outer",
            on="plant_id_egrid",
            suffixes=("_egrid", "_eia923"),
            indicator="source",
            validate="1:1",
        )
        .round(0)
    )
    egrid_eia_comparison[f"{metric}_egrid"] = egrid_eia_comparison[
        f"{metric}_egrid"
    ].fillna(0)
    # calculate an absolute difference and percent difference between the two values
    egrid_eia_comparison["difference"] = (
        egrid_eia_comparison[f"{metric}_egrid"]
        - egrid_eia_comparison[f"{metric}_eia923"]
    )
    egrid_eia_comparison["percent_difference"] = (
        egrid_eia_comparison[f"{metric}_egrid"]
        - egrid_eia_comparison[f"{metric}_eia923"]
    ) / egrid_eia_comparison[f"{metric}_eia923"]
    egrid_eia_comparison.loc[
        egrid_eia_comparison["difference"] == 0, "percent_difference"
    ] = 0

    # add cems data so that we can compare fuel totals
    cems_total = cems.copy()[["plant_id_eia", metric]]
    cems_total["plant_id_egrid"] = cems_total["plant_id_eia"]
    cems_total["plant_id_egrid"].update(
        cems_total["plant_id_egrid"].map(eia_to_egrid_id)
    )
    cems_total = (
        cems_total.groupby("plant_id_egrid")
        .sum()[metric]
        .reset_index()
        .rename(columns={metric: f"{metric}_cems"})
    )

    # merge cems data into egrid
    egrid_eia_comparison = egrid_eia_comparison.merge(
        cems_total, how="outer", on="plant_id_egrid", validate="1:1"
    )

    return egrid_eia_comparison
