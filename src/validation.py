import pandas as pd
import numpy as np

import load_data
import impute_hourly_profiles
from constants import CLEAN_FUELS


# DATA PIPELINE VALIDATION FUNCTIONS
########################################################################################


def validate_year(year):
    """Returns a warning if the year specified is not known to work with the pipeline."""

    earliest_validated_year = 2019
    latest_validated_year = 2021

    if year < earliest_validated_year:
        year_warning = f"""
        ################################################################################
        WARNING: The data pipeline has only been validated to work for years {earliest_validated_year}-{latest_validated_year}.
        Running the pipeline for {year} may cause it to fail or may lead to poor-quality
        or anomalous results. To check on the progress of validating additional years of
        data, see: https://github.com/singularity-energy/open-grid-emissions/issues/117
        ################################################################################
        """
        print(year_warning)
    elif year > latest_validated_year:
        year_warning = f"""
        ################################################################################
        WARNING: The most recent available year of input data is currently {latest_validated_year}.
        Input data for {year} should be available from the EIA in Fall {year+1} and we will
        work to validate that the pipeline works with {year} data as soon as possible
        after the data is released.
        ################################################################################
        """
        raise UserWarning(year_warning)


def check_allocated_gf_matches_input_gf(pudl_out, gen_fuel_allocated):
    """
    Checks that the allocated generation and fuel from EIA-923 matches the input totals.

    We use np.isclose() to identify any values that are off by more than 1e-9% different
    from the total input generation or fuel.
    """
    gf = pudl_out.gf_eia923()
    plant_total_gf = gf.groupby("plant_id_eia")[
        [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ]
    ].sum()
    plant_total_alloc = gen_fuel_allocated.groupby("plant_id_eia")[
        [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ]
    ].sum()
    # calculate the percentage difference between the values
    plant_total_diff = ((plant_total_alloc - plant_total_gf) / plant_total_gf).dropna(
        how="any", axis=0
    )
    # flag rows where the absolute percentage difference is greater than our threshold
    mismatched_allocation = plant_total_diff[
        (~np.isclose(plant_total_diff["fuel_consumed_mmbtu"], 0))
        | (~np.isclose(plant_total_diff["net_generation_mwh"], 0))
    ]
    if len(mismatched_allocation) > 0:
        print(
            "WARNING: Allocated EIA-923 data doesn't match input data for the following plants:"
        )
        print("Percentage Difference:")
        print(mismatched_allocation)
        print("EIA-923 Input Totals:")
        print(plant_total_gf.loc[mismatched_allocation.index, :])
        print("Allocated Totals:")
        print(plant_total_alloc.loc[mismatched_allocation.index, :])


def test_for_negative_values(df, small: bool = False):
    """Checks that there are no unexpected negative values in the data."""
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
        "consumed_co2_rate_lb_per_mwh_for_electricity",
        "consumed_ch4_rate_lb_per_mwh_for_electricity",
        "consumed_n2o_rate_lb_per_mwh_for_electricity",
        "consumed_co2e_rate_lb_per_mwh_for_electricity",
        "consumed_nox_rate_lb_per_mwh_for_electricity",
        "consumed_so2_rate_lb_per_mwh_for_electricity",
        "consumed_co2_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_ch4_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_n2o_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_co2e_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_nox_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_so2_rate_lb_per_mwh_for_electricity_adjusted",
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
                f"WARNING: There are {len(negative_test)} records where {column} is negative."
            )
            negative_warnings += 1
    if negative_warnings > 0:
        if small:
            print(
                " Found negative values during small run, these may be fixed with full data"
            )
        else:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("WARNING: The above negative values are errors and must be fixed")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            # raise UserWarning("The above negative values are errors and must be fixed")
    else:
        print("OK")
    return negative_test


def test_for_missing_values(df, small: bool = False):
    """Checks that there are no unexpected missing values in the output data."""
    print("    Checking that no values are missing...  ", end="")
    columns_that_should_be_complete = [
        "plant_id_eia",
        "fuel_category",
        "datetime_local",
        "datetime_utc",
        "month",
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
        "consumed_co2_rate_lb_per_mwh_for_electricity",
        "consumed_ch4_rate_lb_per_mwh_for_electricity",
        "consumed_n2o_rate_lb_per_mwh_for_electricity",
        "consumed_co2e_rate_lb_per_mwh_for_electricity",
        "consumed_nox_rate_lb_per_mwh_for_electricity",
        "consumed_so2_rate_lb_per_mwh_for_electricity",
        "consumed_co2_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_ch4_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_n2o_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_co2e_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_nox_rate_lb_per_mwh_for_electricity_adjusted",
        "consumed_so2_rate_lb_per_mwh_for_electricity_adjusted",
    ]
    columns_to_test = [
        col for col in columns_that_should_be_complete if col in df.columns
    ]
    missing_warnings = 0
    for column in columns_to_test:
        missing_test = df[df[column].isna()]
        if not missing_test.empty:
            print(" ")
            print(
                f"WARNING: There are {len(missing_test)} records where {column} is missing."
            )
            missing_warnings += 1
    if missing_warnings > 0:
        if small:
            print(
                " Found missing values during small run, these may be fixed with full data"
            )
        else:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("WARNING: The above missing values are errors and must be fixed")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    else:
        print("OK")
    return missing_test


def test_chp_allocation(df):
    """Checks that the CHP allocation didn't create any anomalous values."""
    print(
        "    Checking that total fuel consumed >= fuel consumed for electricity...  ",
        end="",
    )
    chp_allocation_test = df[
        df["fuel_consumed_for_electricity_mmbtu"] > df["fuel_consumed_mmbtu"]
    ]
    if not chp_allocation_test.empty:
        raise UserWarning(
            f"WARNING: There are {len(chp_allocation_test)} records where fuel consumed for electricity is greater than total fuel consumption. Check `chp_allocation_test` for complete list"
        )
    else:
        print("OK")

    return chp_allocation_test


def test_for_missing_energy_source_code(df):
    """Checks that there are no missing energy source codes associated with non-zero fuel consumption."""
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
            f"WARNING: There are {len(missing_esc_test)} records where there is a missing energy source code associated with non-zero fuel consumption. Check `missing_esc_test` for complete list"
        )
    else:
        print("OK")

    return missing_esc_test


def check_non_missing_cems_co2_values_unchanged(cems_original, cems):
    """Checks that no non-missing CO2 values were modified during the process of filling."""
    print(
        "    Checking that original CO2 data in CEMS was not modified by filling missing values...",
        end="",
    )
    # only keep non-zero and non-missing co2 values, since these should have not been modified
    cems_original = cems_original.loc[
        cems_original["co2_mass_lb"] > 0,
        ["plant_id_eia", "emissions_unit_id_epa", "datetime_utc", "co2_mass_lb"],
    ]
    test_fill = cems_original.merge(
        cems[["plant_id_eia", "emissions_unit_id_epa", "datetime_utc", "co2_mass_lb"]],
        how="left",
        on=["plant_id_eia", "emissions_unit_id_epa", "datetime_utc"],
        validate="1:1",
        suffixes=("_original", "_postfill"),
    )
    test_fill["diff"] = (
        test_fill["co2_mass_lb_postfill"] - test_fill["co2_mass_lb_original"]
    )
    num_nonzero_rows = len(test_fill[test_fill["diff"] != 0])
    if num_nonzero_rows > 0:
        print(" ")
        print(
            f"WARNING: There are {num_nonzero_rows} non-missing CO2 CEMS records that were modified by `fill_cems_missing_co2` in error"
        )
    else:
        print("OK")


def check_removed_data_is_empty(cems):
    """Checks that the rows removed by `data_cleaning.remove_cems_with_zero_monthly_data()` don't actually contain non-zero data"""
    check_that_data_is_zero = cems.loc[
        cems["missing_data_flag"] == "remove",
        [
            "gross_generation_mwh",
            "steam_load_1000_lb",
            "fuel_consumed_mmbtu",
            "co2_mass_lb",
            "nox_mass_lb",
            "so2_mass_lb",
        ],
    ].sum(numeric_only=True)
    if check_that_data_is_zero.sum() > 0:
        print("WARNING: Some data being removed has non-zero data associated with it:")
        print(check_that_data_is_zero)


def test_for_missing_subplant_id(df):
    """Checks if any records are missing a `subplant_id`."""
    print("    Checking that all data has an associated `subplant_id`...  ", end="")
    missing_subplant_test = df[df["subplant_id"].isna()]
    if not missing_subplant_test.empty:
        print(" ")
        print(
            f"WARNING: There are {len(missing_subplant_test)} records for {len(missing_subplant_test[['plant_id_eia']].drop_duplicates())} plants without a subplant ID. See `missing_subplant_test` for details"
        )
    else:
        print("OK")
    return missing_subplant_test


def validate_gross_to_net_conversion(cems, eia923_allocated):
    """checks whether the calculated net generation matches the reported net generation from EIA-923 at the annual plant level."""
    print(
        "    Checking that calculated net generation matches reported net generation in EIA-923...  ",
        end="",
    )
    # merge together monthly subplant totals from EIA and calculated from CEMS
    eia_netgen = (
        eia923_allocated.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )[["net_generation_mwh"]]
        .sum(min_count=1)
        .reset_index()
        .dropna(subset="net_generation_mwh")
    )
    calculated_netgen = (
        cems.groupby(["plant_id_eia", "subplant_id", "report_date"], dropna=False)[
            "net_generation_mwh"
        ]
        .sum()
        .reset_index()
    )
    validated_ng = eia_netgen.merge(
        calculated_netgen,
        how="inner",
        on=["plant_id_eia", "subplant_id", "report_date"],
        suffixes=("_eia", "_calc"),
        validate="1:1",
    )

    validated_ng = validated_ng.groupby("plant_id_eia")[
        ["net_generation_mwh_eia", "net_generation_mwh_calc"]
    ].sum()

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
            f"WARNING: There are {len(cems_net_not_equal_to_eia)} plants where calculated annual net generation does not match EIA annual net generation."
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
                f"WARNING: There are {len(bad_adjustment)} records where {pollutant}_mass_lb_for_electricity > {pollutant}_mass_lb"
            )
            bad_adjustment += 1

        # test that mass_lb >= mass_lb_adjusted
        bad_adjustment = df[
            (df[f"{pollutant}_mass_lb"] < df[f"{pollutant}_mass_lb_adjusted"])
        ]
        if len(bad_adjustment) > 0:
            print(
                f"WARNING: There are {len(bad_adjustment)} records where {pollutant}_mass_lb_adjusted > {pollutant}_mass_lb"
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
                f"WARNING: There are {len(bad_adjustment)} records where {pollutant}_mass_lb_for_electricity_adjusted > {pollutant}_mass_lb_for_electricity"
            )
            bad_adjustment += 1

    # if there were any bad adjustments, raise a userwarning.
    if bad_adjustments > 0:
        raise UserWarning("The above issues with emissions adjustments must be fixed.")
    else:
        print("OK")


def ensure_non_overlapping_data_from_all_sources(
    cems, partial_cems_subplant, partial_cems_plant, eia_data
):
    """Ensures that there is no duplicated subplant-months from each of the four sources of cleaned data."""

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

    partial_cems_subplant_data = partial_cems_subplant[
        ["plant_id_eia", "subplant_id", "report_date"]
    ].drop_duplicates()
    partial_cems_subplant_data["in_partial_cems_subplant"] = 1

    partial_cems_plant_data = partial_cems_plant[
        ["plant_id_eia", "subplant_id", "report_date"]
    ].drop_duplicates()
    partial_cems_plant_data["in_partial_cems_plant"] = 1

    data_overlap = eia_only_data.merge(
        cems_data,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="1:1",
    )
    data_overlap = data_overlap.merge(
        partial_cems_subplant_data,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="1:1",
    )
    data_overlap = data_overlap.merge(
        partial_cems_plant_data,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="1:1",
    )
    data_overlap[
        ["in_eia", "in_cems", "in_partial_cems_subplant", "in_partial_cems_plant"]
    ] = data_overlap[
        ["in_eia", "in_cems", "in_partial_cems_subplant", "in_partial_cems_plant"]
    ].fillna(
        0
    )
    data_overlap["number_of_locations"] = (
        data_overlap["in_eia"]
        + data_overlap["in_cems"]
        + data_overlap["in_partial_cems_subplant"]
        + data_overlap["in_partial_cems_plant"]
    )

    if len(data_overlap[data_overlap["number_of_locations"] > 1]) > 0:
        eia_cems_overlap = data_overlap[
            (data_overlap["in_eia"] == 1) & (data_overlap["in_cems"] == 1)
        ]
        if len(eia_cems_overlap) > 0:
            print(" ")
            print(
                f"WARNING: There are {len(eia_cems_overlap)} subplant-months that exist in both shaped EIA data and CEMS"
            )
        eia_pcs_overlap = data_overlap[
            (data_overlap["in_eia"] == 1)
            & (data_overlap["in_partial_cems_subplant"] == 1)
        ]
        if len(eia_pcs_overlap) > 0:
            print(" ")
            print(
                f"WARNING: There are {len(eia_pcs_overlap)} subplant-months that exist in both shaped EIA data and partial CEMS data"
            )
        cems_pcs_overlap = data_overlap[
            (data_overlap["in_cems"] == 1)
            & (data_overlap["in_partial_cems_subplant"] == 1)
        ]
        if len(cems_pcs_overlap) > 0:
            print(" ")
            print(
                f"WARNING: There are {len(cems_pcs_overlap)} subplant-months that exist in both CEMS data and partial CEMS data"
            )
        eia_pcp_overlap = data_overlap[
            (data_overlap["in_eia"] == 1) & (data_overlap["in_partial_cems_plant"] == 1)
        ]
        if len(eia_pcp_overlap) > 0:
            print(" ")
            print(
                f"WARNING: There are {len(eia_pcp_overlap)} subplant-months that exist in both shaped EIA data and partial CEMS data"
            )
        cems_pcp_overlap = data_overlap[
            (data_overlap["in_cems"] == 1)
            & (data_overlap["in_partial_cems_plant"] == 1)
        ]
        if len(cems_pcp_overlap) > 0:
            print(" ")
            print(
                f"WARNING: There are {len(cems_pcp_overlap)} subplant-months that exist in both CEMS data and partial CEMS data"
            )
        pcs_pcp_overlap = data_overlap[
            (data_overlap["in_partial_cems_subplant"] == 1)
            & (data_overlap["in_partial_cems_plant"] == 1)
        ]
        if len(pcs_pcp_overlap) > 0:
            print(" ")
            print(
                f"WARNING: There are {len(pcs_pcp_overlap)} subplant-months that exist in both CEMS data and partial CEMS data"
            )
        all_overlap = data_overlap[data_overlap["number_of_locations"] == 4]
        if len(all_overlap) > 0:
            print(" ")
            print(
                f"WARNING: There are {len(all_overlap)} subplant-months that exist in shaped EIA data, CEMS data, and partial CEMS data."
            )
        raise UserWarning("The above overlaps must be fixed before proceeding.")
    else:
        print("OK")


def validate_shaped_totals(shaped_eia_data, monthly_eia_data_to_shape, group_keys):
    """Checks that any shaped monthly data still adds up to the monthly total after shaping."""

    print("    Checking that shaped hourly data matches monthly totals...  ", end="")

    monthly_group_keys = group_keys + ["report_date"]

    # aggregate data to ba fuel month
    shaped_data_agg = shaped_eia_data.groupby(monthly_group_keys, dropna=False)[
        ["net_generation_mwh", "fuel_consumed_mmbtu"]
    ].sum()
    eia_data_agg = monthly_eia_data_to_shape.groupby(monthly_group_keys, dropna=False)[
        ["net_generation_mwh", "fuel_consumed_mmbtu"]
    ].sum()

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
            "The data shaping process is changing the monthly total values compared to reported EIA values. This process should only shape the data, not alter it."
        )
    else:
        print("OK")


def validate_unique_datetimes(df, df_name, keys):
    """Validates that there are unique datetimes per group in a dataframe.

    Args:
        df: dataframe containing datetime columns
        df_name: a descriptive name for the dataframe
        keys: list of column names that contain the groups within which datetimes should be unique"""

    for datetime_column in ["datetime_utc", "datetime_local"]:
        if datetime_column in list(df.columns):
            duplicate_dt = df[
                df.duplicated(subset=(keys + [datetime_column]), keep=False)
            ]
            if len(duplicate_dt) > 0:
                print(duplicate_dt)
                raise UserWarning(
                    f"The dataframe {df_name} contains duplicate {datetime_column} values within each group of {keys}. See above output"
                )


# DATA QUALITY METRIC FUNCTIONS
########################################################################################


def hourly_profile_source_metric(
    cems, partial_cems_subplant, partial_cems_plant, shaped_eia_data, plant_attributes
):
    """Calculates the percentage of data whose hourly profile was determined by method"""
    data_metrics = [
        "net_generation_mwh",
        "co2_mass_lb",
        "co2_mass_lb_for_electricity",
        "nox_mass_lb",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb",
        "so2_mass_lb_for_electricity",
    ]

    # add ba codes and fuel categories to all of the data
    cems = cems.merge(
        plant_attributes[["plant_id_eia", "ba_code"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )
    partial_cems_subplant = partial_cems_subplant.merge(
        plant_attributes[["plant_id_eia", "ba_code"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )
    partial_cems_plant = partial_cems_plant.merge(
        plant_attributes[["plant_id_eia", "ba_code"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )

    # determine the source of the hourly profile
    profile_from_cems = (
        cems.groupby(["ba_code"], dropna=False)[data_metrics].sum().reset_index()
    )
    profile_from_cems["profile_method"] = "cems_reported"

    profile_from_partial_cems_subplant = (
        partial_cems_subplant.groupby(["ba_code"], dropna=False)[data_metrics]
        .sum()
        .reset_index()
    )
    profile_from_partial_cems_subplant[
        "profile_method"
    ] = "eia_scaled_partial_cems_subplant"

    profile_from_partial_cems_plant = (
        partial_cems_plant.groupby(["ba_code"], dropna=False)[data_metrics]
        .sum()
        .reset_index()
    )
    profile_from_partial_cems_plant["profile_method"] = "eia_shaped_partial_cems_plant"

    profile_from_eia = (
        shaped_eia_data.groupby(["ba_code", "profile_method"], dropna=False)[
            data_metrics
        ]
        .sum()
        .reset_index()
    )
    profile_from_eia["profile_method"] = (
        "eia_shaped_" + profile_from_eia["profile_method"]
    )

    profile_source = pd.concat(
        [
            profile_from_cems,
            profile_from_partial_cems_subplant,
            profile_from_partial_cems_plant,
            profile_from_eia,
        ]
    )

    # groupby and calculate percentages for the entire country
    national_source = profile_source.groupby("profile_method").sum(numeric_only=True)
    national_source = (national_source / national_source.sum(axis=0)).reset_index()
    national_source["ba_code"] = "US Total"

    profile_source = profile_source.set_index(["ba_code", "profile_method"])
    # calculate percentages by ba
    profile_source = (
        (profile_source / profile_source.groupby(["ba_code"]).sum())
        .round(4)
        .reset_index()
    )

    method_order = {
        "cems_reported": "0_cems_reported",
        "eia_scaled_partial_cems_subplant": "1_eia_scaled_partial_cems_subplant",
        "eia_shaped_partial_cems_plant": "2_eia_shaped_partial_cems_plant",
        "eia_shaped_residual_profile": "3_eia_shaped_residual_profile",
        "eia_shaped_shifted_residual_profile": "4_eia_shaped_shifted_residual_profile",
        "eia_shaped_eia930_profile": "5_eia_shaped_eia930_profile",
        "eia_shaped_cems_profile": "6_eia_shaped_cems_profile",
        "eia_shaped_DIBA_average": "7_eia_shaped_DIBA_average",
        "eia_shaped_national_average": "8_eia_shaped_national_average",
        "eia_shaped_assumed_flat": "9_eia_shaped_assumed_flat",
    }
    profile_source["profile_method"] = profile_source["profile_method"].replace(
        method_order
    )

    profile_source = profile_source.sort_values(by=["ba_code", "profile_method"])

    profile_source["profile_method"] = profile_source["profile_method"].str[2:]

    # concat the national data to the ba data
    profile_source = pd.concat([profile_source, national_source], axis=0)

    return profile_source


def identify_percent_of_data_by_input_source(
    cems,
    partial_cems_subplant,
    partial_cems_plant,
    eia_only_data,
    year,
    plant_attributes,
):
    """Identifies what percent of output data comes from each input source (CEMS or EIA)."""

    columns_to_use = [
        "net_generation_mwh",
        "emitting_net_generation_mwh",
        "co2_mass_lb",
        "co2_mass_lb_for_electricity",
        "co2e_mass_lb",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb",
        "so2_mass_lb_for_electricity",
    ]

    # add data resolution column to data that is based on EIA
    eia_only_data = identify_reporting_frequency(eia_only_data, year)
    partial_cems_subplant = identify_reporting_frequency(partial_cems_subplant, year)
    partial_cems_plant = identify_reporting_frequency(partial_cems_plant, year)

    # add ba codes and plant primary fuel to all of the data
    eia_only_data = eia_only_data.merge(
        plant_attributes[["plant_id_eia", "ba_code", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )
    cems = cems.merge(
        plant_attributes[["plant_id_eia", "ba_code", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )
    partial_cems_subplant = partial_cems_subplant.merge(
        plant_attributes[["plant_id_eia", "ba_code", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )
    partial_cems_plant = partial_cems_plant.merge(
        plant_attributes[["plant_id_eia", "ba_code", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )

    # add a column for fossil-based generation
    # this copies the net generation data if the associated fuel is not clean or geothermal, and otherwise adds a zero
    # use the generator-specific energy source code for the eia data, otherwise use the pliant primary fuel
    eia_only_data = eia_only_data.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.energy_source_code.isin(CLEAN_FUELS + ["GEO"]), x.net_generation_mwh, 0
        )
    )
    cems = cems.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.plant_primary_fuel.isin(CLEAN_FUELS + ["GEO"]), x.net_generation_mwh, 0
        )
    )
    partial_cems_subplant = partial_cems_subplant.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.plant_primary_fuel.isin(CLEAN_FUELS + ["GEO"]), x.net_generation_mwh, 0
        )
    )
    partial_cems_plant = partial_cems_plant.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.plant_primary_fuel.isin(CLEAN_FUELS + ["GEO"]), x.net_generation_mwh, 0
        )
    )

    # associate each dataframe with a data source label
    data_sources = {
        "cems": cems,
        "partial_cems_subplant": partial_cems_subplant,
        "partial_cems_plant": partial_cems_plant,
        "eia": eia_only_data,
    }
    # get a count of the number of observations (subplant-hours) from each source
    source_of_input_data = []
    for name, df in data_sources.items():
        if len(df) == 0:  # Empty df. May occur when running `small`
            print(f"WARNING: data source {name} has zero entries")
            continue
        if name == "eia":
            subplant_data = df.groupby(
                ["ba_code", "plant_id_eia", "subplant_id", "eia_data_resolution"],
                dropna=False,
            )[columns_to_use].sum()
            # because EIA data is not hourly, we have to multiply the number of subplants by the number of hours in a year
            if year % 4 == 0:
                hours_in_year = 8784
            else:
                hours_in_year = 8760
            subplant_data["subplant_hours"] = hours_in_year
            # group the data by resolution
            subplant_data = (
                subplant_data.reset_index()
                .groupby(["ba_code", "eia_data_resolution"], dropna=False)[
                    ["subplant_hours"] + columns_to_use
                ]
                .sum()
                .reset_index()
            )
            subplant_data = subplant_data.rename(
                columns={"eia_data_resolution": "source"}
            )
            subplant_data["source"] = subplant_data["source"].replace(
                {
                    "annual": "eia_annual",
                    "monthly": "eia_monthly",
                    "multiple": "eia_multiple",
                }
            )
            source_of_input_data.append(subplant_data)
        # for the partial cems data
        elif (name == "partial_cems_subplant") | (name == "partial_cems_plant"):
            subplant_data = df.groupby(
                [
                    "ba_code",
                    "plant_id_eia",
                    "subplant_id",
                    "datetime_utc",
                    "eia_data_resolution",
                ],
                dropna=False,
            )[columns_to_use].sum()
            subplant_data["subplant_hours"] = 1
            # group the data by resolution
            subplant_data = (
                subplant_data.reset_index()
                .groupby(["ba_code", "eia_data_resolution"], dropna=False)[
                    ["subplant_hours"] + columns_to_use
                ]
                .sum()
                .reset_index()
            )
            subplant_data = subplant_data.rename(
                columns={"eia_data_resolution": "source"}
            )
            subplant_data["source"] = subplant_data["source"].replace(
                {
                    "annual": "eia_annual",
                    "monthly": "eia_monthly",
                    "multiple": "eia_multiple",
                }
            )
            source_of_input_data.append(subplant_data)
        # for the cems data
        else:
            subplant_data = df.groupby(
                ["ba_code", "plant_id_eia", "subplant_id", "datetime_utc"], dropna=False
            )[columns_to_use].sum()
            subplant_data["subplant_hours"] = 1
            subplant_data["source"] = "cems_hourly"
            # group the data by resolution
            subplant_data = (
                subplant_data.reset_index()
                .groupby(["ba_code", "source"], dropna=False)[
                    ["subplant_hours"] + columns_to_use
                ]
                .sum()
                .reset_index()
            )
            source_of_input_data.append(subplant_data)

    # concat the dataframes together
    source_of_input_data = pd.concat(source_of_input_data, axis=0)

    # groupby and calculate percentages for the entire country
    national_source = source_of_input_data.groupby("source").sum(numeric_only=True)
    national_source = (national_source / national_source.sum(axis=0)).reset_index()
    national_source["ba_code"] = "US Total"

    # calculate percentages by ba
    source_of_input_data = (
        source_of_input_data.groupby(["ba_code", "source"]).sum(numeric_only=True)
        / source_of_input_data.groupby(["ba_code"]).sum(numeric_only=True)
    ).reset_index()
    # concat the national data to the ba data
    source_of_input_data = pd.concat([source_of_input_data, national_source], axis=0)

    return source_of_input_data


def identify_reporting_frequency(eia923_allocated, year):
    """Identifies if EIA data was reported as an annual total or monthly totals.
    Returns input dataframe with `eia_data_resolution` column added"""

    # load data about the respondent frequency for each plant and merge into the EIA-923 data
    pudl_out = load_data.initialize_pudl_out(year)
    plant_frequency = pudl_out.plants_eia860()[
        ["plant_id_eia", "reporting_frequency_code"]
    ].copy()
    plant_frequency["reporting_frequency_code"] = plant_frequency[
        "reporting_frequency_code"
    ].fillna("multiple")
    # rename the column and recode the values
    plant_frequency = plant_frequency.rename(
        columns={"reporting_frequency_code": "eia_data_resolution"}
    )
    plant_frequency["eia_data_resolution"] = plant_frequency[
        "eia_data_resolution"
    ].replace({"A": "annual", "AM": "monthly", "M": "monthly"})
    # merge the data resolution column into the EIA data
    eia_data = eia923_allocated.merge(
        plant_frequency, how="left", on="plant_id_eia", validate="m:1"
    )
    return eia_data


def summarize_annually_reported_eia_data(eia923_allocated, year):
    """Creates table summarizing the percent of final data from annually-reported EIA data."""

    columns_to_summarize = [
        "fuel_consumed_mmbtu",
        "net_generation_mwh",
        "co2_mass_lb",
        "co2_mass_lb_for_electricity",
        "nox_mass_lb",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb",
        "so2_mass_lb_for_electricity",
    ]

    eia_data = identify_reporting_frequency(eia923_allocated, year)

    data_from_annual = (
        eia_data.groupby(["eia_data_resolution"], dropna=False)[
            columns_to_summarize
        ].sum()
        / eia_data[columns_to_summarize].sum()
        * 100
    ).reset_index()

    annual_eia_used = (
        eia_data[eia_data["hourly_data_source"] != "cems"]
        .groupby(["eia_data_resolution"], dropna=False)[columns_to_summarize]
        .sum()
        / eia_data[columns_to_summarize].sum()
        * 100
    ).reset_index()

    multi_source_subplants = (
        eia_data[["plant_id_eia", "subplant_id", "hourly_data_source"]]
        .drop_duplicates()
        .drop(columns="hourly_data_source")
    )
    multi_source_subplants = multi_source_subplants[
        multi_source_subplants.duplicated(subset=["plant_id_eia", "subplant_id"])
    ]
    multi_source_subplants = eia_data.merge(
        multi_source_subplants, how="inner", on=["plant_id_eia", "subplant_id"]
    )
    multi_source_summary = (
        multi_source_subplants.groupby(["eia_data_resolution"], dropna=False)[
            columns_to_summarize
        ].sum()
        / eia_data[columns_to_summarize].sum()
        * 100
    ).reset_index()

    annual_data_summary = pd.concat(
        [
            pd.DataFrame(
                data_from_annual.loc[
                    data_from_annual["eia_data_resolution"] == "annual", :
                ]
                .set_index("eia_data_resolution")
                .rename(
                    index={
                        "annual": "% of EIA-923 input data from EIA annual reporters"
                    }
                )
                .round(2)
            ),
            pd.DataFrame(
                annual_eia_used.loc[
                    annual_eia_used["eia_data_resolution"] == "annual", :
                ]
                .set_index("eia_data_resolution")
                .rename(index={"annual": "% of output data from EIA annual reporters"})
                .round(2)
            ),
            pd.DataFrame(
                multi_source_summary.loc[
                    multi_source_summary["eia_data_resolution"] == "annual", :
                ]
                .set_index("eia_data_resolution")
                .rename(
                    index={
                        "annual": "% of output data mixing CEMS and annually-reported EIA data"
                    }
                )
                .round(2)
            ),
        ],
        axis=0,
    )

    annual_data_summary.rename(columns={"eia_data_resolution": "category"})

    annual_data_summary = annual_data_summary.reset_index()

    return annual_data_summary


def summarize_cems_measurement_quality(cems):
    """Creates a table summarizing what percent of CO2, SO2, and NOx mass in CEMS was measured or imputed from other hourly values"""
    cems_quality = cems[
        [
            "co2_mass_lb",
            "co2_mass_measurement_code",
            "so2_mass_lb",
            "so2_mass_measurement_code",
            "nox_mass_lb",
            "nox_mass_measurement_code",
        ]
    ].copy()

    # convert categorical columns to strings
    cems_quality[
        [
            "co2_mass_measurement_code",
            "so2_mass_measurement_code",
            "nox_mass_measurement_code",
        ]
    ] = cems_quality[
        [
            "co2_mass_measurement_code",
            "so2_mass_measurement_code",
            "nox_mass_measurement_code",
        ]
    ].astype(
        str
    )
    # replace the CEMS mass measurement codes with two categories
    measurement_code_map = {
        "Measured": "Measured",
        "Measured and Substitute": "Measured",
        "LME": "Imputed",
        "Substitute": "Imputed",
        "Imputed": "Imputed",
        "Calculated": "Imputed",
        "Other": "Imputed",
    }
    cems_quality[
        [
            "co2_mass_measurement_code",
            "so2_mass_measurement_code",
            "nox_mass_measurement_code",
        ]
    ] = cems_quality[
        [
            "co2_mass_measurement_code",
            "so2_mass_measurement_code",
            "nox_mass_measurement_code",
        ]
    ].replace(
        measurement_code_map
    )

    cems_quality_summary = []
    # calculate the percent of mass for each pollutant that is measured or imputed
    for pollutant in ["co2", "nox", "so2"]:
        percent = (
            cems_quality.groupby([f"{pollutant}_mass_measurement_code"], dropna=False)[
                f"{pollutant}_mass_lb"
            ].sum()
            / cems_quality[f"{pollutant}_mass_lb"].sum()
        )
        cems_quality_summary.append(percent)
    cems_quality_summary = pd.concat(cems_quality_summary, axis=1).round(4)
    # drop NA values
    cems_quality_summary = cems_quality_summary.loc[["Measured", "Imputed"], :]
    cems_quality_summary = cems_quality_summary.reset_index()

    return cems_quality_summary


def identify_cems_gtn_method(cems):
    method_summary = cems.groupby("gtn_method", dropna=False)[
        "gross_generation_mwh"
    ].sum()
    method_summary = method_summary / method_summary.sum(axis=0)
    method_summary = method_summary.reset_index()
    method_summary["gtn_method"] = method_summary["gtn_method"].astype(str)
    method_summary = method_summary.sort_values(by="gtn_method", axis=0)
    return method_summary


def validate_wind_solar_imputation(hourly_profiles, year):
    """Creates a table showing cross-validaton results of the wind and solar profile imputation method"""

    # calculate the results and merge together
    diba_results = validate_diba_imputation_method(hourly_profiles, year)
    nationaal_results = validate_national_imputation_method(hourly_profiles)

    imputation_results = diba_results.merge(
        nationaal_results, how="outer", on=["fuel_category", "ba_code"], validate="1:1"
    )

    return imputation_results


def validate_diba_imputation_method(hourly_profiles, year):
    """Validates the method for imputing missing wind and solar profiles.

    Calculates an imputed profile for regions where we have actual wind and solar profiles,
    then calculates how well each imputed profile is correlated with the actual profile.
    Calculates the correlation for each month, then calculates an annual average correlation coefficient.
    """

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
        profiles_to_impute["report_date"].dt.year == year
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

    # merge the imputed data with the actual data
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

    # calculate the correlation coefficient for each fleet-month
    compare_method = (
        compare_method.groupby(["fuel_category", "report_date", "ba_code"])
        .corr(numeric_only=True)
        .reset_index()
    )
    compare_method = compare_method[compare_method["level_3"] == "eia930_profile"]

    # calculate the annual average correlation coefficent for each month
    compare_method = (
        compare_method.groupby(["fuel_category", "ba_code"])["imputed_profile"]
        .mean()
        .reset_index()
    )

    compare_method = compare_method.rename(
        columns={"imputed_profile": "diba_method_correlation_coefficient"}
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

    # merge the imputed data with the actual data
    compare_method = data_to_validate.merge(
        hourly_profiles_to_add,
        how="left",
        on=["fuel_category", "datetime_utc", "report_date", "ba_code"],
        validate="1:1",
    )

    # calculate the correlation coefficient for each fleet-month
    compare_method = (
        compare_method.groupby(["fuel_category", "report_date", "ba_code"])
        .corr(numeric_only=True)
        .reset_index()
    )
    compare_method = compare_method[compare_method["level_3"] == "eia930_profile"]

    # calculate the annual average correlation coefficent for each month
    compare_method = (
        compare_method.groupby(["fuel_category", "ba_code"])["imputed_profile"]
        .mean()
        .reset_index()
    )

    compare_method = compare_method.rename(
        columns={"imputed_profile": "national_method_correlation_coefficient"}
    )

    return compare_method


def check_for_anomalous_co2_factors(
    df, plant_attributes, min_threshold=10, max_threshold=15000
):
    """This function checks that all co2 factors fall within a specified range.

    This is a sanity check to make sure that there are not any obviously
    nonsensical values. On the upper end, we check for co2 factors greater than
    15,000 lb/MWh, which could indicate abnormally high fuel consumption or abnormally
    low net generation. On the lower end, any plant that does not have a co2
    factor equal to zero (carbon-free plants) should have an emission factor that is no
    less than 10.
    """

    pollutant = "co2"
    factor = f"{pollutant}_rate"

    factor_anomaly = df.copy()

    factor_anomaly[factor] = (
        factor_anomaly[f"{pollutant}_mass_lb_for_electricity"]
        / factor_anomaly["net_generation_mwh"]
    )

    # identify if any plants meet these thresholds
    factor_anomaly = factor_anomaly[
        ((factor_anomaly[factor] > 0) & (factor_anomaly[factor] < min_threshold))
        | (factor_anomaly[factor] > max_threshold)
    ]
    # remove any infinity values
    factor_anomaly = factor_anomaly[factor_anomaly[factor] != np.inf]

    if len(factor_anomaly) > 0:
        # merge in plant primary fuel data
        factor_anomaly = factor_anomaly.merge(
            plant_attributes[["plant_id_eia", "plant_primary_fuel"]],
            how="left",
            on="plant_id_eia",
            validate="m:1",
        )
        print("Potentially anomalous co2 factors detected for the following plants:")
        print(
            factor_anomaly[
                [
                    "plant_id_eia",
                    "plant_primary_fuel",
                    "net_generation_mwh",
                    "fuel_consumed_for_electricity_mmbtu",
                    f"{pollutant}_mass_lb_for_electricity",
                    factor,
                ]
            ].sort_values(by=factor)
        )


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
            f"WARNING: There are {len(missing_fuel_test)} records where {generation_column} is positive but no fuel consumption is reported. Check `missing_fuel_test` for complete list"
        )

    return missing_fuel_test


def test_for_missing_co2(df):
    missing_co2_test = df[df["co2_mass_lb"].isna() & ~df["fuel_consumed_mmbtu"].isna()]
    if not missing_co2_test.empty:
        print(
            f"WARNING: There are {len(missing_co2_test)} records where co2 data is missing. Check `missing_co2_test` for complete list"
        )
    return missing_co2_test


def test_for_missing_data(df, columns_to_test):
    missing_data_test = df[df[columns_to_test].isnull().all(axis=1)]
    if not missing_data_test.empty:
        print(
            f"WARNING: There are {len(missing_data_test)} records for which no data was reported. Check `missing_data_test` for complete list"
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
            f"WARNING: There are {len(incorrect_pm_test)} records for which the allocated prime mover does not match the reported prime mover. Check `incorrect_pm_test` for complete list"
        )

    # check for missing PM code
    missing_pm_test = df[df["prime_mover_code"].isna()]
    if not missing_pm_test.empty:
        print(
            f"WARNING: There are {len(missing_pm_test)} records for which no prime mover was assigned. Check `missing_pm_test` for complete list"
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
                    f"    WARNING: {len(heat_rate_test)} of {len(generators_with_pm)} records for {fuel_type} generators with {pm} prime mover have heat rate of zero or > {outlier_threshold.round(2)} mmbtu/MWh"
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
            f"WARNING: There are {len(zero_data_test)} records where all operating data are zero. Check `zero_data_test` for complete list"
        )
    return zero_data_test


def test_gtn_results(df):
    gtn_test = df[df["net_generation_mwh"] > df["gross_generation_mwh"]]
    if not gtn_test.empty:
        print(
            f"WARNING: There are {round(len(gtn_test)/len(df)*100, 1)}% of records where net generation > gross generation. See `gtn_test` for details"
        )
    return gtn_test
