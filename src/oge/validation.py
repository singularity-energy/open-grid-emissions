import pandas as pd
import numpy as np

import oge.load_data as load_data
import oge.impute_hourly_profiles as impute_hourly_profiles
import oge.emissions as emissions
from oge.column_checks import get_dtypes
from oge.filepaths import downloads_folder, reference_table_folder
from oge.logging_util import get_logger

logger = get_logger(__name__)


# DATA PIPELINE VALIDATION FUNCTIONS
########################################################################################


def validate_year(year):
    """Returns a warning if the year specified is not known to work with the pipeline."""

    earliest_validated_year = 2019
    latest_validated_year = 2021

    if year < earliest_validated_year:
        year_warning = f"""
        ################################################################################
        The data pipeline has only been validated to work for years {earliest_validated_year}-{latest_validated_year}.
        Running the pipeline for {year} may cause it to fail or may lead to poor-quality
        or anomalous results. To check on the progress of validating additional years of
        data, see: https://github.com/singularity-energy/open-grid-emissions/issues/117
        ################################################################################
        """
        logger.warning(year_warning)
    elif year > latest_validated_year:
        year_warning = f"""
        ################################################################################
        The most recent available year of input data is currently {latest_validated_year}.
        Input data for {year} should be available from the EIA in Fall {year+1} and we will
        work to validate that the pipeline works with {year} data as soon as possible
        after the data is released.
        ################################################################################
        """
        raise UserWarning(year_warning)


def check_allocated_gf_matches_input_gf(year, gen_fuel_allocated):
    """
    Checks that the allocated generation and fuel from EIA-923 matches the input totals.

    We use np.isclose() to identify any values that are off by more than 1e-9% different
    from the total input generation or fuel.
    """
    gf = load_data.load_pudl_table("generation_fuel_eia923", year)
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
        logger.warning("Allocated EIA-923 doesn't match input data for plants:")
        logger.warning("Percentage Difference:")
        logger.warning("\n" + mismatched_allocation.to_string())
        logger.warning("EIA-923 Input Totals:")
        logger.warning(
            "\n" + plant_total_gf.loc[mismatched_allocation.index, :].to_string()
        )
        logger.warning("Allocated Totals:")
        logger.warning(
            "\n" + plant_total_alloc.loc[mismatched_allocation.index, :].to_string()
        )


def flag_possible_primary_fuel_mismatches(plant_primary_fuel):
    """
    Since we do not know exactly how plants are assigned to fuel categories in EIA-930,
    it is possible that the primary fuel we assign to the plant may differ from the
    primary fuel category used in EIA-930. The most likely source of this disconnect is
    if these plants are assigned a primary fuel based on the type of generation with the
    highest nameplate capacity (since this is relatively static over time), rather than
    based on fuel consumption (which may change year-to-year).

    This test identifies where the primary fuel that the pipeline assigns would lead to
    the plant being categorized under a different fuel category than if a capacity-based
    method were used.

    Plants flagged by this test are not necessarily incorrectly assigned, but rather
    this is intended to bring this issue to our attention as a potential source of
    inconsistency.
    """
    test = plant_primary_fuel.copy()[
        ["plant_id_eia", "plant_primary_fuel_from_capacity_mw", "plant_primary_fuel"]
    ]

    for esc_column in ["plant_primary_fuel_from_capacity_mw", "plant_primary_fuel"]:
        # load the fuel category table
        energy_source_groups = pd.read_csv(
            reference_table_folder("energy_source_groups.csv"), dtype=get_dtypes()
        )[["energy_source_code", "fuel_category_eia930"]].rename(
            columns={
                "energy_source_code": esc_column,
                "fuel_category_eia930": f"{esc_column}_category",
            }
        )

        # assign a fuel category to the monthly eia data
        test = test.merge(
            energy_source_groups,
            how="left",
            on=esc_column,
            validate="m:1",
        )

    mismatched_primary_fuels = test[
        (
            test["plant_primary_fuel_from_capacity_mw_category"]
            != test["plant_primary_fuel_category"]
        )
        & (~test["plant_primary_fuel_from_capacity_mw_category"].isna())
    ]

    if len(mismatched_primary_fuels) > 0:
        logger.warning(
            f"There are {len(mismatched_primary_fuels)} plants where the assigned primary fuel doesn't match the capacity-based primary fuel.\nIt is possible that these plants will categorized as a different fuel in EIA-930"
        )
        logger.warning("\n" + mismatched_primary_fuels.to_string())


def test_for_negative_values(df, small: bool = False):
    """Checks that there are no unexpected negative values in the data."""
    logger.info("Checking that fuel and emissions values are positive...  ")
    columns_that_can_be_negative = ["net_generation_mwh"]
    negative_warnings = 0
    for column in df.columns:
        # if the column is allowed to be negative, skip the test
        if column in columns_that_can_be_negative:
            pass
        else:
            # if the column is not numeric, skip the test
            if pd.api.types.is_numeric_dtype(df[column].dtype):
                negative_test = df[df[column] < 0]
                if not negative_test.empty:
                    logger.warning(
                        f"There are {len(negative_test)} records where {column} is negative."
                    )
                    negative_warnings += 1
            else:
                pass
    if negative_warnings > 0:
        if small:
            logger.warning(
                " Found negative values during small run, these may be fixed with full data"
            )
        else:
            logger.warning("The above negative values are errors and must be fixed!")
    else:
        logger.info("OK")
    return negative_test


def test_for_missing_values(df, small: bool = False):
    """Checks that there are no unexpected missing values in the output data."""
    logger.info("Checking that no values are missing...  ")
    missing_warnings = 0
    for column in df.columns:
        missing_test = df[df[column].isna()]
        if not missing_test.empty:
            logger.warning(
                f"There are {len(missing_test)} records where {column} is missing."
            )
            missing_warnings += 1
    if missing_warnings > 0:
        if small:
            logger.warning(
                "Found missing values during small run, these may be fixed with full data"
            )
        else:
            logger.warning("The above missing values are errors and must be fixed")
    else:
        logger.info("OK")
    return missing_test


def check_for_orphaned_cc_part_in_subplant(subplant_crosswalk):
    """
    Combined cycle generators contain a steam part (CA) and turbine part (CT) that are
    linked together. Thus, our subplant groups that contain one part of a combined cycle
    plant should always in theory contain the other part as well. This test checks that
    both parts exist in a subplant if one exists.

    Besides CT and CA prime movers, there is also CS prime movers which represent a
    "single shaft" combined cycle unit where the steam and turbine parts share a single
    generator. These prime movers are allowed to be by themselves in a subplant, as are
    CC prime movers, which represent a "total unit."
    """
    cc_pm_codes = ["CA", "CT", "CS", "CC"]
    # keep all rows that contain a combined cycle prime mover part
    cc_subplants = subplant_crosswalk[
        subplant_crosswalk["prime_mover_code"].isin(cc_pm_codes)
    ]
    # for each subplant, identify a list of all CC prime movers in that subplant
    cc_subplants = cc_subplants.groupby(["plant_id_eia", "subplant_id"])[
        "prime_mover_code"
    ].agg(["unique"])
    cc_subplants["unique_cc_pms"] = [
        ",".join(map(str, L)) for L in cc_subplants["unique"]
    ]
    cc_subplants = cc_subplants.drop(columns="unique")
    # identify where there are subplants that only contain a single CC part
    orphaned_cc_parts = cc_subplants[
        (cc_subplants["unique_cc_pms"] == "CA")
        | (cc_subplants["unique_cc_pms"] == "CT")
    ]
    if len(orphaned_cc_parts) > 0:
        logger.warning(
            f"There are {len(orphaned_cc_parts)} subplants that only contain one part of a combined cycle system.\nSubplants that represent combined cycle generation should contain both CA and CT parts."
        )
        logger.warning("\n" + orphaned_cc_parts.to_string())


def test_chp_allocation(df):
    """Checks that the CHP allocation didn't create any anomalous values."""
    logger.info(
        "Checking that total fuel consumed >= fuel consumed for electricity...  "
    )
    chp_allocation_test = df[
        df["fuel_consumed_for_electricity_mmbtu"] > df["fuel_consumed_mmbtu"]
    ]
    if not chp_allocation_test.empty:
        raise UserWarning(
            f"There are {len(chp_allocation_test)} records where fuel consumed for electricity is greater than total fuel consumption. Check `chp_allocation_test` for complete list"
        )
    else:
        logger.info("OK")

    return chp_allocation_test


def test_for_missing_energy_source_code(df):
    """Checks that there are no missing energy source codes associated with non-zero fuel consumption."""
    logger.info(
        "Checking that there are no missing energy source codes associated with non-zero fuel consumption...  "
    )
    missing_esc_test = df[
        (df["energy_source_code"].isna()) & (df["fuel_consumed_mmbtu"] > 0)
    ]
    if not missing_esc_test.empty:
        logger.warning(
            f"There are {len(missing_esc_test)} records where there is a missing energy source code associated with non-zero fuel consumption. Check `missing_esc_test` for complete list"
        )
    else:
        logger.info("OK")

    return missing_esc_test


def check_non_missing_cems_co2_values_unchanged(cems_original, cems):
    """Checks that no non-missing CO2 values were modified during the process of filling."""
    logger.info(
        "Checking that original CO2 data in CEMS was not modified by filling missing values...",
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
        logger.warning(
            f"There are {num_nonzero_rows} non-missing CO2 CEMS records that were modified by `fill_cems_missing_co2` in error"
        )
    else:
        logger.info("OK")


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
        logger.warning("Some data being removed has non-zero data associated with it:")
        logger.warning("\n" + check_that_data_is_zero.to_string())


def test_for_missing_subplant_id(df):
    """Checks if any records are missing a `subplant_id`."""
    logger.info("Checking that all data has an associated `subplant_id`...  ")
    missing_subplant_test = df[df["subplant_id"].isna()]
    if not missing_subplant_test.empty:
        logger.warning(
            f"There are {len(missing_subplant_test)} records for {len(missing_subplant_test[['plant_id_eia']].drop_duplicates())} plants without a subplant ID. See `missing_subplant_test` for details"
        )
    else:
        logger.info("OK")
    return missing_subplant_test


def check_missing_or_zero_generation_matches(combined_gen_data):
    """checks that gross generation is positive when net generation is positive.

    This could indicate an issue with missing data or incorrect subplant matching.
    """

    # identify when there is zero or NA gross generation associated with positive net generation
    missing_gross_gen = combined_gen_data[
        (combined_gen_data["net_generation_mwh"] > 0)
        & (combined_gen_data["gross_generation_mwh"] == 0)
    ]

    # identify when there is zero or NA net generation associated with nonzero gross generation
    missing_net_gen = combined_gen_data[
        (combined_gen_data["gross_generation_mwh"] > 0)
        & (combined_gen_data["net_generation_mwh"] == 0)
    ]

    if len(missing_gross_gen) > 0:
        unique_plants = len(missing_gross_gen[["plant_id_eia"]].drop_duplicates())
        unique_subplants = len(
            missing_gross_gen[["plant_id_eia", "subplant_id"]].drop_duplicates()
        )
        logger.warning(
            f"There are {unique_subplants} subplants at {unique_plants} plants for which there is zero gross generation associated with positive net generation."
        )
        logger.warning(
            "\n"
            + missing_gross_gen[
                [
                    "plant_id_eia",
                    "subplant_id",
                    "report_date",
                    "gross_generation_mwh",
                    "net_generation_mwh",
                    "data_source",
                ]
            ]
            .head(10)
            .to_string()
            + "\n...\n"
            + missing_gross_gen[
                [
                    "plant_id_eia",
                    "subplant_id",
                    "report_date",
                    "gross_generation_mwh",
                    "net_generation_mwh",
                    "data_source",
                ]
            ]
            .tail(10)
            .to_string()
        )

    if len(missing_net_gen) > 0:
        unique_plants = len(missing_net_gen[["plant_id_eia"]].drop_duplicates())
        unique_subplants = len(
            missing_net_gen[["plant_id_eia", "subplant_id"]].drop_duplicates()
        )
        logger.warning(
            f"There are {unique_subplants} subplants at {unique_plants} plants for which there is zero net generation associated with positive gross generation."
        )
        logger.warning(
            "\n"
            + missing_net_gen[
                [
                    "plant_id_eia",
                    "subplant_id",
                    "report_date",
                    "gross_generation_mwh",
                    "net_generation_mwh",
                    "data_source",
                ]
            ]
            .head(10)
            .to_string()
            + "\n...\n"
            + missing_net_gen[
                [
                    "plant_id_eia",
                    "subplant_id",
                    "report_date",
                    "gross_generation_mwh",
                    "net_generation_mwh",
                    "data_source",
                ]
            ]
            .tail(10)
            .to_string()
        )


def identify_anomalous_annual_plant_gtn_ratios(annual_plant_ratio):
    """Identifies when net generation for a plant is substantially higher than gross generation."""

    anomalous_gtn = annual_plant_ratio[annual_plant_ratio["annual_plant_ratio"] > 1.25]

    if len(anomalous_gtn) > 0:
        logger.warning(
            "The following plants have annual net generation that is >125% of annual gross generation:"
        )
        logger.warning(
            "\n"
            + anomalous_gtn[
                [
                    "plant_id_eia",
                    "gross_generation_mwh",
                    "net_generation_mwh",
                    "annual_plant_ratio",
                ]
            ]
            .sort_values(by="annual_plant_ratio", ascending=False)
            .to_string()
        )


def validate_gross_to_net_conversion(cems, eia923_allocated):
    """checks whether the calculated net generation matches the reported net generation from EIA-923 at the annual plant level."""
    logger.info(
        "Checking that calculated net generation matches reported net generation in EIA-923...  "
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

    # get the gtn method used for these plants
    gtn_method = cems.loc[
        cems["plant_id_eia"].isin(list(cems_net_not_equal_to_eia.index)),
        ["plant_id_eia", "gtn_method"],
    ].drop_duplicates()
    gtn_method["gtn_method"] = gtn_method["gtn_method"].astype(str)
    gtn_method = (
        gtn_method.groupby("plant_id_eia").agg(["unique"]).droplevel(level=1, axis=1)
    )
    cems_net_not_equal_to_eia = cems_net_not_equal_to_eia.merge(
        gtn_method, how="left", left_index=True, right_index=True
    )

    if len(cems_net_not_equal_to_eia) > 0:
        logger.warning(
            f"There are {len(cems_net_not_equal_to_eia)} plants where calculated annual net generation does not match EIA annual net generation."
        )
        logger.warning("\n" + cems_net_not_equal_to_eia.to_string())
    else:
        logger.info("OK")


def test_emissions_adjustments(df):
    """For each emission, tests that mass_lb >= mass_lb_for_electricity >= mass_lb_for_electricity_adjusted."""

    logger.info(
        "Checking that adjusted emission values are less than total emissions...  "
    )

    pollutants = ["co2", "ch4", "n2o", "co2e", "nox", "so2"]

    bad_adjustments = 0

    for pollutant in pollutants:
        # test that mass_lb >= mass_lb_for_electricity
        bad_adjustment = df[
            (df[f"{pollutant}_mass_lb"] < df[f"{pollutant}_mass_lb_for_electricity"])
        ]
        if len(bad_adjustment) > 0:
            logger.warning(
                f"There are {len(bad_adjustment)} records where {pollutant}_mass_lb_for_electricity > {pollutant}_mass_lb"
            )
            bad_adjustment += 1

        # test that mass_lb >= mass_lb_adjusted
        bad_adjustment = df[
            (df[f"{pollutant}_mass_lb"] < df[f"{pollutant}_mass_lb_adjusted"])
        ]
        if len(bad_adjustment) > 0:
            logger.warning(
                f"There are {len(bad_adjustment)} records where {pollutant}_mass_lb_adjusted > {pollutant}_mass_lb"
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
            logger.warning(
                f"There are {len(bad_adjustment)} records where {pollutant}_mass_lb_for_electricity_adjusted > {pollutant}_mass_lb_for_electricity"
            )
            bad_adjustment += 1

    # if there were any bad adjustments, raise a userwarning.
    if bad_adjustments > 0:
        raise UserWarning("The above issues with emissions adjustments must be fixed.")
    else:
        logger.info("OK")


def ensure_non_overlapping_data_from_all_sources(
    cems, partial_cems_subplant, partial_cems_plant, eia_data
):
    """Ensures that there is no duplicated subplant-months from each of the four sources of cleaned data."""

    logger.info("Checking that all data to be combined is unique...  ")

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
    ].fillna(0)
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
            logger.warning(
                f"There are {len(eia_cems_overlap)} subplant-months that exist in both shaped EIA data and CEMS"
            )
        eia_pcs_overlap = data_overlap[
            (data_overlap["in_eia"] == 1)
            & (data_overlap["in_partial_cems_subplant"] == 1)
        ]
        if len(eia_pcs_overlap) > 0:
            logger.warning(
                f"There are {len(eia_pcs_overlap)} subplant-months that exist in both shaped EIA data and partial CEMS data"
            )
        cems_pcs_overlap = data_overlap[
            (data_overlap["in_cems"] == 1)
            & (data_overlap["in_partial_cems_subplant"] == 1)
        ]
        if len(cems_pcs_overlap) > 0:
            logger.warning(
                f"There are {len(cems_pcs_overlap)} subplant-months that exist in both CEMS data and partial CEMS data"
            )
        eia_pcp_overlap = data_overlap[
            (data_overlap["in_eia"] == 1) & (data_overlap["in_partial_cems_plant"] == 1)
        ]
        if len(eia_pcp_overlap) > 0:
            logger.warning(
                f"There are {len(eia_pcp_overlap)} subplant-months that exist in both shaped EIA data and partial CEMS data"
            )
        cems_pcp_overlap = data_overlap[
            (data_overlap["in_cems"] == 1)
            & (data_overlap["in_partial_cems_plant"] == 1)
        ]
        if len(cems_pcp_overlap) > 0:
            logger.warning(
                f"There are {len(cems_pcp_overlap)} subplant-months that exist in both CEMS data and partial CEMS data"
            )
        pcs_pcp_overlap = data_overlap[
            (data_overlap["in_partial_cems_subplant"] == 1)
            & (data_overlap["in_partial_cems_plant"] == 1)
        ]
        if len(pcs_pcp_overlap) > 0:
            logger.warning(
                f"There are {len(pcs_pcp_overlap)} subplant-months that exist in both CEMS data and partial CEMS data"
            )
        all_overlap = data_overlap[data_overlap["number_of_locations"] == 4]
        if len(all_overlap) > 0:
            logger.warning(
                f"There are {len(all_overlap)} subplant-months that exist in shaped EIA data, CEMS data, and partial CEMS data."
            )
        raise UserWarning("The above overlaps must be fixed before proceeding.")
    else:
        logger.info("OK")


def validate_shaped_totals(shaped_eia_data, monthly_eia_data_to_shape, group_keys):
    """Checks that any shaped monthly data still adds up to the monthly total after shaping."""

    logger.info("Checking that shaped hourly data matches monthly totals...  ")

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
        logger.warning(
            "\n"
            + compare[
                (compare["net_generation_mwh"] != 0)
                | (compare["fuel_consumed_mmbtu"] != 0)
            ].to_string()
        )
        raise UserWarning(
            "The data shaping process is changing the monthly total values compared to reported EIA values. This process should only shape the data, not alter it."
        )
    else:
        logger.info("OK")


def validate_unique_datetimes(df, df_name, keys):
    """Validates that there are unique datetimes per group in a dataframe.

    Args:
        df: dataframe containing datetime columns
        df_name: a descriptive name for the dataframe
        keys: list of column names that contain the groups within which datetimes should be unique
    """

    for datetime_column in ["datetime_utc", "datetime_local"]:
        if datetime_column in list(df.columns):
            duplicate_dt = df[
                df.duplicated(subset=(keys + [datetime_column]), keep=False)
            ]
            if len(duplicate_dt) > 0:
                logger.warning("\n" + duplicate_dt.to_string())
                raise UserWarning(
                    f"The dataframe {df_name} contains duplicate {datetime_column} values within each group of {keys}. See above output"
                )


def check_for_complete_timeseries(df, df_name, keys, period):
    """Validates that a timeseries contains complete hourly data.

    If the `period` is a 'year', checks that the length of the timeseries is 8760 (for a
    non-leap year) or 8784 (for a leap year). If the `period` is a 'month', checks that
    the length of the timeseries is equal to the length of the complete date_range
    between the earliest and latest timestamp in a month.

    Args:
        df: dataframe containing datetime columns
        df_name: a descriptive name for the dataframe
        year
        keys: list of column names that contain the groups within which datetimes should be unique
        period: either 'month' or 'year'. Period within which to ensure complete hourly data
    """

    if period == "year":
        # identify the year of the data
        year = df.datetime_utc.dt.year.mode()[0]
        # count the number of timestamps in each group
        test = df.groupby(keys)[["datetime_utc"]].count()
        # if the year is divisible by 4, it is a leap year
        if year % 4 == 0:
            hours_in_year = 8784
        else:
            hours_in_year = 8760
        test["expected_num_hours"] = hours_in_year
        # identify any rows where the number of timestamps is not equal to the total number of hours in the year
        test = test[test["datetime_utc"] != test["expected_num_hours"]]
        if len(test) > 0:
            logger.warning(
                f"There are incomplete timeseries for the following {keys} groups in {df_name}"
            )
            logger.warning("\n" + test.to_string())
    elif period == "month":
        # count the number of timestamps in each group-month
        test = (
            df.groupby(keys + ["report_date"])[["datetime_utc"]]
            .agg(["count", "min", "max"])
            .droplevel(level=0, axis=1)
        )
        # identify the number of hours in a complete date range for that month
        test["expected_num_hours"] = test.apply(
            lambda row: len(pd.date_range(row["min"], row["max"], freq="H")), axis=1
        )
        # identify any rows where the number of timestamps is not equal to the total number of hours in the month
        test = test[test["count"] != test["expected_num_hours"]]
        if len(test) > 0:
            logger.warning(
                f"There are incomplete timeseries for the following {keys} groups in {df_name}"
            )
            logger.warning("\n" + test.to_string())
    else:
        raise UserWarning(
            f"{period} is not a valid value for the `period` argument in `check_for_complete_timeseries`. Value must be 'year' or 'month'"
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
            ~x.energy_source_code.isin(emissions.CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
        )
    )
    cems = cems.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.plant_primary_fuel.isin(emissions.CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
        )
    )
    partial_cems_subplant = partial_cems_subplant.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.plant_primary_fuel.isin(emissions.CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
        )
    )
    partial_cems_plant = partial_cems_plant.assign(
        emitting_net_generation_mwh=lambda x: np.where(
            ~x.plant_primary_fuel.isin(emissions.CLEAN_FUELS + ["GEO"]),
            x.net_generation_mwh,
            0,
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
            logger.warning(f"data source {name} has zero entries")
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
    plant_frequency = load_data.load_pudl_table(
        "plants_eia860", year, columns=["plant_id_eia", "reporting_frequency_code"]
    )
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
    ].astype(str)
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
    ].replace(measurement_code_map)

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
        logger.warning(
            "Potentially anomalous co2 factors detected for the following plants:"
        )
        logger.warning(
            "\n"
            + factor_anomaly[
                [
                    "plant_id_eia",
                    "plant_primary_fuel",
                    "net_generation_mwh",
                    "fuel_consumed_for_electricity_mmbtu",
                    f"{pollutant}_mass_lb_for_electricity",
                    factor,
                ]
            ]
            .sort_values(by=factor)
            .to_string()
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
        logger.warning(
            f"There are {len(missing_fuel_test)} records where {generation_column} is positive but no fuel consumption is reported. Check `missing_fuel_test` for complete list"
        )

    return missing_fuel_test


def test_for_missing_co2(df):
    missing_co2_test = df[df["co2_mass_lb"].isna() & ~df["fuel_consumed_mmbtu"].isna()]
    if not missing_co2_test.empty:
        logger.warning(
            f"There are {len(missing_co2_test)} records where co2 data is missing. Check `missing_co2_test` for complete list"
        )
    return missing_co2_test


def test_for_missing_data(df, columns_to_test):
    missing_data_test = df[df[columns_to_test].isnull().all(axis=1)]
    if not missing_data_test.empty:
        logger.warning(
            f"There are {len(missing_data_test)} records for which no data was reported. Check `missing_data_test` for complete list"
        )
    return missing_data_test


def test_for_missing_incorrect_prime_movers(df, year):
    # cehck for incorrect PM by comparing to EIA-860 data
    pms_in_eia860 = load_data.load_pudl_table(
        "generators_eia860",
        year,
        columns=["plant_id_eia", "generator_id", "prime_mover_code"],
    )
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
        logger.warning(
            f"There are {len(incorrect_pm_test)} records for which the allocated prime mover does not match the reported prime mover. Check `incorrect_pm_test` for complete list"
        )

    # check for missing PM code
    missing_pm_test = df[df["prime_mover_code"].isna()]
    if not missing_pm_test.empty:
        logger.warning(
            f"There are {len(missing_pm_test)} records for which no prime mover was assigned. Check `missing_pm_test` for complete list"
        )

    return incorrect_pm_test, missing_pm_test


def test_for_outlier_heat_rates(df):
    # check heat rates
    logger.warning("Heat Rate Test")
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
                logger.warning(
                    f"{len(heat_rate_test)} of {len(generators_with_pm)} records for {fuel_type} generators with {pm} prime mover have heat rate of zero or > {outlier_threshold.round(2)} mmbtu/MWh"
                )
                logger.warning(
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
        logger.warning(
            f"There are {len(zero_data_test)} records where all operating data are zero. Check `zero_data_test` for complete list"
        )
    return zero_data_test


def test_gtn_results(df):
    gtn_test = df[df["net_generation_mwh"] > df["gross_generation_mwh"]]
    if not gtn_test.empty:
        logger.warning(
            f"There are {round(len(gtn_test)/len(df)*100, 1)}% of records where net generation > gross generation. See `gtn_test` for details"
        )
    return gtn_test


# EGRID VALIDATION METRIC FUNCTIONS
########################################################################################


def load_egrid_plant_file(year):
    # load plant level data from egrid
    egrid_plant = pd.read_excel(
        downloads_folder(f"egrid/egrid{year}_data.xlsx"),
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
            "UNNOX",
            "UNSO2",
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
            "PNAME": "plant_name_eia",
            "UNHTIT": "fuel_consumed_mmbtu",
            "PLHTIANT": "fuel_consumed_for_electricity_mmbtu",
            "UNCO2": "co2_mass_lb",  # this is actually in tons, but we are converting in the next step
            "UNNOX": "nox_mass_lb",  # this is actually in tons, but we are converting in the next step
            "UNSO2": "so2_mass_lb",  # this is actually in tons, but we are converting in the next step
            "PLCO2AN": "co2_mass_lb_for_electricity_adjusted",  # this is actually in tons, but we are converting in the next step
            "CHPFLAG": "chp_flag",
            "UNHTIOZT": "fuel_consumed_mmbtu_ozone_season",
            "UNHTISRC": "fuel_data_source_annual",
            "UNHOZSRC": "fuel_data_source_ozone",
        }
    )

    # convert co2 mass tons to lb
    egrid_plant["co2_mass_lb"] = egrid_plant["co2_mass_lb"] * 2000
    egrid_plant["nox_mass_lb"] = egrid_plant["nox_mass_lb"] * 2000
    egrid_plant["so2_mass_lb"] = egrid_plant["so2_mass_lb"] * 2000
    egrid_plant["co2_mass_lb_for_electricity_adjusted"] = (
        egrid_plant["co2_mass_lb_for_electricity_adjusted"] * 2000
    )

    # if egrid has a missing value for co2 for a clean plant, replace with zero
    egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(emissions.CLEAN_FUELS),
        "co2_mass_lb_for_electricity_adjusted",
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(emissions.CLEAN_FUELS),
        "co2_mass_lb_for_electricity_adjusted",
    ].fillna(0)
    egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(emissions.CLEAN_FUELS), "co2_mass_lb"
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(emissions.CLEAN_FUELS), "co2_mass_lb"
    ].fillna(0)

    # reorder the columns
    egrid_plant = egrid_plant[
        [
            "ba_code",
            "state",
            "plant_id_egrid",
            "plant_name_eia",
            "plant_primary_fuel",
            "chp_flag",
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
            "co2_mass_lb",
            "co2_mass_lb_for_electricity_adjusted",
            "nox_mass_lb",
            "so2_mass_lb",
            "fuel_consumed_mmbtu_ozone_season",
            "fuel_data_source_annual",
            "fuel_data_source_ozone",
        ]
    ]

    # We also want to remove any plants that are located in Puerto Rico
    egrid_plant = egrid_plant[(egrid_plant["state"] != "PR")]

    # create a column for eia id
    egrid_plant = add_egrid_plant_id(egrid_plant, from_id="egrid", to_id="eia")

    return egrid_plant


def load_egrid_ba_file(year):
    # load egrid BA totals
    egrid_ba = pd.read_excel(
        downloads_folder(f"egrid/egrid{year}_data.xlsx"),
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
        reference_table_folder("eGRID2020_crosswalk_of_EIA_ID_to_EPA_ID.csv"),
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
    columns_to_compare = [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb_for_electricity_adjusted",
        "co2_mass_lb",
        "so2_mass_lb",
        "nox_mass_lb",
    ]
    # standardize column names and index so that the two dfs can be divided
    calculated_to_compare = (
        plant_data.groupby("plant_id_egrid", dropna=False)
        .sum()
        .drop(columns=["plant_id_eia"])
    )

    # drop the plants that have no data in eGRID
    plants_with_no_data_in_egrid = list(
        egrid_plant[egrid_plant[columns_to_compare].sum(axis=1) == 0]["plant_id_egrid"]
    )
    egrid_plant = egrid_plant[
        ~egrid_plant["plant_id_eia"].isin(plants_with_no_data_in_egrid)
    ]

    egrid_to_compare = egrid_plant.set_index(["plant_id_egrid"]).drop(
        columns=["ba_code", "state", "plant_name_eia", "plant_id_eia"]
    )
    # only keep plants that are in the comparison data
    egrid_to_compare = egrid_to_compare[
        egrid_to_compare.index.isin(list(calculated_to_compare.index.unique()))
    ]

    # divide calculated value by egrid value
    compared = (
        calculated_to_compare.div(egrid_to_compare)
        .merge(
            egrid_plant[["plant_id_egrid", "plant_name_eia", "ba_code", "state"]],
            how="left",
            left_index=True,
            right_on="plant_id_egrid",
            validate="1:1",
        )
        .set_index("plant_id_egrid")
    )
    compared["plant_name_eia"] = compared["plant_name_eia"].fillna("unknown")

    # create a dataframe that merges the two sources of data together
    compared_merged = calculated_to_compare.merge(
        egrid_to_compare,
        how="left",
        on="plant_id_egrid",
        suffixes=("_calc", "_egrid"),
        validate="1:1",
    )

    # for each column, change missing values to zero if both values are zero (only nan b/c divide by zero)
    for col in columns_to_compare:
        # identify plants with zero values for both
        plant_ids = list(
            compared_merged[
                (compared_merged[f"{col}_calc"] == 0)
                & (compared_merged[f"{col}_egrid"] == 0)
            ].index
        )
        compared.loc[compared.index.isin(plant_ids), col] = 1

    # for each column, categorize the data based on how far it is off from egrid
    for col in columns_to_compare:
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
    for col in columns_to_compare:
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
    for col in columns_to_compare:
        count = (
            compared.groupby(f"{col}_status", dropna=False)
            .count()["plant_name_eia"]
            .rename(col)
        )
        count.index = count.index.rename("status")
        comparison_count.append(count)

    comparison_count = pd.concat(comparison_count, axis=1).fillna(0).astype(int)
    comparison_count = pd.concat(
        [comparison_count, pd.DataFrame(comparison_count.sum().rename("Total")).T],
        axis=0,
    )

    compared = compared_merged.merge(
        compared[
            [
                "plant_name_eia",
                "ba_code",
                "state",
                "net_generation_mwh_status",
                "fuel_consumed_mmbtu_status",
                "fuel_consumed_for_electricity_mmbtu_status",
                "co2_mass_lb_for_electricity_adjusted_status",
                "co2_mass_lb_status",
                "so2_mass_lb_status",
                "nox_mass_lb_status",
            ]
        ],
        how="left",
        left_index=True,
        right_index=True,
    )

    compared = compared[
        [
            "plant_name_eia",
            "ba_code",
            "state",
            "net_generation_mwh_status",
            "net_generation_mwh_calc",
            "net_generation_mwh_egrid",
            "fuel_consumed_mmbtu_status",
            "fuel_consumed_mmbtu_calc",
            "fuel_consumed_mmbtu_egrid",
            "fuel_consumed_for_electricity_mmbtu_status",
            "fuel_consumed_for_electricity_mmbtu_calc",
            "fuel_consumed_for_electricity_mmbtu_egrid",
            "co2_mass_lb_status",
            "co2_mass_lb_calc",
            "co2_mass_lb_egrid",
            "nox_mass_lb_status",
            "nox_mass_lb_calc",
            "nox_mass_lb_egrid",
            "so2_mass_lb_status",
            "so2_mass_lb_calc",
            "so2_mass_lb_egrid",
            "co2_mass_lb_for_electricity_adjusted_status",
            "co2_mass_lb_for_electricity_adjusted_calc",
            "co2_mass_lb_for_electricity_adjusted_egrid",
        ]
    ]

    return comparison_count, compared


def identify_plants_missing_from_our_calculations(
    egrid_plant, annual_plant_results, year
):
    # remove any plants that have no reported data in egrid
    # NOTE: it seems that egrid includes a lot of proposed projects that are not yet operating, but just has missing data for them
    plants_with_no_data_in_egrid = list(
        egrid_plant[
            egrid_plant[
                [
                    "net_generation_mwh",
                    "fuel_consumed_mmbtu",
                    "fuel_consumed_for_electricity_mmbtu",
                    "co2_mass_lb",
                    "co2_mass_lb_for_electricity_adjusted",
                ]
            ].sum(axis=1)
            == 0
        ]["plant_id_egrid"]
    )
    egrid_plant_no_missing = egrid_plant.copy()[
        ~egrid_plant["plant_id_egrid"].isin(plants_with_no_data_in_egrid)
    ]
    # identify any plants that are in egrid but not our totals, and any plants that are in our totals, but not egrid
    PLANTS_MISSING_FROM_CALCULATION = list(
        set(egrid_plant_no_missing["plant_id_eia"].unique())
        - set(annual_plant_results["plant_id_eia"].unique())
    )

    # Which plants are included in eGRID but are missing from our calculations?
    missing_from_calc = egrid_plant_no_missing[
        egrid_plant_no_missing["plant_id_egrid"].isin(PLANTS_MISSING_FROM_CALCULATION)
    ]

    # see if any of these plants are retired
    generator_status = load_data.load_pudl_table(
        "generators_eia860",
        year,
        columns=[
            "plant_id_eia",
            "operational_status",
            "current_planned_generator_operating_date",
            "generator_retirement_date",
        ],
    ).drop_duplicates()
    missing_from_calc.merge(
        generator_status,
        how="left",
        on="plant_id_eia",
        validate="m:m",
    )

    return missing_from_calc, PLANTS_MISSING_FROM_CALCULATION


def identify_plants_missing_from_egrid(egrid_plant, annual_plant_results):
    # Which plants are in our calculations, but are missing from eGRID?
    PLANTS_MISSING_FROM_EGRID = list(
        set(annual_plant_results["plant_id_egrid"].unique())
        - set(egrid_plant["plant_id_egrid"].unique())
    )

    plant_names = load_data.load_pudl_table(
        "plants_entity_eia", columns=["plant_id_eia", "plant_name_eia"]
    )
    missing_from_egrid = annual_plant_results[
        annual_plant_results["plant_id_egrid"].isin(PLANTS_MISSING_FROM_EGRID)
    ].merge(plant_names, how="left", on="plant_id_eia", validate="m:1")

    return missing_from_egrid, PLANTS_MISSING_FROM_EGRID


def segment_plants_by_known_issues(
    annual_plant_results,
    egrid_plant,
    eia923_allocated,
    year,
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
    gens_eia860 = load_data.load_pudl_table("generators_eia860", year)
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
    bf_reporter = list(
        load_data.load_pudl_table(
            "boiler_fuel_eia923", year, columns=["plant_id_eia"]
        ).unique()
    )
    gen_reporter = list(
        load_data.load_pudl_table(
            "generation_eia923", year, columns=["plant_id_eia"]
        ).unique()
    )
    annual_plant_results_segmented["flag_bf_gen_reporter"] = 0
    annual_plant_results_segmented.loc[
        (
            annual_plant_results_segmented["plant_id_eia"].isin(bf_reporter)
            | annual_plant_results_segmented["plant_id_eia"].isin(gen_reporter)
        ),
        "flag_bf_gen_reporter",
    ] = 1

    # identify plants with proposed generators
    status = load_data.load_pudl_table(
        "generators_eia860",
        year,
        columns=["plant_id_eia", "generator_id", "operational_status"],
    )
    plants_with_proposed_gens = list(
        status.loc[
            status["operational_status"] == "proposed",
            "plant_id_eia",
        ].unique()
    )
    plants_with_proposed_generators = status.loc[
        status["plant_id_eia"].isin(plants_with_proposed_gens),
        ["plant_id_eia", "operational_status"],
    ].drop_duplicates()
    entirely_new_plants = list(
        plants_with_proposed_generators.loc[
            (
                ~plants_with_proposed_generators.duplicated(
                    subset="plant_id_eia", keep=False
                )
            )
            & (plants_with_proposed_generators["operational_status"] == "proposed"),
            "plant_id_eia",
        ].unique()
    )
    annual_plant_results_segmented["flag_plant_w_proposed_gen"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_id_eia"].isin(plants_with_proposed_gens),
        "flag_plant_w_proposed_gen",
    ] = 1
    annual_plant_results_segmented["flag_proposed_plant"] = 0
    annual_plant_results_segmented.loc[
        annual_plant_results_segmented["plant_id_eia"].isin(entirely_new_plants),
        "flag_proposed_plant",
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
        columns=["ba_code", "state", "plant_name_eia", "plant_id_eia"]
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


def identify_potential_missing_fuel_in_egrid(year, egrid_plant, cems):
    # load the EIA generator fuel data
    IDX_PM_ESC = [
        "report_date",
        "plant_id_eia",
        "energy_source_code",
        "prime_mover_code",
    ]
    gf = load_data.load_pudl_table(
        "generation_fuel_eia923",
        year,
        columns=IDX_PM_ESC
        + [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ],
    )

    # add egrid plant ids
    egrid_crosswalk = pd.read_csv(
        reference_table_folder("eGRID2020_crosswalk_of_EIA_ID_to_EPA_ID.csv")
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
            [
                "plant_id_egrid",
                "plant_name_eia",
                "ba_code",
                "plant_primary_fuel",
                metric,
            ]
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
        cems_total.groupby("plant_id_egrid")[metric]
        .sum()
        .reset_index()
        .rename(columns={metric: f"{metric}_cems"})
    )

    # merge cems data into egrid
    egrid_eia_comparison = egrid_eia_comparison.merge(
        cems_total, how="outer", on="plant_id_egrid", validate="1:1"
    )

    return egrid_eia_comparison
