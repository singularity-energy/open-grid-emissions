import pandas as pd

import load_data
from emissions import CLEAN_FUELS
from column_checks import get_dtypes
from filepaths import downloads_folder, manual_folder


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
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS),
        "co2_mass_lb_for_electricity_adjusted",
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS),
        "co2_mass_lb_for_electricity_adjusted",
    ].fillna(
        0
    )
    egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS), "co2_mass_lb"
    ] = egrid_plant.loc[
        egrid_plant["plant_primary_fuel"].isin(CLEAN_FUELS), "co2_mass_lb"
    ].fillna(
        0
    )

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
        manual_folder("eGRID2020_crosswalk_of_EIA_ID_to_EPA_ID.csv"),
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
        validate="m:m",
    )

    return missing_from_calc, PLANTS_MISSING_FROM_CALCULATION


def identify_plants_missing_from_egrid(egrid_plant, annual_plant_results):
    # Which plants are in our calculations, but are missing from eGRID?
    PLANTS_MISSING_FROM_EGRID = list(
        set(annual_plant_results["plant_id_egrid"].unique())
        - set(egrid_plant["plant_id_egrid"].unique())
    )

    plant_names = load_data.load_pudl_table(table_name="plants_entity_eia")[
        ["plant_id_eia", "plant_name_eia"]
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

    # identify plants with proposed generators
    status = pudl_out.gens_eia860()[
        ["plant_id_eia", "generator_id", "operational_status"]
    ]
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
        manual_folder("eGRID2020_crosswalk_of_EIA_ID_to_EPA_ID.csv")
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
