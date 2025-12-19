import pandas as pd
import numpy as np
from itertools import product

import pudl.analysis.allocate_gen_fuel as allocate_gen_fuel

import oge.load_data as load_data
import oge.validation as validation
import oge.emissions as emissions
from oge.constants import CLEAN_FUELS, earliest_data_year
from oge.column_checks import get_dtypes, apply_dtypes, DATA_COLUMNS
from oge.filepaths import reference_table_folder, outputs_folder
from oge.helpers import (
    create_plant_ba_table,
    add_subplant_ids_to_df,
    assign_fleet_to_subplant_data,
)
from oge.logging_util import get_logger

logger = get_logger(__name__)


def clean_eia923(
    year: int,
    calculate_nox_emissions: bool = True,
    calculate_so2_emissions: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """This is the coordinating function for cleaning and allocating generation and
    fuel data in EIA-923.

    Args:
        year (int): year to consider.
        calculate_nox_emissions (bool, optional): whether or not clculate NOx emission
            from fuel consumption. Defaults to True.
        calculate_so2_emissions (bool, optional): whether or not clculate So2 emission
            from fuel consumption. Defaults to True.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: the first data frame has the
            generation, fuel alocated and GHG emissions data, and only NOX and SO2
            emissions if requested. The second data frame is the primary fuel table.
            Finally, the third data frame has subplant-specific emission factors, which
            is used for filling missing emissions data in CEMS later in the pipeline.
            All three data frames are at the subplant level.
    """
    # Allocate fuel and generation across each generator-pm-energy source
    gf = load_data.load_pudl_table("out_eia923__monthly_generation_fuel_combined", year)
    bf = load_data.load_pudl_table("out_eia923__monthly_boiler_fuel", year)
    gen = load_data.load_pudl_table("out_eia923__monthly_generation", year)
    gens = load_data.load_pudl_table("out_eia__yearly_generators", year)
    bga = load_data.load_pudl_table("core_eia860__assn_boiler_generator", year)

    # NOTE: As of 12/7/2024, there is a bug in the pudl data where incorrect generators
    # are getting introduced.
    # See: https://github.com/catalyst-cooperative/pudl/issues/3987
    # To fix this, we need to filter `gens` to remove data with a missing
    # "data_maturity" column
    # (As of 11/27/25) this is no longer needed, since these generators should be filtered out
    # in the cleaning proces
    # gens = gens[~gens["data_maturity"].isna()]

    gf, bf, gen, bga, gens = allocate_gen_fuel.select_input_data(
        gf=gf, bf=bf, gen=gen, bga=bga, gens=gens
    )
    gen_fuel_allocated = allocate_gen_fuel.allocate_gen_fuel_by_generator_energy_source(
        gf,
        bf,
        gen,
        bga,
        gens,
        freq="MS",
    )
    # NOTE: instead of running this in pudl, we can load the data directly from pudl.
    # however, we have some changes to this code in the oge_dev
    """
    gen_fuel_allocated = load_data.load_pudl_table(
        "out_eia923__monthly_generation_fuel_by_generator_energy_source", year
    )
    """

    # Remove plant/generator combinations not in EIA-860.
    # This situation can occur because some generators present in the allocated EIA-923
    # data are filtered out of EIA-860 based on their operational status codes (e.g.,
    # 'CN', 'IP', 'P', 'L'). These codes are excluded from the EIA-860 generator list
    # used for subplant assignments. As a result, plant/generator pairs may exist in
    # the allocated EIA-923 data but not in EIA-860, so we remove them here to ensure
    # consistency.
    to_remove = set(
        gen_fuel_allocated[["plant_id_eia", "generator_id"]]
        .drop_duplicates()
        .apply(tuple, axis=1)
    ).difference(
        set(
            load_data.load_complete_eia_generators_for_subplants()[
                ["plant_id_eia", "generator_id"]
            ]
            .drop_duplicates()
            .apply(tuple, axis=1)
        )
    )
    if len(to_remove) > 0:
        # Check that we are not removing any generators with non-zero generation or
        # fuel consumption
        for p, g in to_remove:
            gen_fuel_allocated_to_remove = gen_fuel_allocated[
                (
                    (gen_fuel_allocated["plant_id_eia"] == p)
                    & (gen_fuel_allocated["generator_id"] == g)
                )
            ]
            if (gen_fuel_allocated_to_remove["net_generation_mwh"] > 0).any() or (
                gen_fuel_allocated_to_remove["fuel_consumed_mmbtu"] > 0
            ).any():
                continue
            else:
                logger.warning(
                    f"Removing ({p},{g}) not in EIA-860 with zero generation and fuel "
                    "consumption"
                )
                gen_fuel_allocated = gen_fuel_allocated[
                    ~(
                        (gen_fuel_allocated["plant_id_eia"] == p)
                        & (gen_fuel_allocated["generator_id"] == g)
                    )
                ]

    # drop bad data where there is negative fuel consumption
    # NOTE(greg) this is in response to a specific issue with the input data for
    # plant 10613 in May 2022, where the data is reported incorrectly in the source
    # data from EIA. EIA has been notified to fix this as of 12/15/2023
    for column in ["fuel_consumed_mmbtu", "fuel_consumed_for_electricity_mmbtu"]:
        bad_fuel_data = gen_fuel_allocated[gen_fuel_allocated[column] < 0]
        if len(bad_fuel_data) > 0:
            logger.warning("Bad input fuel data detected for the following generators:")
            logger.warning(
                bad_fuel_data[
                    [
                        "report_date",
                        "plant_id_eia",
                        "generator_id",
                        "energy_source_code",
                        "prime_mover_code",
                        column,
                    ]
                ]
                .merge(
                    create_plant_ba_table(year)[["plant_id_eia", "ba_code"]],
                    how="left",
                    on="plant_id_eia",
                    validate="m:1",
                )
                .to_string()
            )
            logger.warning("These values will be treated as missing values")
            gen_fuel_allocated.loc[gen_fuel_allocated[column] < 0, column] = np.NaN

    # test to make sure allocated totals match input totals
    validation.check_allocated_gf_matches_input_gf(year, gen_fuel_allocated)

    # manually update energy source code when OTH
    gen_fuel_allocated = update_energy_source_codes(gen_fuel_allocated, year)

    # round all values to the nearest tenth of a unit
    gen_fuel_allocated.loc[
        :,
        [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ],
    ] = gen_fuel_allocated.loc[
        :,
        [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ],
    ].round(1)

    validation.test_for_missing_energy_source_code(gen_fuel_allocated)
    validation.test_for_negative_values(gen_fuel_allocated, year)

    # create a table that identifies the primary fuel of each generator and plant
    primary_fuel_table = create_primary_fuel_table(gen_fuel_allocated, year)

    # calculate co2 emissions for each generator-fuel based on allocated fuel
    # consumption
    gen_fuel_allocated = emissions.calculate_ghg_emissions_from_fuel_consumption(
        df=gen_fuel_allocated,
        year=year,
        include_co2=True,
        include_ch4=True,
        include_n2o=True,
    )

    # Calculate NOx and SO2 emissions
    if calculate_nox_emissions:
        gen_fuel_allocated = emissions.calculate_nox_from_fuel_consumption(
            gen_fuel_allocated, year
        )
    if calculate_so2_emissions:
        gen_fuel_allocated = emissions.calculate_so2_from_fuel_consumption(
            gen_fuel_allocated, year
        )

    # adjust total emissions for biomass
    gen_fuel_allocated = emissions.adjust_emissions_for_biomass(gen_fuel_allocated)

    # adjust emissions for CHP
    logger.info("Adding subplant_id to gen_fuel_allocated for CHP adjustment")
    gen_fuel_allocated = add_subplant_ids_to_df(
        gen_fuel_allocated,
        year,
        plant_part_to_map="generator_id",
        how_merge="left",
        validate_merge="m:1",
    )
    gen_fuel_allocated = emissions.adjust_fuel_and_emissions_for_chp(gen_fuel_allocated)

    gen_fuel_allocated = emissions.calculate_co2e_mass(
        gen_fuel_allocated, year, gwp_horizon=100, ar5_climate_carbon_feedback=True
    )

    validation.test_emissions_adjustments(gen_fuel_allocated, year)

    # calculate weighted emission factors for each subplant-month
    subplant_emission_factors = calculate_subplant_efs(gen_fuel_allocated)

    # before aggregating, output a table of allocated fuels for each subplant (for GRETA)
    subplant_923 = (
        gen_fuel_allocated.groupby(
            by=["plant_id_eia", "subplant_id", "report_date", "energy_source_code"]
        )[["net_generation_mwh", "fuel_consumed_mmbtu"]]
        .sum()
        .reset_index()
    )
    subplant_923 = validation.identify_reporting_frequency(subplant_923, year)

    # aggregate the allocated data to the generator level
    gen_fuel_allocated = (
        gen_fuel_allocated.groupby(by=["report_date", "plant_id_eia", "generator_id"])[
            DATA_COLUMNS
        ]
        .sum(min_count=1)
        .reset_index()
        .pipe(apply_dtypes)
    )

    # remove any plants that we don't want in the data
    gen_fuel_allocated = remove_plants(
        gen_fuel_allocated,
        year,
        non_grid_connected=True,
        remove_states=["PR"],
        remove_steam_only=False,
        distribution_connected_plants=False,
    )

    # round all values to the nearest tenth of a unit
    gen_fuel_allocated.loc[:, DATA_COLUMNS] = gen_fuel_allocated.loc[
        :, DATA_COLUMNS
    ].round(1)

    # map subplant ids to the data
    logger.info("Adding subplant_id to gen_fuel_allocated")
    gen_fuel_allocated = add_subplant_ids_to_df(
        gen_fuel_allocated,
        year,
        plant_part_to_map="generator_id",
        how_merge="left",
        validate_merge="m:1",
    )

    # add the cleaned prime mover code to the data
    gen_pm = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year,
        columns=["plant_id_eia", "generator_id", "prime_mover_code"],
    )
    gen_fuel_allocated = gen_fuel_allocated.merge(
        gen_pm, how="left", on=["plant_id_eia", "generator_id"], validate="m:1"
    )

    gen_fuel_allocated = apply_dtypes(gen_fuel_allocated)
    primary_fuel_table = apply_dtypes(primary_fuel_table)

    # run validation checks on EIA-923 data
    validation.test_for_negative_values(gen_fuel_allocated, year)

    return (
        gen_fuel_allocated,
        primary_fuel_table,
        subplant_emission_factors,
        subplant_923,
    )


def update_energy_source_codes(df, year):
    """Manually update fuel source codes."""
    # load the table of updated fuel types
    updated_esc = pd.read_csv(
        reference_table_folder("updated_oth_energy_source_codes.csv")
    )

    for index, row in updated_esc.iterrows():
        plant_id = row["plant_id_eia"]
        updated_code = row["updated_energy_source_code"]
        df.loc[
            (df["plant_id_eia"] == plant_id) & (df["energy_source_code"] == "OTH"),
            "energy_source_code",
        ] = updated_code

    # print warning if any plants are still other and have nonzero fuel consumption
    plants_with_other_fuel = df[
        (df["energy_source_code"] == "OTH") & (df["fuel_consumed_mmbtu"] > 0)
    ]
    if len(plants_with_other_fuel) > 0:
        id2ba = (
            create_plant_ba_table(year).set_index("plant_id_eia")["ba_code"].to_dict()
        )
        plants_to_check = {
            i: id2ba[i] for i in plants_with_other_fuel["plant_id_eia"].unique()
        }
        logger.warning(
            f"""
            After cleaning energy source codes, some fuel consumption is still 
            associated with an 'OTH' fuel type. This will lead to incorrect emissions 
            calculations. Check the following plants: {plants_to_check}. Assign a fuel 
            type in `data_cleaning.update_energy_source_codes`"""
        )

    return df


def create_primary_fuel_table(
    gen_fuel_allocated: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Identifies the primary fuel for each generator and plant Gen primary fuel is
    identified based on the "energy source code 1" identified in EIA-860. Plant primary
    fuel is based on the most-consumed fuel at a plant based on allocated heat input.

    Args:
        gen_fuel_allocated (pd.DataFrame): data frame of generation and fuel allocated.
        year (int): year under consideration.

    Returns:
        pd.DataFrame: the primary fuel table.
    """
    logger.info("Creating Primary Fuel Table")
    # add under construction generators to the dataframe
    gen_fuel_allocated = add_under_construction_generator_ids_to_df(
        gen_fuel_allocated, year
    )
    gen_fuel_allocated = add_recently_retired_generator_ids_to_df(
        gen_fuel_allocated, year
    )

    # add subplant ids so that we can create subplant-specific primary fuels
    logger.info("Adding subplant_id to gen_fuel_allocated for primary_fuel_table")
    gen_fuel_allocated = add_subplant_ids_to_df(
        gen_fuel_allocated,
        year,
        plant_part_to_map="generator_id",
        how_merge="left",
        validate_merge="m:1",
    )

    # flag and remove any missing ESCs
    gen_fuel_allocated = validation.test_for_missing_energy_source_code(
        gen_fuel_allocated, drop_missing=True
    )

    # get a table of primary energy source codes by generator
    # this will be used in `calculate_aggregated_primary_fuel()` to determine the
    # mode of energy source codes by plant
    # sum the fuel consumption by ESC within each generator
    gen_primary_fuel = (
        gen_fuel_allocated.groupby(
            ["plant_id_eia", "subplant_id", "generator_id", "energy_source_code"],
            dropna=False,
        )["fuel_consumed_mmbtu"]
        .sum()
        .reset_index()
    )

    # only keep the ESC associated with the highest fuel consumption for each gen
    gen_primary_fuel = gen_primary_fuel.sort_values(
        by=["plant_id_eia", "subplant_id", "generator_id", "fuel_consumed_mmbtu"],
        ascending=True,
    ).drop_duplicates(
        subset=["plant_id_eia", "subplant_id", "generator_id"], keep="last"
    )[["plant_id_eia", "subplant_id", "generator_id", "energy_source_code"]]

    # calculate the subplant primary fuel
    subplant_primary_fuel = calculate_aggregated_primary_fuel(
        gen_fuel_allocated,
        gen_primary_fuel,
        "subplant",
        year,
    )
    primary_fuel_table = gen_primary_fuel.merge(
        subplant_primary_fuel,
        how="outer",
        on=["plant_id_eia", "subplant_id"],
        validate="many_to_one",
    )

    plant_primary_fuel = calculate_aggregated_primary_fuel(
        gen_fuel_allocated, gen_primary_fuel, "plant", year
    )

    validation.flag_possible_primary_fuel_mismatches(plant_primary_fuel, year)

    # merge the plant primary fuel into the gen primary fuel
    primary_fuel_table = primary_fuel_table.merge(
        plant_primary_fuel,
        how="outer",
        on="plant_id_eia",
        validate="many_to_one",
    )
    # load manual manual assignment of primary fuel
    # TODO: temporary fix. Need to be removed when EIA adds missing generators at
    # plants listed in file
    primary_fuel_manual = pd.read_csv(
        reference_table_folder("temporary_primary_fuel_manual.csv")
    )
    primary_fuel_manual = primary_fuel_manual[primary_fuel_manual["year"] == year]

    # Drop columns from manual table in-place before concatenation
    primary_fuel_manual.drop(columns=["year", "notes"], inplace=True)

    # Only concatenate if primary_fuel_manual is not empty
    if not primary_fuel_manual.empty:
        primary_fuel_table = pd.concat(
            [primary_fuel_table, primary_fuel_manual],
            axis=0,
            ignore_index=True,
            copy=False,
        ).sort_values(by=["plant_id_eia", "subplant_id"])
    else:
        primary_fuel_table = primary_fuel_table.sort_values(
            by=["plant_id_eia", "subplant_id"]
        )

    return primary_fuel_table


def calculate_aggregated_primary_fuel(
    gen_fuel_allocated, gen_primary_fuel, agg_level, year
):
    """
    Takes generator-level fuel data and calculates primary fuel for the subplant or plant level.

    Args:
        gen_fuel_allocated: dataframe of allocated fuel, generation, and emissions data by generator
        gen_primary_fuel: dataframe of primary fuel by generator
        agg_level: either "plant" or "subplant"
    """
    if agg_level == "plant":
        agg_keys = ["plant_id_eia"]
    elif agg_level == "subplant":
        agg_keys = ["plant_id_eia", "subplant_id"]
    else:
        raise UserWarning(
            f"Argument agg_level must be 'plant' or 'subplant', not '{agg_level}'"
        )

    primary_fuel_from_capacity = calculate_capacity_based_primary_fuel(
        gen_fuel_allocated, agg_level, agg_keys, year
    )

    # NOTE: In some rare cases, a plant will have no fuel specified by
    # energy_source_code_1, and will have zero fuel consumption and net generation for
    # all fuel types. When that happens, we simply assign a plant to have the same fuel
    # type as the majority of its generators.
    primary_fuel_from_mode = (
        gen_primary_fuel.groupby(agg_keys, dropna=False)["energy_source_code"]
        .agg(lambda x: pd.Series.mode(x)[0])
        .to_frame()
        .reset_index()
        .rename(columns={"energy_source_code": f"{agg_level}_primary_fuel_from_mode"})
    )

    # create a blank dataframe with all of the plant ids to hold primary fuel data
    agg_primary_fuel = gen_fuel_allocated[agg_keys].drop_duplicates()

    # calculate the total annual fuel consumption, generation, and capacity by fuel type
    #  for each plant
    agg_totals_by_fuel = (
        gen_fuel_allocated.groupby(agg_keys + ["energy_source_code"], dropna=False)[
            ["fuel_consumed_for_electricity_mmbtu", "net_generation_mwh"]
        ]
        .sum()
        .reset_index()
    )

    # we will calculate primary fuel based on the fuel with the most consumption,
    # generation, and capacity
    for source in ["fuel_consumed_for_electricity_mmbtu", "net_generation_mwh"]:
        # only keep values greater than zero so that these can be filled by other
        # methods if non-zero
        primary_fuel_calc = agg_totals_by_fuel[agg_totals_by_fuel[source] > 0]

        # identify the fuel type with the maximum value for each plant
        primary_fuel_calc = primary_fuel_calc[
            primary_fuel_calc.groupby(agg_keys, dropna=False)[source].transform("max")
            == primary_fuel_calc[source]
        ][agg_keys + ["energy_source_code"]]

        # remove duplicate values (if two fuels are both the maximum)
        primary_fuel_calc = primary_fuel_calc.drop_duplicates(
            subset=agg_keys, keep=False
        )
        primary_fuel_calc = primary_fuel_calc.rename(
            columns={"energy_source_code": f"{agg_level}_primary_fuel_from_{source}"}
        )

        # merge the primary fuel into the main table
        agg_primary_fuel = agg_primary_fuel.merge(
            primary_fuel_calc, how="left", on=agg_keys, validate="1:1"
        )

    # merge the primary fuel into the main table
    agg_primary_fuel = agg_primary_fuel.merge(
        primary_fuel_from_capacity, how="outer", on=agg_keys, validate="1:1"
    )

    agg_primary_fuel = agg_primary_fuel.merge(
        primary_fuel_from_mode, how="left", on=agg_keys, validate="1:1"
    )

    # Use the fuel consumption-based primary fuel first, then fill using capacity-based
    # primary fuel, then generation based. Finally, to break all ties, use the energy
    # source code that appears most often for generators of a plant (mode).
    agg_primary_fuel[f"{agg_level}_primary_fuel"] = agg_primary_fuel[
        f"{agg_level}_primary_fuel_from_fuel_consumed_for_electricity_mmbtu"
    ]
    agg_primary_fuel[f"{agg_level}_primary_fuel"] = agg_primary_fuel[
        f"{agg_level}_primary_fuel"
    ].fillna(agg_primary_fuel[f"{agg_level}_primary_fuel_from_capacity_mw"])
    agg_primary_fuel[f"{agg_level}_primary_fuel"] = agg_primary_fuel[
        f"{agg_level}_primary_fuel"
    ].fillna(agg_primary_fuel[f"{agg_level}_primary_fuel_from_net_generation_mwh"])
    agg_primary_fuel[f"{agg_level}_primary_fuel"] = agg_primary_fuel[
        f"{agg_level}_primary_fuel"
    ].fillna(agg_primary_fuel[f"{agg_level}_primary_fuel_from_mode"])

    # sometimes nuclear generators report 0 fuel consumption in EIA-923.
    # To ensure that a nuclear plant does not get assigned a fuel code of
    # a backup generator, we use the nameplate capacity to assign the
    # primary fuel for any plants that contain a nuclear generator
    agg_primary_fuel.loc[
        agg_primary_fuel[f"{agg_level}_primary_fuel_from_capacity_mw"] == "NUC",
        f"{agg_level}_primary_fuel",
    ] = "NUC"

    # check that there are no missing primary fuels
    if len(agg_primary_fuel[agg_primary_fuel[f"{agg_level}_primary_fuel"].isna()]) > 0:
        plants_with_no_primary_fuel = agg_primary_fuel[
            agg_primary_fuel[f"{agg_level}_primary_fuel"].isna()
        ]
        id2ba = (
            create_plant_ba_table(year).set_index("plant_id_eia")["ba_code"].to_dict()
        )
        plants_to_check = {
            i: id2ba[i] for i in plants_with_no_primary_fuel["plant_id_eia"].unique()
        }
        logger.warning(f"Check the following plants: {plants_to_check}")
        raise UserWarning(
            f"{agg_level} primary fuel table contains missing primary fuels.\
            Update method of `create_primary_fuel_table()` to fix"
        )

    return agg_primary_fuel


def calculate_capacity_based_primary_fuel(
    gen_fuel_allocated, agg_level, agg_keys, year
):
    # create a table of primary fuel by nameplate capacity
    gen_capacity = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year=max(earliest_data_year, year - 4),
        end_year=year,
        columns=[
            "report_date",
            "plant_id_eia",
            "generator_id",
            "capacity_mw",
            "energy_source_code_1",
        ],
    )
    # only keep data for the most recent availble data year
    gen_capacity = gen_capacity[
        gen_capacity["report_date"]
        == gen_capacity.groupby(["plant_id_eia", "generator_id"])[
            "report_date"
        ].transform("max")
    ]
    # drop the report date column in-place
    gen_capacity.drop(columns=["report_date"], inplace=True)

    # only keep keys that exist in gen_fuel_allocated
    gen_capacity = gen_capacity.merge(
        gen_fuel_allocated[["plant_id_eia", "generator_id"]].drop_duplicates(),
        how="right",
        on=["plant_id_eia", "generator_id"],
        validate="1:1",
    )

    if "subplant_id" in agg_keys:
        logger.info("Adding subplant_id to gen_capacity")
        gen_capacity = add_subplant_ids_to_df(
            gen_capacity,
            year,
            plant_part_to_map="generator_id",
            how_merge="left",
            validate_merge="m:1",
        )

    gen_capacity = (
        gen_capacity.groupby(agg_keys + ["energy_source_code_1"], dropna=False)[
            "capacity_mw"
        ]
        .sum()
        .reset_index()
    )

    # drop the battery portion of any hybrid plants so that we don't accidentally
    # identify the primary fuel as storage
    gen_capacity = gen_capacity[
        ~(
            (gen_capacity.duplicated(subset="plant_id_eia", keep=False))
            & (gen_capacity.energy_source_code_1 == "MWH")
        )
    ]

    # find the fuel with the greatest capacity
    gen_capacity = gen_capacity[
        gen_capacity.groupby(agg_keys, dropna=False)["capacity_mw"].transform("max")
        == gen_capacity["capacity_mw"]
    ][agg_keys + ["energy_source_code_1"]].rename(
        columns={"energy_source_code_1": f"{agg_level}_primary_fuel_from_capacity_mw"}
    )

    # drop any duplicate entries (if two fuel types have the same nameplate capacity)
    gen_capacity = gen_capacity[~(gen_capacity.duplicated(subset=agg_keys, keep=False))]

    return gen_capacity


def add_under_construction_generator_ids_to_df(
    df: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Adds rows to df for generators that are under construction. Used to ensure
    complete coverage when a generator starts reporting data before coming online

    NOTE: this function may result in the addition of duplicate generator_ids. If that
    is an issue for the context, run drop_duplicates after this function.

    Args:
        df (pd.DataFrame): the df to add generator IDs to
        year (int): the data year

    Returns:
        pd.DataFrame: df with new rows for under construction generator ids added
    """
    # create a table of primary fuel by nameplate capacity
    gen_capacity = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year,
        columns=[
            "plant_id_eia",
            "generator_id",
            "capacity_mw",
            "energy_source_code_1",
            "operational_status",
            "operational_status_code",
        ],
    ).rename(columns={"energy_source_code_1": "energy_source_code"})

    # keep operating generators and proposed generators that are already under construction
    under_construction_status_codes = ["U", "V", "TS", "OT"]
    gen_cap_under_construction = gen_capacity[
        (
            (gen_capacity["operational_status"] == "proposed")
            & (
                gen_capacity["operational_status_code"].isin(
                    under_construction_status_codes
                )
            )
        )
    ]

    # add subplant_ids
    logger.info("Adding subplant_id to gen_cap_under_construction")
    gen_cap_under_construction = add_subplant_ids_to_df(
        gen_cap_under_construction,
        year,
        plant_part_to_map="generator_id",
        how_merge="left",
        validate_merge="m:1",
    )

    columns_to_append = [
        col for col in gen_cap_under_construction.columns if col in df.columns
    ]

    # check that none of the generators to be added are already in the dataframe
    unique_gens = df[["plant_id_eia", "generator_id"]].drop_duplicates()
    gen_cap_under_construction = gen_cap_under_construction.merge(
        unique_gens,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="1:1",
        indicator="copy",
    )
    gen_cap_under_construction = gen_cap_under_construction[
        gen_cap_under_construction["copy"] != "both"
    ]
    gen_cap_under_construction.drop(columns="copy", inplace=True)

    # add under construction plants to this
    df = pd.concat(
        [df, gen_cap_under_construction[columns_to_append]], axis=0, copy=False
    )

    return df


def add_recently_retired_generator_ids_to_df(
    df: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Adds rows to df for generators that are recently retired. Used to ensure
    complete coverage when a generator continues reporting CEMS data even after retiring

    NOTE: this function may result in the addition of duplicate generator_ids. If that
    is an issue for the context, run drop_duplicates after this function.

    Args:
        df (pd.DataFrame): the df to add generator IDs to
        year (int): the data year

    Returns:
        pd.DataFrame: df with new rows for under construction generator ids added
    """
    # create a table of primary fuel by nameplate capacity
    # load data for up to the past 5 years
    gen_capacity = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year=max(earliest_data_year, year - 4),
        end_year=year,
        columns=[
            "report_date",
            "plant_id_eia",
            "generator_id",
            "capacity_mw",
            "energy_source_code_1",
            "operational_status",
            "operational_status_code",
            "generator_retirement_date",
        ],
    ).rename(columns={"energy_source_code_1": "energy_source_code"})

    # keep generators that have retired in teh past 5 years
    # this is based on data reported in the data year
    gen_cap_recently_retired = gen_capacity[
        (
            (gen_capacity["report_date"] == year)
            & (gen_capacity["operational_status"] == "retired")
            & (gen_capacity["operational_status_code"] == "RE")
            & (gen_capacity["generator_retirement_date"].dt.year >= (year - 4))
        )
    ]

    # some generators "silently retire", meaning they just disappear from the data one
    # year rather than being marked as retired. We also want to include these in the
    # recently retired generators data
    # identify which generators have a most recent report date that is prior to the current year
    silent_retirers = gen_capacity[
        gen_capacity.groupby(["plant_id_eia", "generator_id"])["report_date"]
        .transform("max")
        .dt.year
        < year
    ]
    # filter to only include the most recent year of reported data
    silent_retirers = silent_retirers[
        silent_retirers["report_date"]
        == silent_retirers.groupby(["plant_id_eia", "generator_id"])[
            "report_date"
        ].transform("max")
    ]
    # now, only keep units that were existing in that year (not proposed gens that disappear)
    silent_retirers = silent_retirers[
        silent_retirers["operational_status"] == "existing"
    ]

    gen_cap_recently_retired = pd.concat(
        [gen_cap_recently_retired, silent_retirers], axis=0, copy=False
    )
    gen_cap_recently_retired.drop(columns=["report_date"], inplace=True)

    # add subplant_ids
    logger.info("Adding subplant_id to gen_cap_recently_retired")
    gen_cap_recently_retired = add_subplant_ids_to_df(
        gen_cap_recently_retired,
        year,
        plant_part_to_map="generator_id",
        how_merge="left",
        validate_merge="m:1",
    )

    columns_to_append = [
        col for col in gen_cap_recently_retired.columns if col in df.columns
    ]

    # check that none of the generators to be added are already in the dataframe
    unique_gens = df[["plant_id_eia", "generator_id"]].drop_duplicates()
    gen_cap_recently_retired = gen_cap_recently_retired.merge(
        unique_gens,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="1:1",
        indicator="copy",
    )
    gen_cap_recently_retired = gen_cap_recently_retired[
        gen_cap_recently_retired["copy"] != "both"
    ]
    gen_cap_recently_retired.drop(columns="copy", inplace=True)

    # add under construction plants to this
    df = pd.concat(
        [df, gen_cap_recently_retired[columns_to_append]], axis=0, copy=False
    )

    return df


def calculate_subplant_efs(gen_fuel_allocated):
    """
    Calculates weighted emission factors for each subplant-month for filling in missing data.
    """

    # Select only needed columns from gen_fuel_allocated to reduce memory
    subplant_efs = gen_fuel_allocated[
        [
            "plant_id_eia",
            "subplant_id",
            "report_date",
            "fuel_consumed_mmbtu",
            "co2_mass_lb",
            "ch4_mass_lb",
            "n2o_mass_lb",
            "co2e_mass_lb",
            "nox_mass_lb",
            "so2_mass_lb",
        ]
    ].copy()

    # calculate the total emissions and fuel consumption by subplant-month
    subplant_efs = subplant_efs.groupby(
        ["plant_id_eia", "subplant_id", "report_date"], dropna=False
    )[
        [
            "fuel_consumed_mmbtu",
            "co2_mass_lb",
            "ch4_mass_lb",
            "n2o_mass_lb",
            "co2e_mass_lb",
            "nox_mass_lb",
            "so2_mass_lb",
        ]
    ].sum()

    # drop any observations with no fuel consumption
    subplant_efs = subplant_efs[subplant_efs["fuel_consumed_mmbtu"] > 0]

    # calculate the fuel emission factor
    for pollutant in ["co2", "ch4", "n2o", "co2e", "nox", "so2"]:
        subplant_efs[f"{pollutant}_lb_per_mmbtu"] = (
            subplant_efs[f"{pollutant}_mass_lb"] / subplant_efs["fuel_consumed_mmbtu"]
        )

    # only keep relevant columns
    subplant_efs = subplant_efs.reset_index()[
        [
            "plant_id_eia",
            "subplant_id",
            "report_date",
            "co2_lb_per_mmbtu",
            "ch4_lb_per_mmbtu",
            "n2o_lb_per_mmbtu",
            "co2e_lb_per_mmbtu",
            "nox_lb_per_mmbtu",
            "so2_lb_per_mmbtu",
        ]
    ]

    return subplant_efs


def remove_plants(
    df,
    year,
    non_grid_connected=False,
    remove_states=[],
    remove_steam_only=False,
    distribution_connected_plants=False,
):
    """
    Coordinating function to remove specific plants based on specified options
    Each function should identify how many plants are being removed
    Args:
        df: dataframe containing plant_id_eia column
        non_grid_connected: if True, remove all plants that are not grid connected
        remove_states: list of two-letter state codes for which plants should be
        removed if located within remove_steam_only: if True, remove plants that only
        generate heat and no electricity (not yet implemented)
        distribution_connected_plants: if True, remove plants that are connected to the
        distribution grid (not yet implemented)
    """
    if non_grid_connected:
        df = remove_non_grid_connected_plants(df, year)
    if len(remove_states) > 0:
        plant_states = load_data.load_pudl_table(
            "core_eia__entity_plants", columns=["plant_id_eia", "state"]
        )
        plants_in_states_to_remove = list(
            plant_states[
                plant_states["state"].isin(remove_states)
            ].plant_id_eia.unique()
        )

        logger.info(
            f"Removing {len(plants_in_states_to_remove)} plants located in the following states: {remove_states}"
        )
        df = df[~df["plant_id_eia"].isin(plants_in_states_to_remove)]
    if remove_steam_only:
        # df = manually_remove_steam_units(df)
        df = remove_unmapped_fuel(df, year)
        df = identify_and_remove_steam_only_units(df, year)

    if distribution_connected_plants:
        pass

    return df


def remove_non_grid_connected_plants(df: pd.DataFrame, year: int) -> pd.DataFrame:
    """Removes any records from a dataframe associated with plants that are not
    connected to the electricity grid.

    Args:
        df (pd.DataFrame): any pandas dataframe containing the column 'plant_id_eia'
        year (int): The data year

    Returns:
        pd.DataFrame: The df with NGCs removed
    """

    # get the list of plant_id_eia from the static table
    ngc_plants = list(
        pd.read_csv(
            reference_table_folder("plants_not_connected_to_grid.csv"),
            dtype=get_dtypes(),
        )["Plant ID"]
    )

    num_plants = len(
        df[df["plant_id_eia"].isin(ngc_plants)]["plant_id_eia"].unique()
    ) + len(
        df[(df["plant_id_eia"] >= 880000) & (df["plant_id_eia"] < 890000)][
            "plant_id_eia"
        ].unique()
    )
    logger.info(f"Removing {num_plants} plants that are not grid-connected")

    df = df[~df["plant_id_eia"].isin(ngc_plants)]

    # for CEMS data, remove any units marked as "Manual CAMD Excluded".
    # Their documentation notes: "Any CAMD units in the manual match file that should be
    # excluded from the matching process, mostly due to the lack of a connection to the
    # electricity grid (e.g., industrial boilers), are added to the crosswalk with an
    # indicator that they were manually excluded."
    if "emissions_unit_id_epa" in df.columns:
        epa_crosswalk = load_data.load_epa_eia_crosswalk_from_raw(year)[
            ["plant_id_eia", "emissions_unit_id_epa", "epa_match_type"]
        ]
        epa_ngc = epa_crosswalk[
            epa_crosswalk["epa_match_type"] == "Manual CAMD Excluded"
        ].drop_duplicates()
        # merge into dataframe and exclude these plants
        df = df.merge(
            epa_ngc,
            how="left",
            on=["plant_id_eia", "emissions_unit_id_epa"],
            validate="m:1",
        )
        df = df[df["epa_match_type"] != "Manual CAMD Excluded"]
        df.drop(columns=["epa_match_type"], inplace=True)

    # according to the egrid documentation, any plants that have an id of 88XXXX are
    # not grid connected only keep plants that dont have an id of 88XXXX
    df = df[(df["plant_id_eia"] < 880000) | (df["plant_id_eia"] >= 890000)]

    return df


def remove_unmapped_fuel(cems: pd.DataFrame, year: int) -> pd.DataFrame:
    """Removes units that only report fuel input but no outputs, and are unmapped to
    an EIA generator

    Args:
        cems (pd.DataFrame): The CEMS data to assess and remove steam units from
        year (int): The data year

    Returns:
        pd.DataFrame: The CEMS data with any unmapped steam-only units removed
    """
    # calculate annual totals by unit
    annual_cems = (
        cems.groupby(["plant_id_eia", "emissions_unit_id_epa"], dropna=False)[
            ["gross_generation_mwh", "steam_load_1000_lb", "fuel_consumed_mmbtu"]
        ]
        .sum()
        .reset_index()
    )

    # merge in mapped generator ids so that we can check which CEMS data is mapped to
    # an EIA generator, and thus will not be double counted
    subplant_crosswalk = (
        pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv.zip"),
            dtype=get_dtypes(),
        )[["plant_id_eia", "emissions_unit_id_epa", "generator_id"]].drop_duplicates()
    ).dropna(subset=["emissions_unit_id_epa"])
    annual_cems = annual_cems.merge(
        subplant_crosswalk,
        how="left",
        on=["plant_id_eia", "emissions_unit_id_epa"],
        validate="1:m",
    )
    # raise a warning about units that have non-zero generation but are missing a
    # generator_id. This may result in double-counting of generation and emissions
    potential_missing_map = annual_cems[
        (annual_cems["generator_id"].isna()) & (annual_cems["gross_generation_mwh"] > 0)
    ]
    if len(potential_missing_map) > 0:
        logger.error(
            "The following CEMS units report non-zero generation but are missing a mapping to an EIA generator_id"
        )
        logger.error(
            "These may need to be mapped to prevent double-counting of generation"
        )
        logger.error(
            validation.limit_error_output_df(potential_missing_map).to_string(
                index=False
            )
        )

    # flag units that report fuel input but no generation or steam output, and which
    # are not mapped to an EIA unit. This is potentially anomalous data that we want
    # to remove
    fuel_only_unmapped = annual_cems[
        (annual_cems["generator_id"].isna())
        & (annual_cems["gross_generation_mwh"] == 0)
        & (annual_cems["steam_load_1000_lb"] == 0)
        & (annual_cems["fuel_consumed_mmbtu"] > 0)
    ]
    if len(fuel_only_unmapped) > 0:
        logger.warning(
            "The following CEMS units report fuel consumption but no generation or steam output"
        )
        logger.warning(
            "These units also are missing a mapping to an EIA generator, meaning they will have their own subplant ID"
        )
        logger.warning(
            "To prevent these fuel and emissions from being counted, they will be removed"
        )
        logger.warning(
            validation.limit_error_output_df(fuel_only_unmapped).to_string(index=False)
        )

        cems = cems.merge(
            fuel_only_unmapped[
                ["plant_id_eia", "emissions_unit_id_epa"]
            ].drop_duplicates(),
            how="left",
            on=["plant_id_eia", "emissions_unit_id_epa"],
            validate="m:1",
            indicator="fuel_only_unmapped",
        )
        cems = cems[cems["fuel_only_unmapped"] != "both"]
        cems.drop(columns="fuel_only_unmapped", inplace=True)

    return cems


def clean_cems(year: int, primary_fuel_table, subplant_emission_factors):
    """
    Coordinating function for all of the cems data cleaning
    """
    # load the CEMS data
    cems = load_data.load_cems_data(year)
    validation.validate_unique_datetimes(
        year, cems, "cems", ["plant_id_eia", "emissions_unit_id_epa"]
    )

    cems = remove_negative_cems_data(cems, year)

    # remove non-grid connected plants
    cems = remove_plants(
        cems,
        year,
        non_grid_connected=True,
        remove_states=["PR"],
        remove_steam_only=True,
        distribution_connected_plants=False,
    )

    # add a report date
    cems = load_data.add_report_date(cems)

    # remove data for any unit-months where there are incomplete data reported
    # this is generally when there is a single observation reported for an entire month
    cems = remove_incomplete_unit_months(cems)

    # create an inventory of which plant-months have input data from any source
    inventory_input_data_sources(cems, year)

    # Flag outliers in generation, fuel consumption and CO2 emission data
    validation.check_for_outliers_in_cems_generation_fuel_and_co2_time_series(cems)

    # TODO: identify and remove any hourly values that appear to be outliers
    # See: https://github.com/singularity-energy/open-grid-emissions/issues/50

    # add subplant id
    logger.info("Adding subplant_id to cems")
    cems = add_subplant_ids_to_df(
        cems,
        year,
        plant_part_to_map="emissions_unit_id_epa",
        how_merge="left",
        validate_merge="m:1",
    )

    # add a fuel type to each observation
    cems = assign_fuel_type_to_cems(cems, year, primary_fuel_table)

    # fill in missing hourly emissions data using the fuel type and heat input
    validation.test_for_missing_energy_source_code(cems)
    cems = emissions.fill_cems_missing_co2(cems, year, subplant_emission_factors)

    # TODO: Add functions for filling missing NOx and SOx
    # See: https://github.com/singularity-energy/open-grid-emissions/issues/153

    # calculate ch4 and n2o emissions
    cems = emissions.calculate_ghg_emissions_from_fuel_consumption(
        df=cems, year=year, include_co2=False, include_ch4=True, include_n2o=True
    )

    # remove any observations from cems where zero operation is reported for an entire
    # month. Although this data could be considered to be accurately reported, let's
    # remove it so that we can double check against the eia data
    # NOTE(12/22/23): We will treat reported zeros as actual data and not attempt to
    # fill reported zeros with EIA-923 data due to the issues that exist with the method
    # EIA uses to allocate annual data to months. In the future, we could use monthly-
    # reported EIA data since this is directly reported by the generator.
    # NOTE (11/27/25) re-adding this filter to reduce memory issues in the pipeline.
    # cems = remove_cems_with_zero_monthly_data(cems)

    validation.test_for_negative_values(cems, year)
    validation.validate_unique_datetimes(
        year, cems, "cems", ["plant_id_eia", "emissions_unit_id_epa"]
    )

    cems = apply_dtypes(cems)

    return cems


def remove_negative_cems_data(cems: pd.DataFrame, year: int) -> pd.DataFrame:
    """Removes anomalous data, defined as negative fuel or emisisons data.

    Args:
        cems (pd.DataFrame): The CEMS data to screen and remove data from
        year (int): The data year

    Returns:
        pd.DataFrame: CEMS data with anomalous data removed
    """

    # treat negative emissions as bad data and replace with missing values
    for column in ["co2_mass_lb", "nox_mass_lb", "so2_mass_lb"]:
        negative_emissions_data = cems[cems[column] < 0]
        if len(negative_emissions_data) > 0:
            logger.warning(
                f"Bad input {column} data detected for the following plants:"
            )
            logger.warning(
                negative_emissions_data[
                    [
                        "datetime_utc",
                        "plant_id_eia",
                        "plant_id_epa",
                        "gross_generation_mwh",
                        column,
                        f"{column.split('_')[0]}_mass_measurement_code",
                    ]
                ]
                .merge(
                    create_plant_ba_table(year)[["plant_id_eia", "ba_code"]],
                    how="left",
                    on="plant_id_eia",
                    validate="m:1",
                )
                .to_string()
            )
            logger.warning("These values will be treated as missing values")
            cems.loc[cems[column] < 0, column] = np.NaN

    return cems


def manually_remove_steam_units(df):
    """
    Removes any records from CEMS that we've identified as being steam only plants that need to be removed
    """

    # get the list of plant_id_eia from the static table
    units_to_remove = pd.read_csv(
        reference_table_folder("steam_units_to_remove.csv"),
        dtype=get_dtypes(),
    )[["plant_id_eia", "emissions_unit_id_epa"]]

    logger.info(
        f"Removing {len(units_to_remove)} units that only produce steam and do not report to EIA"
    )

    df = df.merge(
        units_to_remove,
        how="outer",
        on=["plant_id_eia", "emissions_unit_id_epa"],
        indicator="source",
        validate="m:1",
    )

    df = df[df["source"] == "left_only"].drop(columns=["source"])

    return df


def identify_and_remove_steam_only_units(cems: pd.DataFrame, year: int) -> pd.DataFrame:
    """Identifies CEMS units that report steam output data and either removes or flags
    the data.

    Where steam-only output data is unmapped to an EIA generator, it is removed.
    Where steam-only output data is mapped to an EIA generator, it is flagged.

    Args:
        cems (pd.DataFrame): The CEMS data to assess and remove steam units from
        year (int): The data year

    Returns:
        pd.DataFrame: The CEMS data with any unmapped steam-only units removed
    """
    # calculate annual totals by unit
    annual_cems = (
        cems.groupby(["plant_id_eia", "emissions_unit_id_epa"], dropna=False)[
            ["gross_generation_mwh", "steam_load_1000_lb", "fuel_consumed_mmbtu"]
        ]
        .sum()
        .reset_index()
    )

    # merge in mapped generator ids so that we can check which CEMS data is mapped to
    # an EIA generator, and thus will not be double counted
    subplant_crosswalk = (
        pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv.zip"),
            dtype=get_dtypes(),
        )[["plant_id_eia", "emissions_unit_id_epa", "generator_id"]].drop_duplicates()
    ).dropna(subset=["emissions_unit_id_epa"])
    annual_cems = annual_cems.merge(
        subplant_crosswalk,
        how="left",
        on=["plant_id_eia", "emissions_unit_id_epa"],
        validate="1:m",
    )

    # flag units that report only steam output but no generation, and are not mapped to
    # an EIA generator. We want to remove these so that the steam only generation is not
    # counted
    steam_only_unmapped = annual_cems[
        (annual_cems["generator_id"].isna())
        & (annual_cems["gross_generation_mwh"] == 0)
        & (annual_cems["steam_load_1000_lb"] > 0)
        & (annual_cems["fuel_consumed_mmbtu"] > 0)
    ]
    if len(steam_only_unmapped) > 0:
        logger.warning(
            "The following CEMS units report only steam output and are missing a generator ID mapping"
        )
        logger.warning(
            "This means they would have their own subplant ID and potentially be double counted, so will be removed"
        )
        logger.warning(
            validation.limit_error_output_df(steam_only_unmapped).to_string(index=False)
        )

        cems = cems.merge(
            steam_only_unmapped[
                ["plant_id_eia", "emissions_unit_id_epa"]
            ].drop_duplicates(),
            how="left",
            on=["plant_id_eia", "emissions_unit_id_epa"],
            validate="m:1",
            indicator="unmapped_steam",
        )
        cems = cems[cems["unmapped_steam"] != "both"]
        cems.drop(columns="unmapped_steam", inplace=True)

    # flag units that report non-zero steam output and are mapped to a generator. We
    # are still not entirely certain how to interpret steam data in CEMS, so its
    # inclusion could result in potentially anomalous results
    mapped_steam = annual_cems[
        (annual_cems["steam_load_1000_lb"] > 0) & (~annual_cems["generator_id"].isna())
    ]
    if len(mapped_steam) > 0:
        logger.warning(
            "The following CEMS units report non-zero steam output and are mapped to an EIA generator"
        )
        logger.warning(
            "This may result in anomalous results for these units until steam output data is handled in the pipeline"
        )
        logger.warning(
            validation.limit_error_output_df(mapped_steam).to_string(index=False)
        )

    return cems


def remove_incomplete_unit_months(cems):
    # get a count of how many hours are reported in each month for each unit
    unit_hours_in_month = (
        cems[["plant_id_eia", "report_date", "emissions_unit_id_epa", "datetime_utc"]]
        .groupby(["plant_id_eia", "report_date", "emissions_unit_id_epa"], dropna=False)
        .count()
        .reset_index()
    )

    # identify months where there is less than a single day of data
    # The fewest number of hours in a month is 28*24 = 672
    unit_months_to_remove = unit_hours_in_month[
        unit_hours_in_month["datetime_utc"] < 24
    ].drop(columns="datetime_utc")

    logger.info(
        f"Removing {len(unit_months_to_remove)} unit-months with incomplete hourly data"
    )

    cems = cems.merge(
        unit_months_to_remove,
        how="outer",
        on=["plant_id_eia", "report_date", "emissions_unit_id_epa"],
        validate="m:1",
        indicator="to_remove",
    )

    cems = cems[cems["to_remove"] != "both"].drop(columns="to_remove")

    return cems


def assign_fuel_type_to_cems(cems, year, primary_fuel_table):
    "Assigns a fuel type to each observation in CEMS"

    # merge in the subplant primary fuel type
    cems = cems.merge(
        primary_fuel_table[
            ["plant_id_eia", "subplant_id", "subplant_primary_fuel"]
        ].drop_duplicates(),
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )
    cems = cems.rename(columns={"subplant_primary_fuel": "energy_source_code"})

    # fill missing fuel codes for plants that only have a single fuel type
    if len(cems[cems["energy_source_code"].isna()]) > 0:
        single_fuel_plants = (
            primary_fuel_table.drop_duplicates(
                subset=["plant_id_eia", "energy_source_code"]
            ).drop_duplicates(subset=["plant_id_eia"], keep=False)
        )[["plant_id_eia", "energy_source_code"]]
        cems = cems.merge(
            single_fuel_plants,
            how="left",
            on=["plant_id_eia"],
            suffixes=(None, "_plant"),
            validate="m:1",
        )
        cems = fillna_with_missing_strings(
            cems,
            column_to_fill="energy_source_code",
            filler_column="energy_source_code_plant",
        )

    # Fill fuel codes for plants that only have a single fossil type identified in EIA
    if len(cems[cems["energy_source_code"].isna()]) > 0:
        cems = fill_missing_fuel_for_single_fuel_plant_months(cems, year)

    # fill any remaining missing fuel codes with the plant primary fuel identified from EIA-923
    if len(cems[cems["energy_source_code"].isna()]) > 0:
        cems = cems.merge(
            primary_fuel_table[
                ["plant_id_eia", "plant_primary_fuel"]
            ].drop_duplicates(),
            how="left",
            on="plant_id_eia",
            validate="m:1",
        )
        cems = fillna_with_missing_strings(
            cems,
            column_to_fill="energy_source_code",
            filler_column="plant_primary_fuel",
        )

    # if there are still missing fuels, the plant might be proposed and not yet in EIA-923
    # in this case, load data from EIA-860 to see if the plant exists in the proposed category
    if len(cems[cems["energy_source_code"].isna()]) > 0:
        gen_fuel = load_data.load_pudl_table(
            "core_eia860__scd_generators",
            year,
            columns=["plant_id_eia", "generator_id", "energy_source_code_1"],
        ).drop_duplicates()
        logger.info("Adding subplant_id to gen_fuel for cems fuel assignment")
        gen_fuel = add_subplant_ids_to_df(
            gen_fuel, year, "generator_id", how_merge="inner", validate_merge="1:1"
        )

        # make sure there are no duplicate unit entries
        gen_fuel = gen_fuel.drop_duplicates(subset=["plant_id_eia", "subplant_id"])
        cems = cems.merge(
            gen_fuel[["plant_id_eia", "subplant_id", "energy_source_code_1"]],
            how="left",
            on=["plant_id_eia", "subplant_id"],
            validate="m:1",
        )
        cems = fillna_with_missing_strings(
            cems,
            column_to_fill="energy_source_code",
            filler_column="energy_source_code_1",
        )

    # if we are still missing fuel codes, merge in from the epa-assigned fuel code
    if (len(cems[cems["energy_source_code"].isna()]) > 0) and (
        "emissions_unit_id_epa" in cems.columns
    ):
        crosswalk = load_data.load_epa_eia_crosswalk_from_raw(year)[
            ["plant_id_eia", "emissions_unit_id_epa", "energy_source_code_epa"]
        ].drop_duplicates(subset=["plant_id_eia", "emissions_unit_id_epa"])
        cems = cems.merge(
            crosswalk,
            how="left",
            on=["plant_id_eia", "emissions_unit_id_epa"],
            validate="m:1",
        )
        cems = fillna_with_missing_strings(
            cems,
            column_to_fill="energy_source_code",
            filler_column="energy_source_code_epa",
        )

    # if we are still missing fuel codes, use energy_source_code_1 of generator with
    # greatest nameplate capacity
    if len(cems[cems["energy_source_code"].isna()]) > 0:
        plant_backstop_fuel = load_backstop_energy_source_codes_for_plant()
        cems = cems.merge(
            plant_backstop_fuel,
            how="left",
            on=["plant_id_eia"],
            validate="m:1",
        )
        cems = fillna_with_missing_strings(
            cems,
            column_to_fill="energy_source_code",
            filler_column="energy_source_code_plant",
        )

    # update
    cems = update_energy_source_codes(cems, year)

    return cems


def load_backstop_energy_source_codes_for_plant() -> pd.DataFrame:
    """Assign an energy_source_code to plant using energy_source_code_1 of generator
    with the greatest nameplate capacity. This backstop was primarily written to
    address issues where a plant no longer exists in EIA due to retirement but
    continues to report data to CEMS, and/or has an emissions_unit_id that is not
    mapped to a generator_id that exists in EIA.

    Returns:
        pd.DataFrame: a data frame relating a plant to an energy source code.
    """

    gens = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        columns=[
            "plant_id_eia",
            "generator_id",
            "report_date",
            "capacity_mw",
            "energy_source_code_1",
        ],
    )
    # keep the set of records for the most recent year for which data is available for
    # that plant
    plant_backstop_fuel = gens[
        (
            gens["report_date"].dt.year
            == gens.groupby("plant_id_eia")["report_date"].transform("max").dt.year
        )
    ]
    # calculate the total capacity associated with each ESC at each plant
    plant_backstop_fuel = (
        plant_backstop_fuel.groupby(["plant_id_eia", "energy_source_code_1"])[
            "capacity_mw"
        ]
        .sum()
        .reset_index()
    )
    # identify the ESC associated with the greatest amount of capacity
    plant_backstop_fuel = plant_backstop_fuel[
        plant_backstop_fuel["capacity_mw"]
        == plant_backstop_fuel.groupby("plant_id_eia")["capacity_mw"].transform("max")
    ]
    # prepare the table for merging by renaming the column and dropping any duplicate
    # primary fuels
    plant_backstop_fuel = plant_backstop_fuel.drop_duplicates(
        subset=["plant_id_eia"], keep=False
    )
    plant_backstop_fuel = plant_backstop_fuel.rename(
        columns={"energy_source_code_1": "energy_source_code_plant"}
    ).drop(columns="capacity_mw")

    return plant_backstop_fuel


def inventory_input_data_sources(cems: pd.DataFrame, year: int):
    """
    Exports a csv that identifies for each plant-month whether input data exists from
    CEMS or EIA-923. This will be used when checking for missing timestamps in the
    output data.

    This function expects CEMS data after it has been loaded, non-grid connected plants
    removed, report date added, and incomplete unit-months removed.
    """

    # load EIA-923 generation and fuel data
    gf = load_data.load_pudl_table("out_eia923__generation_fuel_combined", year)

    # sum generation and fuel by plant-month
    plant_months_in_eia = (
        gf.groupby(["plant_id_eia", "report_date"], dropna=False)[
            ["net_generation_mwh", "fuel_consumed_mmbtu"]
        ]
        .sum(min_count=1)
        .reset_index()
    )
    # add a flag if the data exists
    plant_months_in_eia["data_in_eia"] = 1
    # add a flag if the data is nonzero
    plant_months_in_eia["nonzero_data_in_eia"] = 0
    plant_months_in_eia.loc[
        plant_months_in_eia[["net_generation_mwh", "fuel_consumed_mmbtu"]].sum(axis=1)
        > 0,
        "nonzero_data_in_eia",
    ] = 1

    plant_months_in_eia = plant_months_in_eia[
        ["plant_id_eia", "report_date", "data_in_eia", "nonzero_data_in_eia"]
    ].sort_values(by=["plant_id_eia", "report_date"])

    # sum the generation and fuel data for each plant-month,
    # if there is at least one non-na value
    plant_months_in_cems = (
        cems.groupby(["plant_id_eia", "report_date"], dropna=False)[
            ["gross_generation_mwh", "fuel_consumed_mmbtu"]
        ]
        .sum(min_count=1)
        .reset_index()
    )
    # add an indicator column if any data exists
    plant_months_in_cems["data_in_cems"] = 0
    plant_months_in_cems.loc[
        (
            ~plant_months_in_cems[["gross_generation_mwh", "fuel_consumed_mmbtu"]]
            .isna()
            .all(axis=1)
        ),
        "data_in_cems",
    ] = 1

    # add a flag if the data is zero
    plant_months_in_cems["nonzero_data_in_cems"] = 0
    plant_months_in_cems.loc[
        plant_months_in_cems[["gross_generation_mwh", "fuel_consumed_mmbtu"]].sum(
            axis=1
        )
        > 0,
        "nonzero_data_in_cems",
    ] = 1

    plant_months_in_cems = plant_months_in_cems[
        ["plant_id_eia", "report_date", "data_in_cems", "nonzero_data_in_cems"]
    ].sort_values(by=["plant_id_eia", "report_date"])

    # outer merge cems and eia inventories together
    input_data_exists = plant_months_in_eia.merge(
        plant_months_in_cems,
        how="outer",
        on=["plant_id_eia", "report_date"],
        validate="1:1",
    ).sort_values(by=["plant_id_eia", "report_date"])

    # create a dataframe with a complete set of plant-months in case some months
    # are missing from the inventory
    complete_plant_months = pd.DataFrame(
        list(
            product(
                list(input_data_exists.plant_id_eia.unique()),
                list(input_data_exists.report_date.unique()),
            )
        ),
        columns=["plant_id_eia", "report_date"],
    )
    # make sure the datetime dtypes match before merging
    complete_plant_months["report_date"] = complete_plant_months.report_date.astype(
        input_data_exists.report_date.dtype
    )
    # complete the report dates
    input_data_exists = input_data_exists.merge(
        complete_plant_months,
        how="outer",
        on=["plant_id_eia", "report_date"],
        validate="1:1",
    ).sort_values(by=["plant_id_eia", "report_date"])

    # fill missing values with zeros
    input_data_exists["data_in_eia"] = input_data_exists["data_in_eia"].fillna(0)
    input_data_exists["data_in_cems"] = input_data_exists["data_in_cems"].fillna(0)
    input_data_exists["nonzero_data_in_eia"] = input_data_exists[
        "nonzero_data_in_eia"
    ].fillna(0)
    input_data_exists["nonzero_data_in_cems"] = input_data_exists[
        "nonzero_data_in_cems"
    ].fillna(0)

    # create a column to indicate if data exists from any source
    input_data_exists["input_data_exists"] = 0
    input_data_exists.loc[
        input_data_exists[["data_in_eia", "data_in_cems"]].sum(axis=1) > 0,
        "input_data_exists",
    ] = 1
    # create a column to indicate if nonzero data exists from any source
    input_data_exists["nonzero_input_data_exists"] = 0
    input_data_exists.loc[
        input_data_exists[["nonzero_data_in_eia", "nonzero_data_in_cems"]].sum(axis=1)
        > 0,
        "nonzero_input_data_exists",
    ] = 1

    # export to csv
    input_data_exists.to_csv(
        outputs_folder(f"{year}/input_data_inventory_{year}.csv.zip"),
        index=False,
        compression="zip",
    )


def fillna_with_missing_strings(df, column_to_fill, filler_column):
    """
    Fills missing values in string columns using another string column with missing values.

    If the column you are using to fill missing values in a string column also has missing values,
    pandas raises a `ValueError: Must provide strings.` To get around this, we must convert NA
    values to strings, fillna, then convert back to NA.

    Args:
        df: dataframe containing column_to_fill and filler_column
        column_to_fill: name of the column that you want to fill
        filler_column: name of the column from which you will fill values
    """
    # replace missing values with a string
    df[filler_column] = df[filler_column].fillna("MISSING")
    # fill missing values
    df[column_to_fill] = df[column_to_fill].fillna(df[filler_column])
    # drop the filler column in-place
    df.drop(columns=[filler_column], inplace=True)
    # convert the missing string back into a missing value
    df[column_to_fill] = df[column_to_fill].replace("MISSING", np.NaN)

    return df


def fill_missing_fuel_for_single_fuel_plant_months(df, year):
    """
    Identifies all plant-months where a single fuel was burned based on the EIA-923 generation fuel table
    Uses this to fill in the energy source code if a match was not made based on the PSDC
    """

    # identify plant-months for which there is a single fossil fuel type reported
    gf = load_data.load_pudl_table(
        "out_eia923__generation_fuel_combined",
        year,
        columns=[
            "plant_id_eia",
            "report_date",
            "energy_source_code",
            "fuel_consumed_mmbtu",
        ],
    )

    # remove any rows for clean fuels
    gf = gf[~gf["energy_source_code"].isin(CLEAN_FUELS)]

    # group the data by plant, month, and fuel
    gf = (
        gf.groupby(["plant_id_eia", "report_date", "energy_source_code"], dropna=False)
        .sum()
        .reset_index()
    )

    # only keep rows with >0 fuel consumed reported
    gf = gf[gf["fuel_consumed_mmbtu"] > 0]

    # identify which plant months have multiple fuels reported
    multi_fuels = (
        gf.groupby(["plant_id_eia", "report_date"], dropna=False)["energy_source_code"]
        .count()
        .reset_index()
        .rename(columns={"energy_source_code": "num_fuels"})
    )

    # merge this information back into the other dataframe
    gf = gf.merge(
        multi_fuels, how="left", on=["plant_id_eia", "report_date"], validate="m:1"
    )

    # only keep rows that have a single fuel type in a month
    gf = gf[gf["num_fuels"] == 1]

    # clean up the columns
    gf = gf.rename(columns={"energy_source_code": "energy_source_code_single"})
    gf.drop(columns=["num_fuels", "fuel_consumed_mmbtu"], inplace=True)
    gf["report_date"] = pd.to_datetime(gf["report_date"]).astype("datetime64[s]")

    # merge this data into the df
    df = df.merge(gf, how="left", on=["plant_id_eia", "report_date"], validate="m:1")

    df = fillna_with_missing_strings(
        df,
        column_to_fill="energy_source_code",
        filler_column="energy_source_code_single",
    )

    return df


def remove_cems_with_zero_monthly_data(cems, column_to_check: list[str]):
    """
    Identifies months where zero generation or heat inputare reported.
    from each unit and removes associated hours from CEMS so that these can be filled using the eia923 data
    Inputs:
        cems: pandas dataframe of hourly cems data containing columns "plant_id_eia", "emissions_unit_id_epa" and "report_date"
    Returns:
        cems df with hourly observations for months when no emissions reported removed
    """
    # calculate the totals reported in each month
    # NOTE: columns_to_check = ["gross_generation_mwh", "fuel_consumed_mmbtu"]
    data_with_zero_monthly_values = cems.groupby(
        ["plant_id_eia", "subplant_id", "report_date"], dropna=False
    )[column_to_check].sum()
    # identify unit-months where zero emissions reported
    data_with_zero_monthly_values = data_with_zero_monthly_values[
        data_with_zero_monthly_values.sum(axis=1) == 0
    ]
    # add a flag to these observations
    data_with_zero_monthly_values["zero_data_flag"] = "remove"

    # merge the missing data flag into the cems data
    cems = cems.merge(
        data_with_zero_monthly_values.reset_index()[
            [
                "plant_id_eia",
                "subplant_id",
                "report_date",
                "zero_data_flag",
            ]
        ],
        how="left",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="m:1",
    )
    # get a count of the number of observations with the missing data flag
    num_observations_with_zero_data = len(cems[cems["zero_data_flag"] == "remove"])
    # get a count of the total number of zero observations in the data
    total_zero_observations = len(cems[cems[column_to_check].sum(axis=1) == 0])
    # remove any observations with the missing data flag
    logger.info(
        f"{total_zero_observations / len(cems) * 100:.2f}% of rows in the data have zero values, and {num_observations_with_zero_data / total_zero_observations * 100:.2f}% of these represent complete months of zero data"
    )

    validation.check_removed_data_is_empty(cems)
    pre_memory_usage_gb = cems.memory_usage().sum() / 1_000_000_000
    # cems = cems[cems["zero_data_flag"] != "remove"]
    # remove all zero obeservations
    cems = cems[~(cems[column_to_check].sum(axis=1) == 0)]
    post_memory_usage_gb = cems.memory_usage().sum() / 1_000_000_000
    logger.info(
        f"Memory usage after removing zero data: {post_memory_usage_gb:.2f} GB ({(post_memory_usage_gb / pre_memory_usage_gb * 100):.2f}% of original)"
    )
    # drop the missing data flag column in-place
    cems.drop(columns="zero_data_flag", inplace=True)

    return cems


def complete_hourly_timeseries(
    df: pd.DataFrame,
    year: int,
    group_cols: list[str] = [],
    columns_to_fill_with_zero: list[str] = [],
    columns_to_bffill: list[str] = [],
) -> pd.DataFrame:
    """
    Completes a hourly timeseries for each plantgroup in a given dataframe.


    This function will create a complete timeseries for each group in the dataframe,
    assuming the data is supposed to represent a complete year. It will create a complete
    timeseries for each group in UTC time.

    Because we are repairing the "datetime_utc" column, but the complete set of utc
    timestamps for a year depends on the local timezone of the underlying data, we
    have to consider the timezone of each plant_id_eia in the data, and assign the
    appropriate complete timeseries to each. This approach is more robust than just
    creating a complete date_range between the min and max datetimes already in the
    timeseries, in case that timeseries is missing timestamps at the beginning or end
    of the year.

    If "report_date" is passed in as one of the `group_cols`, the behavior of this
    function is that it will only complete timeseries for "report_date"s that already
    exist in `df`. For example, if `df` only contained data for January-June, if
    "report_date" is included, then this function will only repair hourly timeseries for
    January-June, but will not add timestamps for July-December. If that same df were
    passed in without "report_date" in the `group_cols`, the function would also add
    hourly timestamps for July-December. This functionality exists so that when we repair
    `combined_cems_subplant_data` before combining it with the shaped eia data in step
    18, we don't end up with overlapping timeseries.

    Args:
        df: dataframe to complete the timeseries for
        year: year to create a complete timeseries for (should match the local year of the data)
        group_cols: columns to group by (one column must be "plant_id_eia")
        columns_to_fill_with_zero: a list of columns (generally containing numeric data)
            that should be filled with zero values for missing timestamps that are filled
        columns_to_bffill: a list of columns that should be filled based on the values in the previous and next rows within each group (bbffill refers to both bfill and ffill)
            in the previous and next rows within each group (bbffill refers to both bfill
            and ffill)
    Returns:
        dataframe with complete timeseries
    """

    # check if there are any missing timestamps
    expected_hours = 8784 if (year % 4 == 0) else 8760
    test = df.groupby(group_cols)[["datetime_utc"]].count()
    # only repair if it is needed, otherwise, skip this
    if len(test[test["datetime_utc"] < expected_hours]) > 0:
        # get all unique groups for which to create complete timeseries
        complete_timeseries = df[group_cols].drop_duplicates()

        # merge in timezone data for each plant
        plant_timezone = load_data.load_pudl_table(
            "core_eia__entity_plants", columns=["plant_id_eia", "timezone"]
        )
        complete_timeseries = complete_timeseries.merge(
            plant_timezone, on="plant_id_eia", how="left"
        )

        # localize the datetime_local column using the timezone column
        # get a list of timezones to iterate through. Becuase pandas has trouble with
        # tz localization and conversion if multiple tz's exist in a single column, we
        # need to do these operations one at a time for each timezone, then concat them
        timezones = complete_timeseries["timezone"].unique()
        timeseries_for_timezones = []
        for timezone in timezones:
            # first get the complete set of timestamps in local time
            timezone_df = pd.DataFrame(
                data=pd.date_range(
                    start=f"{year}-01-01 00:00:00",
                    end=f"{year}-12-31 23:00:00",
                    freq="h",
                    tz=timezone,
                    name="datetime_local",
                )
            )
            # now convert to UTC
            timezone_df["datetime_utc"] = timezone_df["datetime_local"].dt.tz_convert(
                "UTC"
            )
            timezone_df["report_date"] = (
                timezone_df["datetime_local"].dt.to_period("M").dt.to_timestamp()
            )
            timezone_df = timezone_df.drop(columns=["datetime_local"])
            timezone_df["timezone"] = timezone
            timeseries_for_timezones.append(timezone_df)
        timeseries_for_timezones = pd.concat(timeseries_for_timezones)

        # merge the complete timeseries for each timezone into the groups to create
        # complete timeseries for each group
        if "report_date" in group_cols:
            complete_timeseries = complete_timeseries.merge(
                timeseries_for_timezones,
                on=["timezone", "report_date"],
                how="left",
                validation="m:m",
            )
        else:
            complete_timeseries = complete_timeseries.merge(
                timeseries_for_timezones.drop(columns=["report_date"]),
                on=["timezone"],
                how="left",
                validation="m:m",
            )
        complete_timeseries = complete_timeseries.drop(columns=["timezone"])

        # complete the timeseries in the original dataframe
        pre_completion_size = len(df)
        df = df.merge(
            complete_timeseries, on=(group_cols + ["datetime_utc"]), how="outer"
        ).sort_values(by=group_cols + ["datetime_utc"], ascending=True)
        post_completion_size = len(df)

        if post_completion_size != pre_completion_size:
            logger.info(
                f"complete_hourly_timeseries() added {post_completion_size - pre_completion_size} missing rows to the dataframe"
            )

        # fill values for missing timestamps
        df[columns_to_fill_with_zero] = df[columns_to_fill_with_zero].fillna(0.0)

        # forward and backfill columns within each group. This can be used for
        # non-numeric columns
        df[columns_to_bffill] = (
            df.groupby(group_cols)[columns_to_bffill].ffill().bfill()
        )

        return df
    else:
        return df


def adjust_cems_for_chp(cems, eia923_allocated):
    """
    Adjusts CEMS fuel consumption and emissions data for CHP.

    Steps:
        1. Calculate the ratio between `fuel_consumed_for_electricity_mmbtu` and `fuel_consumed_mmbtu` in EIA-923
        2. Use this ratio to calculate a `fuel_consumed_for_electricity_mmbtu` from the `fuel_consumed_mmbtu` data reported in CEMS
        3. Calculate an electric allocation factor using the fuel and net generation data
        4. Use the allocation factor to adjust emissions and fuel consumption
    Args:
        cems: dataframe of hourly cems data after cleaning and gross to net calculations
        eia923_allocated: dataframe of EIA-923 data after allocation
    """
    # calculate a subplant fuel ratio
    subplant_fuel_ratio = (
        eia923_allocated.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )[["fuel_consumed_mmbtu", "fuel_consumed_for_electricity_mmbtu"]]
        .sum()
        .reset_index()
    )
    subplant_fuel_ratio["subplant_fuel_ratio"] = (
        subplant_fuel_ratio["fuel_consumed_for_electricity_mmbtu"]
        / subplant_fuel_ratio["fuel_consumed_mmbtu"]
    )
    subplant_fuel_ratio.loc[
        (subplant_fuel_ratio["fuel_consumed_for_electricity_mmbtu"] == 0)
        & (subplant_fuel_ratio["fuel_consumed_mmbtu"] == 0),
        "subplant_fuel_ratio",
    ] = 1
    # calculate a plant fuel ratio to fill missing values where there is not a matching subplant in CEMS
    plant_fuel_ratio = (
        eia923_allocated.groupby(["plant_id_eia", "report_date"], dropna=False)[
            ["fuel_consumed_mmbtu", "fuel_consumed_for_electricity_mmbtu"]
        ]
        .sum()
        .reset_index()
    )
    plant_fuel_ratio["plant_fuel_ratio"] = (
        plant_fuel_ratio["fuel_consumed_for_electricity_mmbtu"]
        / plant_fuel_ratio["fuel_consumed_mmbtu"]
    )
    plant_fuel_ratio.loc[
        (plant_fuel_ratio["fuel_consumed_for_electricity_mmbtu"] == 0)
        & (plant_fuel_ratio["fuel_consumed_mmbtu"] == 0),
        "plant_fuel_ratio",
    ] = 1

    # merge the fuel ratios into cems and fill missing subplant ratios with plant ratios
    cems = cems.merge(
        subplant_fuel_ratio[
            ["plant_id_eia", "subplant_id", "report_date", "subplant_fuel_ratio"]
        ],
        how="left",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="m:1",
    )
    cems = cems.merge(
        plant_fuel_ratio[["plant_id_eia", "report_date", "plant_fuel_ratio"]],
        how="left",
        on=["plant_id_eia", "report_date"],
        validate="m:1",
    )
    cems["subplant_fuel_ratio"] = cems["subplant_fuel_ratio"].fillna(
        cems["plant_fuel_ratio"]
    )

    # if there are any missing ratios, assume that the ratio is 1
    cems["subplant_fuel_ratio"] = cems["subplant_fuel_ratio"].fillna(1)

    # calculate fuel_consumed_for_electricity_mmbtu
    cems["fuel_consumed_for_electricity_mmbtu"] = (
        cems["fuel_consumed_mmbtu"] * cems["subplant_fuel_ratio"]
    )

    # remove intermediate columns in-place
    cems.drop(columns=["subplant_fuel_ratio", "plant_fuel_ratio"], inplace=True)

    # add adjusted emissions columns
    cems = emissions.adjust_fuel_and_emissions_for_chp(cems)

    return cems


def identify_hourly_data_source(eia923_allocated, cems, year):
    """Identifies whether there is hourly CEMS data available for each subplant-month.
    Possible categories:
        1. `cems`: For subplant-months for which we have hourly CEMS data for all CEMS
            units that make up that subplant, we will use the hourly values reported in
            CEMS.
        2. `partial_cems_subplant`: For subplant-months for which we have hourly CEMS data
            for only some of the CEMS units that make up a subplant, we will use the reported
            EIA-923 values to scale the partial hourly CEMS data from the other units to
            match the total value for the entire subplant. We will also calculate a partial
            subplant scaling factor for each data column (e.g. net generation, fuel
            consumption) by comparing the total monthly CEMS data to the monthly EIA-923 data.
        3. `partial_cems_plant`: when some subplants in a plant have CEMS data, use the shape
            from those subplants to shape the remaining generation from that plant.
            NOTE: This method only applies if the subplants missing CEMS data have a primary
            fuel that would report to CEMS (ie no clean generators). This prevents a diesel
            backup generator that reports to CEMS from being used to shape the generation of
            a nuclear plant (for example) that does not report to CEMS.
        4. `eia`: for subplant-months for which no hourly data is reported in CEMS, we will
            attempt to use EIA-930 data to assign an hourly profile to the monthly EIA-923 data
    Inputs:
        eia923_allocated:
        cems:
        year:
    Returns:
        eia923_allocated with new column `hourly_data_source`
    """

    # Note: We copy here to avoid modifying the original dataframe passed by the caller
    all_data = eia923_allocated.copy()

    # create a binary column indicating whether any data was reported in 923
    columns_to_test = [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
    ]
    all_data = all_data.assign(
        reported_eia923=lambda x: np.where(
            x[columns_to_test].notnull().all(axis=1), 1, 0
        )
    )

    # identify which cems data only represents part of a subplant
    cems_status = identify_partial_cems_subplants(year, cems, eia923_allocated)

    # merge in the data source column from CEMS
    all_data = all_data.merge(
        cems_status,
        how="left",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="m:1",
    )

    # for the remaining plants, identify the hourly data source as EIA
    all_data["hourly_data_source"] = all_data["hourly_data_source"].fillna("eia")
    all_data.drop(columns=["reported_eia923"], inplace=True)

    # identify the partial cems plants
    all_data = identify_partial_cems_plants(all_data)

    return all_data


def identify_partial_cems_subplants(year, cems, eia923_allocated):
    """Identifies subplants for which data for only some units is reported in CEMS.

    If the number of unique emissions_unit_id_epa reported in cems for a subplant is less than the
    number of units that make up the subplant per the subplant crosswalk, and the
    amount of fuel reported for that subplant in CEMS is less than 95% of the fuel
    reported for that subplant in EIA, then we mark the subplant-month as "partial cems"

    This means that we will only use the cems data for this subplant month to shape the
    EIA data.
    """
    units_in_subplant = count_total_units_in_subplant(year)

    # aggregate cems data to plant-unit-month
    cems_unit_month_agg = (
        cems.groupby(
            ["plant_id_eia", "subplant_id", "emissions_unit_id_epa", "report_date"],
            dropna=False,
        )[["fuel_consumed_mmbtu"]]
        .sum()
        .reset_index()
    )

    # create a dataframe that counts the number of units reported in CEMS in each subplant-month
    cems_units_reported = count_reported_units_in_subplant(cems_unit_month_agg)

    # merge in the total number of units that exist in each subplant
    # this will allow us to compare where a subplant-month is missing data from one or more units
    cems_units_reported = cems_units_reported.merge(
        units_in_subplant,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )

    # aggregate fuel reported in EIA by subplant month and merge into our cems units reported
    # This will allow us to compare whether reported fuel consumption is significantly less
    eia_subplant_month_agg = (
        eia923_allocated.groupby(
            ["report_date", "plant_id_eia", "subplant_id"], dropna=False
        )[["fuel_consumed_mmbtu"]]
        .sum(min_count=1)
        .reset_index()
    )
    cems_units_reported = cems_units_reported.merge(
        eia_subplant_month_agg,
        how="left",
        on=["plant_id_eia", "subplant_id", "report_date"],
        suffixes=("_cems", "_eia"),
        validate="1:1",
    )

    # identify which subplant-months have complete or partial cems data
    # we identify this as subplant months where the number of CEMS reported units is
    # less than the total number of units in that subplant, AND the total, non-zero fuel
    # reported in CEMS is less than 95% of the non-zero fuel reported in EIA
    cems_status = cems_units_reported.assign(
        hourly_data_source=lambda x: np.where(
            (
                (x.reported_units_in_subplant < x.units_in_subplant)
                & (x.fuel_consumed_mmbtu_cems < (x.fuel_consumed_mmbtu_eia * 0.95))
                & (x.fuel_consumed_mmbtu_cems > 0)
                & (x.fuel_consumed_mmbtu_eia > 0)
            ),
            "partial_cems_subplant",
            "cems",
        )
    )

    cems_status = cems_status[
        ["plant_id_eia", "subplant_id", "report_date", "hourly_data_source"]
    ]

    return cems_status


def count_total_units_in_subplant(year):
    # load the subplant crosswalk and identify unique emissions_unit_id_epas in each subplant
    units_in_subplant = (
        pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv.zip"),
            dtype=get_dtypes(),
            parse_dates=[
                "current_planned_generator_operating_date",
                "generator_retirement_date",
            ],
        )[
            [
                "plant_id_eia",
                "emissions_unit_id_epa",
                "subplant_id",
                "generator_retirement_date",
            ]
        ]
        .drop_duplicates()
        .dropna(subset="emissions_unit_id_epa")
    )

    # remove units that retired before the current year
    units_in_subplant = units_in_subplant[
        ~(units_in_subplant["generator_retirement_date"].dt.year < year)
    ]

    # get a count of the number of CEMS units in each subplant
    units_in_subplant = (
        units_in_subplant.groupby(["plant_id_eia", "subplant_id"], dropna=False)
        .count()["emissions_unit_id_epa"]
        .reset_index()
        .rename(columns={"emissions_unit_id_epa": "units_in_subplant"})
    )

    return units_in_subplant


def count_reported_units_in_subplant(cems_monthly):
    """Counts the number of units in each subplant-month in CEMS.

    Also sums the total fuel reported in each subplant-month to compare with EIA
    """

    reported_units_in_subplant = (
        cems_monthly.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .agg({"emissions_unit_id_epa": "count", "fuel_consumed_mmbtu": "sum"})
        .reset_index()
        .rename(columns={"emissions_unit_id_epa": "reported_units_in_subplant"})
    )

    return reported_units_in_subplant


def identify_partial_cems_plants(all_data):
    """Identifies subplants for which CEMS data is reported for at least one other subplant at the plant.

    Args:
        all_data: dataframe identifying the hourly data source for each subplant-month
    Returns:
        all_data with updated hourly_data_source column indicating partial cems plants
    """

    # create a column that indicates EIA-only subplant data
    all_data = all_data.assign(
        eia_data=lambda x: np.where(x.hourly_data_source == "eia", 1, 0)
    )

    # create a column that indicates whether cems data is reported
    all_data = all_data.assign(
        cems_data=lambda x: np.where(
            x.hourly_data_source.isin(["cems", "partial_cems_subplant"]), 1, 0
        )
    )

    # count the number of records for each plant-month that have EIA-only data and CEMS data
    partial_plant = (
        all_data.groupby(["plant_id_eia", "report_date"], dropna=False)[
            ["eia_data", "cems_data"]
        ]
        .sum()
        .reset_index()
    )

    # we want to identify plants where there is some EIA data and some cems data
    # if eia_data = 0, then all of the data for that plant-month is from cems
    # if if cems_data = 0, then all of the data for that plant-month is from EIA
    partial_plant = partial_plant[
        (partial_plant["eia_data"] > 0) & (partial_plant["cems_data"] > 0)
    ]

    # drop intermediate columns in-place
    partial_plant.drop(columns=["eia_data", "cems_data"], inplace=True)

    # merge this data into all_data
    all_data = all_data.merge(
        partial_plant,
        how="outer",
        on=["plant_id_eia", "report_date"],
        indicator="partial_plant",
        validate="m:1",
    )

    # Plants with these primary fuel types shouldn't be shaped using CEMS data --
    # any cems data is likely to be a backup generator and therefore not representative of overall generation.
    non_cems_fuels = CLEAN_FUELS + ["GEO"]

    # for plants where the data source is eia only, and it exists in partial plant,
    # change the source to partial_cems_plant
    # but only if the fuel type is consistent with CEMS fuel types
    all_data.loc[
        (all_data["hourly_data_source"] == "eia")
        & (all_data["partial_plant"] == "both")
        & (
            ~all_data["plant_primary_fuel"].isin(non_cems_fuels)
        ),  # CEMS data should only be used to shape if plant is primarilly fossil
        "hourly_data_source",
    ] = "partial_cems_plant"

    # It is possible to have generators within a subplant assigned different hourly methods
    # due to the source code check above if subplant_id is mixing fuel types.
    # This was an issue before PR #239; this check will catch it if it reoccurs
    mixed_method_subplants = (
        all_data.groupby(
            [
                "plant_id_eia",
                "subplant_id",
                "report_date",
            ],
            dropna=False,
        )
        .nunique()
        .reset_index()
    )
    mixed_method_subplants = mixed_method_subplants[
        mixed_method_subplants.hourly_data_source > 1
    ]
    if len(mixed_method_subplants) > 0:
        # Check for subplants with mixed hourly data sources,
        # likely resulting from mixed fuel types.
        # If subplant_id assignment is working, there shouldn't be any
        print(mixed_method_subplants)
        raise Exception(
            f"ERROR: {len(mixed_method_subplants)} subplant-months have multiple hourly methods assigned."
        )

    # remove the intermediate indicator columns in-place
    all_data.drop(
        columns=[
            "partial_plant",
            "eia_data",
            "cems_data",
        ],
        inplace=True,
    )

    return all_data


def filter_unique_cems_data(cems, partial_cems):
    """Removes subplant-months from cems that also appear in partial_cems."""
    # filter the cems data to remove data that was scaled in partial_cems
    partial_cems_subplant_months = partial_cems[
        ["plant_id_eia", "subplant_id", "report_date"]
    ].drop_duplicates()
    filtered_cems = cems.merge(
        partial_cems_subplant_months,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        indicator="source",
        validate="m:1",
    )

    filtered_cems = filtered_cems[filtered_cems["source"] == "left_only"]
    filtered_cems.drop(columns=["source"], inplace=True)

    return filtered_cems


def aggregate_subplant_data_to_fleet(
    combined_subplant_data: pd.DataFrame,
    plant_attributes_table: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Group plant data by BA and fuel category (fleet).

    Args:
        combined_subplant_data (pd.DataFrame): the combined subplant data.
        plant_attributes_table (pd.DataFrame): the plant attributes table with the
            BA code of each plant.
        primary_fuel_table (pd.DataFrame): table with subplant and plant-level fuels
        year (int): the data year

    Raises:
        UserWarning: if no BA or fuel category can be assigned to any subplant

    Returns:
        pd.DataFrame: a data frame grouped by fuel category and BA
    """

    # Assign fleet to subplant data
    ba_fuel_data = assign_fleet_to_subplant_data(
        combined_subplant_data,
        plant_attributes_table,
        primary_fuel_table,
        year,
        ba_col="ba_code",
        primary_fuel_col="subplant_primary_fuel",
        fuel_category_col="fuel_category",
    )

    # drop subplants that have missing fuel category and no generation or fuel data
    # this prevents them from creating blank entries in the power sector results data
    ba_fuel_data = ba_fuel_data[
        ~(
            ba_fuel_data["fuel_category"].isna()
            & (ba_fuel_data["net_generation_mwh"] == 0)
            & (ba_fuel_data["fuel_consumed_for_electricity_mmbtu"] == 0)
        )
    ]

    # if the input data is hourly, aggregate at the hourly level
    if "datetime_utc" in ba_fuel_data.columns:
        agg_cols = ["ba_code", "fuel_category", "datetime_utc", "report_date"]
    else:
        agg_cols = ["ba_code", "fuel_category", "report_date"]
    ba_fuel_data = (
        ba_fuel_data.groupby(
            agg_cols,
            dropna=False,
        )[DATA_COLUMNS]
        .sum()
        .reset_index()
    )

    return ba_fuel_data


def aggregate_cems_to_subplant(cems):
    GROUPBY_COLUMNS = ["plant_id_eia", "subplant_id", "datetime_utc", "report_date"]

    cems_columns_to_aggregate = [
        "gross_generation_mwh",
        # "steam_load_1000_lb",
        "fuel_consumed_mmbtu",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
    ]

    cems = (
        cems.groupby(GROUPBY_COLUMNS, dropna=False)[cems_columns_to_aggregate]
        .sum()
        .reset_index()
    )

    cems = apply_dtypes(cems)

    return cems
