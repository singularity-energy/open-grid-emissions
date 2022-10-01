import pandas as pd
import numpy as np
import os
import sqlalchemy as sa

import pudl.analysis.allocate_net_gen as allocate_gen_fuel
import pudl.analysis.epa_crosswalk as epa_crosswalk
import pudl.output.pudltabl

import load_data
import validation
import emissions
from emissions import CLEAN_FUELS
from column_checks import get_dtypes, apply_dtypes
from filepaths import manual_folder, outputs_folder, downloads_folder


DATA_COLUMNS = [
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
]


def identify_subplants(year, number_of_years):
    """This is the coordinating function for loading and calculating subplant IDs, GTN regressions, and GTN ratios."""
    start_year = year - (number_of_years - 1)
    end_year = year

    # load 5 years of monthly data from CEMS
    cems_ids = load_data.load_cems_ids(start_year, end_year)

    # add subplant ids to the data
    generate_subplant_ids(start_year, end_year, cems_ids)


def generate_subplant_ids(start_year, end_year, cems_ids):
    """
    Groups units and generators into unique subplant groups.

    This function consists of three primary parts:
    1. Identify a list of all unique plant-units that exist in the CEMS data
        for the years in question. This will be used to filter the crosswalk.
    2. Load the EPA-EIA crosswalk and filter it based on the units that exist
        in the CEMS data for the years in question
    3. Use graph analysis to identify distinct groupings of EPA units and EIA
        generators based on 1:1, 1:m, m:1, or m:m relationships.

    Returns:
        exports the subplant crosswalk to a csv file
        cems_ids and gen_fuel_allocated with subplant_id added
    """

    # load the crosswalk and filter it by the data that actually exists in cems
    crosswalk = pudl.output.epacems.epa_crosswalk()

    # update the crosswalk with manual matches
    crosswalk = manual_crosswalk_updates(crosswalk)

    # strip leading zeros
    crosswalk["CAMD_UNIT_ID"] = crosswalk["CAMD_UNIT_ID"].str.lstrip("0")
    crosswalk["EIA_GENERATOR_ID"] = crosswalk["EIA_GENERATOR_ID"].str.lstrip("0")

    # filter the crosswalk to drop any units that don't exist in CEMS
    filtered_crosswalk = epa_crosswalk.filter_crosswalk(crosswalk, cems_ids)[
        [
            "plant_id_eia",
            "unitid",
            "CAMD_PLANT_ID",
            "CAMD_UNIT_ID",
            "CAMD_GENERATOR_ID",
            "EIA_PLANT_ID",
            "EIA_GENERATOR_ID",
        ]
    ]

    # change the plant id to an int
    filtered_crosswalk["EIA_PLANT_ID"] = filtered_crosswalk["EIA_PLANT_ID"].astype(int)
    # use graph analysis to identify subplants
    crosswalk_with_subplant_ids = epa_crosswalk.make_subplant_ids(filtered_crosswalk)

    # fix the column names
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.drop(
        columns=["plant_id_eia", "unitid", "CAMD_GENERATOR_ID"]
    )
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.rename(
        columns={
            "CAMD_PLANT_ID": "plant_id_epa",
            "EIA_PLANT_ID": "plant_id_eia",
            "CAMD_UNIT_ID": "unitid",
            "EIA_GENERATOR_ID": "generator_id",
        }
    )
    # change the eia plant id to int
    crosswalk_with_subplant_ids["plant_id_eia"] = crosswalk_with_subplant_ids[
        "plant_id_eia"
    ].astype(int)

    # change the order of the columns
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids[
        ["plant_id_epa", "unitid", "plant_id_eia", "generator_id", "subplant_id"]
    ]

    # update the subplant_crosswalk to ensure completeness
    # prepare the subplant crosswalk by adding a complete list of generators and adding the unit_id_pudl column
    pudl_out = load_data.initialize_pudl_out(year=end_year)
    complete_generator_ids = pudl_out.gens_eia860()[
        ["plant_id_eia", "generator_id", "unit_id_pudl"]
    ].drop_duplicates()
    subplant_crosswalk_complete = crosswalk_with_subplant_ids.merge(
        complete_generator_ids,
        how="outer",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )
    # also add a complete list of cems unitid
    subplant_crosswalk_complete = subplant_crosswalk_complete.merge(
        cems_ids[["plant_id_eia", "unitid"]].drop_duplicates(),
        how="outer",
        on=["plant_id_eia", "unitid"],
        validate="m:1",
    )
    # update the subplant ids for each plant
    subplant_crosswalk_complete = subplant_crosswalk_complete.groupby(
        "plant_id_eia"
    ).apply(update_subplant_ids)

    # remove the intermediate columns created by update_subplant_ids
    subplant_crosswalk_complete["subplant_id"].update(
        subplant_crosswalk_complete["new_subplant"]
    )
    subplant_crosswalk_complete = subplant_crosswalk_complete.reset_index(drop=True)[
        [
            "plant_id_epa",
            "unitid",
            "plant_id_eia",
            "generator_id",
            "subplant_id",
            "unit_id_pudl",
        ]
    ]

    subplant_crosswalk_complete = subplant_crosswalk_complete.drop_duplicates(
        subset=[
            "plant_id_epa",
            "unitid",
            "plant_id_eia",
            "generator_id",
            "subplant_id",
        ],
        keep="last",
    )

    # add proposed operating dates and retirements to the subplant id crosswalk
    subplant_crosswalk_complete = add_operating_and_retirement_dates(
        subplant_crosswalk_complete, start_year, end_year
    )

    os.makedirs(outputs_folder(f"{end_year}"), exist_ok=True)

    # export the crosswalk to csv
    subplant_crosswalk_complete.to_csv(
        outputs_folder(f"{end_year}/subplant_crosswalk_{end_year}.csv"),
        index=False,
    )


def update_subplant_ids(subplant_crosswalk):
    """
    Ensures a complete and accurate subplant_id mapping for all generators.

    NOTE:
        1. This function is a temporary placeholder until the `pudl.analysis.epa_crosswalk` code is updated.
        2. This function is meant to be applied using a .groupby("plant_id_eia").apply() function. This function
        will only properly work when applied to a single plant_id_eia at a time.

    Data Preparation
        Because the existing subplant_id crosswalk was only meant to map CAMD units to EIA generators, it
        is missing a large number of subplant_ids for generators that do not report to CEMS. Before applying this
        function to the subplant crosswalk, the crosswalk must be completed with all generators by outer
        merging in the complete list of generators from EIA-860 (specifically the gens_eia860 table from pudl).
        This dataframe also contains the complete list of `unit_id_pudl` mappings that will be necessary.

    This function follows several steps:
        1. Because the current subplant_id code does not take boiler-generator associations into account,
        there may be instances where the code assigns generators to different subplants when in fact, according
        to the boiler-generator association table, these generators are grouped into a single unit based on their
        boiler associations. The first step of this function is thus to identify if multiple subplant_id have
        been assigned to a single unit_id_pudl. If so, we replace the existing subplant_ids with a single subplant_id.
        For example, if a generator A was assigned subplant_id 0 and generator B was assigned subplant_id 1, but
        both generators A and B are part of unit_id_pudl 1, we would re-assign the subplant_id to both generators to
        0 (we always use the lowest number subplant_id in each unit_id_pudl group). This may result in some subplant_id
        being skipped, but this is okay because we will later renumber all subplant ids (i.e. if there were also a
        generator C with subplant_id 2, there would no be no subplant_id 1 at the plant)
        3. The second issue is that there are many NA subplant_id that we should fill. To do this, we first look at
        unit_id_pudl. If a group of generators are assigned a unit_id_pudl but have NA subplant_ids, we assign a single
        new subplant_id to this group of generators. If there are still generators at a plant that have both NA subplant_id
        and NA unit_id_pudl, we for now assume that each of these generators consitutes its own subplant. We thus assign a unique
        subplant_id to each generator that is unique from any existing subplant_id already at the plant.

        Args:
            subplant_crosswalk: a dataframe containing the output of `epa_crosswalk.make_subplant_ids` with
    """
    # Step 1: Create corrected versions of subplant_id and unit_id_pudl
    # if multiple unit_id_pudl are connected by a single subplant_id, unit_id_pudl_connected groups these unit_id_pudl together
    subplant_crosswalk = connect_ids(
        subplant_crosswalk, id_to_update="unit_id_pudl", connecting_id="subplant_id"
    )
    # if multiple subplant_id are connected by a single unit_id_pudl, group these subplant_id together
    subplant_crosswalk = connect_ids(
        subplant_crosswalk, id_to_update="subplant_id", connecting_id="unit_id_pudl"
    )

    # Step 2: Fill missing subplant_id
    # We will use unit_id_pudl to fill missing subplant ids, so first we need to fill any missing unit_id_pudl
    # We do this by assigning a new unit_id_pudl to each generator that isn't already grouped into a unit

    # create a numeric version of each generator_id
    # ngroup() creates a unique number for each element in the group
    subplant_crosswalk["numeric_generator_id"] = subplant_crosswalk.groupby(
        ["plant_id_eia", "generator_id"], dropna=False
    ).ngroup()
    # when filling in missing unit_id_pudl, we don't want these numeric_generator_id to overlap existing unit_id
    # to ensure this, we will add 1000 to each of these numeric generator ids to ensure they are unique
    subplant_crosswalk["numeric_generator_id"] = (
        subplant_crosswalk["numeric_generator_id"] + 1000
    )
    # fill any missing unit_id_pudl with a number for each unique generator
    subplant_crosswalk["unit_id_pudl_filled"] = (
        subplant_crosswalk["unit_id_pudl_connected"]
        .fillna(subplant_crosswalk["subplant_id_connected"] + 100)
        .fillna(subplant_crosswalk["numeric_generator_id"])
    )
    # create a new unique subplant_id based on the connected subplant ids and the filled unit_id
    subplant_crosswalk["new_subplant"] = subplant_crosswalk.groupby(
        ["plant_id_eia", "subplant_id_connected", "unit_id_pudl_filled"],
        dropna=False,
    ).ngroup()

    return subplant_crosswalk


def connect_ids(df, id_to_update, connecting_id):
    """Corrects an id value if it is connected by an id value in another column.

    if multiple subplant_id are connected by a single unit_id_pudl, this groups these subplant_id together
    if multiple unit_id_pudl are connected by a single subplant_id, this groups these unit_id_pudl together

    Args:
        df: dataframe containing columns with id_to_update and connecting_id columns
        subplant_unit_pairs
    """

    # get a table with all unique subplant to unit pairs
    subplant_unit_pairs = df[
        ["plant_id_eia", "subplant_id", "unit_id_pudl"]
    ].drop_duplicates()

    # identify if any non-NA id_to_update are duplicated, indicated that it is associated with multiple connecting_id
    duplicates = subplant_unit_pairs[
        (subplant_unit_pairs.duplicated(subset=id_to_update, keep=False))
        & (~subplant_unit_pairs[id_to_update].isna())
    ].copy()

    # if there are any duplicate units, indicating an incorrect id_to_update, fix the id_to_update
    df[f"{connecting_id}_connected"] = df[connecting_id]
    if len(duplicates) > 0:
        # find the lowest number subplant id associated with each duplicated unit_id_pudl
        duplicates.loc[:, f"{connecting_id}_to_replace"] = (
            duplicates.groupby(["plant_id_eia", id_to_update])[connecting_id]
            .min()
            .iloc[0]
        )
        # merge this replacement subplant_id into the dataframe and use it to update the existing subplant id
        df = df.merge(
            duplicates,
            how="left",
            on=["plant_id_eia", id_to_update, connecting_id],
            validate="m:1",
        )
        df[f"{connecting_id}_connected"].update(df[f"{connecting_id}_to_replace"])
    return df


def add_operating_and_retirement_dates(df, start_year, end_year):
    """Adds columns listing a generator's planned operating date or retirement date to a dataframe."""
    pudl_db = f"sqlite:///{downloads_folder()}pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)
    # get values starting with the year prior to teh start year so that we can get proposed operating dates for the start year (which are reported in year -1)
    pudl_out_status = pudl.output.pudltabl.PudlTabl(
        pudl_engine,
        freq="MS",
        start_date=f"{start_year-1}-01-01",
        end_date=f"{end_year}-12-31",
    )
    generator_status = pudl_out_status.gens_eia860().loc[
        :,
        [
            "plant_id_eia",
            "generator_id",
            "report_date",
            "operational_status",
            "current_planned_operating_date",
            "retirement_date",
        ],
    ]
    # only keep values that have a planned operating date or retirement date
    generator_status = generator_status[
        (~generator_status["current_planned_operating_date"].isna())
        | (~generator_status["retirement_date"].isna())
    ]
    # drop any duplicate entries
    generator_status = generator_status.sort_values(
        by=["plant_id_eia", "generator_id", "report_date"]
    ).drop_duplicates(
        subset=[
            "plant_id_eia",
            "generator_id",
            "current_planned_operating_date",
            "retirement_date",
        ],
        keep="last",
    )
    # for any generators that have different retirement or planned dates reported in different years, keep the most recent value
    generator_status = generator_status.sort_values(
        by=["plant_id_eia", "generator_id", "report_date"]
    ).drop_duplicates(subset=["plant_id_eia", "generator_id"], keep="last")

    # merge the dates into the crosswalk
    df = df.merge(
        generator_status[
            [
                "plant_id_eia",
                "generator_id",
                "current_planned_operating_date",
                "retirement_date",
            ]
        ],
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )

    return df


def manual_crosswalk_updates(crosswalk):
    # load manual matches
    crosswalk_manual = pd.read_csv(
        manual_folder("epa_eia_crosswalk_manual.csv"),
        dtype=get_dtypes(),
    ).drop(columns=["notes"])
    crosswalk_manual = crosswalk_manual.rename(
        columns={
            "plant_id_epa": "CAMD_PLANT_ID",
            "unitid": "CAMD_UNIT_ID",
            "plant_id_eia": "EIA_PLANT_ID",
            "generator_id": "EIA_GENERATOR_ID",
        }
    )

    # The EPA's crosswalk document incorrectly maps plant_id_epa 55248 to plant_id_eia 55248
    # the correct plant_id_eia is 2847
    crosswalk.loc[crosswalk["CAMD_PLANT_ID"] == 55248, "EIA_PLANT_ID"] = 2847

    # move missing crosswalk matches to a different dataframe
    unmatched = crosswalk.copy()[crosswalk["EIA_GENERATOR_ID"].isna()]
    crosswalk = crosswalk[~crosswalk["EIA_GENERATOR_ID"].isna()]

    # append the manual data to the crosswalk
    crosswalk = pd.concat([crosswalk, crosswalk_manual], axis=0)

    # filter the list of unmatched generators to those that were not in our manual list
    unmatched = unmatched.merge(
        crosswalk_manual,
        how="outer",
        on=["CAMD_PLANT_ID", "CAMD_UNIT_ID"],
        indicator="source",
        suffixes=(None, "_manual"),
        validate="m:m",
    )
    unmatched = unmatched[unmatched["source"] == "left_only"].drop(
        columns=["EIA_PLANT_ID_manual", "EIA_GENERATOR_ID_manual", "source"]
    )

    # add these back to the crosswalk
    crosswalk = pd.concat([crosswalk, unmatched], axis=0)

    return crosswalk


def clean_eia923(
    year: int,
    small: bool,
    add_subplant_id: bool = True,
    calculate_nox_emissions: bool = True,
    calculate_so2_emissions: bool = True,
):
    """
    This is the coordinating function for cleaning and allocating generation and fuel data in EIA-923.
    """
    # Distribute net generation and heat input data reported by the three different EIA-923 tables

    pudl_out = load_data.initialize_pudl_out(year=year)

    # allocate net generation and heat input to each generator-fuel grouping
    gen_fuel_allocated = allocate_gen_fuel.allocate_gen_fuel_by_generator_energy_source(
        pudl_out, drop_interim_cols=True
    )

    # manually update energy source code when OTH
    gen_fuel_allocated = update_energy_source_codes(gen_fuel_allocated)

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
    ].round(
        1
    )

    validation.test_for_missing_energy_source_code(gen_fuel_allocated)
    validation.test_for_negative_values(gen_fuel_allocated)

    # create a table that identifies the primary fuel of each generator and plant
    primary_fuel_table = create_primary_fuel_table(gen_fuel_allocated, pudl_out)

    if small:
        gen_fuel_allocated = smallerize_test_data(df=gen_fuel_allocated, random_seed=42)

    # calculate co2 emissions for each generator-fuel based on allocated fuel consumption
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
            gen_fuel_allocated, pudl_out, year
        )
    if calculate_so2_emissions:
        gen_fuel_allocated = emissions.calculate_so2_from_fuel_consumption(
            gen_fuel_allocated, pudl_out, year
        )

    # adjust total emissions for biomass
    gen_fuel_allocated = emissions.adjust_emissions_for_biomass(gen_fuel_allocated)

    # adjust emissions for CHP
    gen_fuel_allocated = emissions.adjust_fuel_and_emissions_for_CHP(gen_fuel_allocated)

    gen_fuel_allocated = emissions.calculate_co2e_mass(
        gen_fuel_allocated, year, gwp_horizon=100, ar5_climate_carbon_feedback=True
    )

    validation.test_emissions_adjustments(gen_fuel_allocated)

    # aggregate the allocated data to the generator level
    gen_fuel_allocated = allocate_gen_fuel.agg_by_generator(
        gen_fuel_allocated,
        sum_cols=DATA_COLUMNS,
    )

    # remove any plants that we don't want in the data
    gen_fuel_allocated = remove_plants(
        gen_fuel_allocated,
        non_grid_connected=True,
        remove_states=["PR"],
        steam_only_plants=False,
        distribution_connected_plants=False,
    )

    # round all values to the nearest tenth of a unit
    gen_fuel_allocated.loc[:, DATA_COLUMNS] = gen_fuel_allocated.loc[
        :, DATA_COLUMNS
    ].round(1)

    # NOTE(milo): Subplant IDs are required for the hourly emissions pipeline, but not
    # needed for exporting standalone EIA-923 data. We allow the user to skip the merge
    # below with a flag.
    if add_subplant_id:
        subplant_crosswalk = pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
            dtype=get_dtypes(),
        )[["plant_id_eia", "generator_id", "subplant_id"]].drop_duplicates()
        gen_fuel_allocated = gen_fuel_allocated.merge(
            subplant_crosswalk,
            how="left",
            on=["plant_id_eia", "generator_id"],
            validate="m:1",
        )

    # add the cleaned prime mover code to the data
    gen_pm = pudl_out.gens_eia860()[
        ["plant_id_eia", "generator_id", "prime_mover_code"]
    ]
    gen_fuel_allocated = gen_fuel_allocated.merge(
        gen_pm, how="left", on=["plant_id_eia", "generator_id"], validate="m:1"
    )

    gen_fuel_allocated = apply_dtypes(gen_fuel_allocated)
    primary_fuel_table = apply_dtypes(primary_fuel_table)

    # run validation checks on EIA-923 data
    validation.test_for_negative_values(gen_fuel_allocated)

    return gen_fuel_allocated, primary_fuel_table


def update_energy_source_codes(df):
    """
    Manually update fuel source codes
    """
    # load the table of updated fuel types
    updated_esc = pd.read_csv(manual_folder("updated_oth_energy_source_codes.csv"))

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
        print(
            "WARNING: After cleaning energy source codes, some fuel consumption is still associated with an 'OTH' fuel type."
        )
        print("This will lead to incorrect emissions calculations.")
        print(
            f"Check the following plants: {list(plants_with_other_fuel.plant_id_eia.unique())}"
        )
        print("Assign a fuel type in `data_cleaning.update_energy_source_codes`")

    return df


def create_primary_fuel_table(gen_fuel_allocated, pudl_out):
    """
    Identifies the primary fuel for each generator and plant
    Gen primary fuel is identified based on the "energy source code 1" identified in EIA-860
    Plant primary fuel is based on the most-consumed fuel at a plant based on allocated heat input
    """

    primary_fuel_from_capacity = calculate_capacity_based_primary_fuel(pudl_out)

    # get a table of primary energy source codes
    gen_primary_fuel = gen_fuel_allocated[
        gen_fuel_allocated["energy_source_code_num"] == "energy_source_code_1"
    ].drop_duplicates(subset=["plant_id_eia", "generator_id"])[
        ["plant_id_eia", "generator_id", "energy_source_code"]
    ]

    # NOTE(milo): In some rare cases, a plant will have no fuel specified by
    # energy_source_code_1, and will have zero fuel consumption and net generation for
    # all fuel types. When that happens, we simply assign a plant to have the same fuel
    # type as the majority of its generators.
    primary_fuel_from_mode = (
        gen_primary_fuel.groupby("plant_id_eia", dropna=False)["energy_source_code"]
        .agg(lambda x: pd.Series.mode(x)[0])
        .to_frame()
        .reset_index()
        .rename(columns={"energy_source_code": "primary_fuel_from_mode"})
    )

    # create a blank dataframe with all of the plant ids to hold primary fuel data
    plant_primary_fuel = gen_fuel_allocated[["plant_id_eia"]].drop_duplicates()

    # calculate the total annual fuel consumption, generation, and capacity by fuel type
    #  for each plant
    plant_totals_by_fuel = (
        gen_fuel_allocated.groupby(["plant_id_eia", "energy_source_code"], dropna=False)
        .sum()[["fuel_consumed_mmbtu", "net_generation_mwh"]]
        .reset_index()
    )

    # we will calculate primary fuel based on the fuel with the most consumption,
    # generation, and capacity
    for source in ["fuel_consumed_mmbtu", "net_generation_mwh"]:

        # only keep values greater than zero so that these can be filled by other
        # methods if non-zero
        primary_fuel_calc = plant_totals_by_fuel[plant_totals_by_fuel[source] > 0]

        # identify the fuel type with the maximum value for each plant
        primary_fuel_calc = primary_fuel_calc[
            primary_fuel_calc.groupby("plant_id_eia", dropna=False)[source].transform(
                max
            )
            == primary_fuel_calc[source]
        ][["plant_id_eia", "energy_source_code"]]

        # remove duplicate values (if two fuels are both the maximum)
        primary_fuel_calc = primary_fuel_calc.drop_duplicates(
            subset="plant_id_eia", keep=False
        )
        primary_fuel_calc = primary_fuel_calc.rename(
            columns={"energy_source_code": f"primary_fuel_from_{source}"}
        )

        # merge the primary fuel into the main table
        plant_primary_fuel = plant_primary_fuel.merge(
            primary_fuel_calc, how="left", on="plant_id_eia", validate="1:1"
        )

    # merge the primary fuel into the main table
    plant_primary_fuel = plant_primary_fuel.merge(
        primary_fuel_from_capacity, how="left", on="plant_id_eia", validate="1:1"
    )

    plant_primary_fuel = plant_primary_fuel.merge(
        primary_fuel_from_mode, how="left", on="plant_id_eia", validate="1:1"
    )

    # Use the fuel consumption-based primary fuel first, then fill using capacity-based
    # primary fuel, then generation based. Finally, to break all ties, use the energy
    # source code that appears most often for generators of a plant (mode).
    plant_primary_fuel["plant_primary_fuel"] = plant_primary_fuel[
        "primary_fuel_from_fuel_consumed_mmbtu"
    ]
    plant_primary_fuel["plant_primary_fuel"] = plant_primary_fuel[
        "plant_primary_fuel"
    ].fillna(plant_primary_fuel["primary_fuel_from_capacity_mw"])
    plant_primary_fuel["plant_primary_fuel"] = plant_primary_fuel[
        "plant_primary_fuel"
    ].fillna(plant_primary_fuel["primary_fuel_from_net_generation_mwh"])
    plant_primary_fuel["plant_primary_fuel"] = plant_primary_fuel[
        "plant_primary_fuel"
    ].fillna(plant_primary_fuel["primary_fuel_from_mode"])

    if len(plant_primary_fuel[plant_primary_fuel["plant_primary_fuel"].isna()]) > 0:
        plants_with_no_primary_fuel = plant_primary_fuel[
            plant_primary_fuel["plant_primary_fuel"].isna()
        ]
        print(
            f"Check the following plants: {list(plants_with_no_primary_fuel.plant_id_eia.unique())}"
        )
        raise UserWarning(
            "Plant primary fuel table contains missing primary fuels.\
            Update method of `create_primary_fuel_table()` to fix"
        )

    # merge the plant primary fuel into the gen primary fuel
    primary_fuel_table = gen_primary_fuel.merge(
        plant_primary_fuel[["plant_id_eia", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
        validate="many_to_one",
    )

    return primary_fuel_table


def calculate_capacity_based_primary_fuel(pudl_out):
    # create a table of primary fuel by nameplate capacity
    gen_capacity = pudl_out.gens_eia860().loc[
        :, ["plant_id_eia", "generator_id", "capacity_mw", "energy_source_code_1"]
    ]

    gen_capacity = (
        gen_capacity.groupby(["plant_id_eia", "energy_source_code_1"], dropna=False)
        .sum()["capacity_mw"]
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
        gen_capacity.groupby("plant_id_eia", dropna=False)["capacity_mw"].transform(max)
        == gen_capacity["capacity_mw"]
    ][["plant_id_eia", "energy_source_code_1"]].rename(
        columns={"energy_source_code_1": "primary_fuel_from_capacity_mw"}
    )

    # drop any duplicate entries (if two fuel types have the same nameplate capacity)
    gen_capacity = gen_capacity[
        ~(gen_capacity.duplicated(subset="plant_id_eia", keep=False))
    ]

    return gen_capacity


def remove_plants(
    df,
    non_grid_connected=False,
    remove_states=[],
    steam_only_plants=False,
    distribution_connected_plants=False,
):
    """
    Coordinating function to remove specific plants based on specified options
    Each function should identify how many plants are being removed
    Args:
        df: dataframe containing plant_id_eia column
        non_grid_connected: if True, remove all plants that are not grid connected
        remove_states: list of two-letter state codes for which plants should be removed if located within
        steam_only_plants: if True, remove plants that only generate heat and no electricity (not yet implemented)
        distribution_connected_plants: if True, remove plants that are connected to the distribution grid (not yet implemented)
    """
    if non_grid_connected:
        df = remove_non_grid_connected_plants(df)
    if len(remove_states) > 0:
        plant_states = load_data.load_pudl_table("plants_entity_eia").loc[
            :, ["plant_id_eia", "state"]
        ]
        plants_in_states_to_remove = list(
            plant_states[
                plant_states["state"].isin(remove_states)
            ].plant_id_eia.unique()
        )
        print(
            f"    Removing {len(plants_in_states_to_remove)} plants located in the following states: {remove_states}"
        )
        df = df[~df["plant_id_eia"].isin(plants_in_states_to_remove)]
    if steam_only_plants:
        pass
    if distribution_connected_plants:
        pass

    return df


def remove_non_grid_connected_plants(df):
    """
    Removes any records from a dataframe associated with plants that are not connected to the electricity grid
    Inputs:
        df: any pandas dataframe containing the column 'plant_id_eia'
    Returns:
        df: pandas dataframe with non-grid connected plants removed
    """

    # get the list of plant_id_eia from the static table
    ngc_plants = list(
        pd.read_csv(
            manual_folder("plants_not_connected_to_grid.csv"),
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
    print(f"    Removing {num_plants} plants that are not grid-connected")

    df = df[~df["plant_id_eia"].isin(ngc_plants)]

    # according to the egrid documentation, any plants that have an id of 88XXXX are not grid connected
    # only keep plants that dont have an id of 88XXXX
    df = df[(df["plant_id_eia"] < 880000) | (df["plant_id_eia"] >= 890000)]

    return df


def clean_cems(year: int, small: bool, primary_fuel_table):
    """
    Coordinating function for all of the cems data cleaning
    """
    # load the CEMS data
    cems = load_data.load_cems_data(year)

    if small:
        cems = smallerize_test_data(df=cems, random_seed=42)

    # remove non-grid connected plants
    cems = remove_plants(
        cems,
        non_grid_connected=True,
        remove_states=["PR"],
        steam_only_plants=False,
        distribution_connected_plants=False,
    )

    # manually remove steam-only units
    cems = manually_remove_steam_units(cems)

    # add a report date
    cems = load_data.add_report_date(cems)

    # remove data for any unit-months where there are incomplete data reported
    # this is generally when there is a single observation reported for an entire month
    cems = remove_incomplete_unit_months(cems)

    # TODO: identify and remove any hourly values that appear to be outliers

    # add a fuel type to each observation
    cems = assign_fuel_type_to_cems(cems, year, primary_fuel_table)

    # fill in missing hourly emissions data using the fuel type and heat input
    validation.test_for_missing_energy_source_code(cems)
    cems = emissions.fill_cems_missing_co2(cems, year)

    # TODO: Add functions for filling missing NOx and SOx

    # calculate ch4 and n2o emissions
    cems = emissions.calculate_ghg_emissions_from_fuel_consumption(
        df=cems, year=year, include_co2=False, include_ch4=True, include_n2o=True
    )

    # remove any observations from cems where zero operation is reported for an entire month
    # although this data could be considered to be accurately reported, let's remove it so that we can double check against the eia data
    # TODO: check if any of these observations are from geothermal generators
    cems = remove_cems_with_zero_monthly_data(cems)

    validation.test_for_negative_values(cems)

    # add subplant id
    subplant_crosswalk = (
        pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
            dtype=get_dtypes(),
        )[["plant_id_eia", "unitid", "subplant_id"]]
        .drop_duplicates()
        .dropna(subset="unitid")
    )
    cems = cems.merge(
        subplant_crosswalk, how="left", on=["plant_id_eia", "unitid"], validate="m:1"
    )
    validation.test_for_missing_subplant_id(cems)

    cems = apply_dtypes(cems)

    return cems


def smallerize_test_data(df, random_seed=None):
    print("    Randomly selecting 5% of plants for faster test run.")
    # Select 5% of plants
    selected_plants = df.plant_id_eia.unique()
    if random_seed is not None:
        np.random.seed(random_seed)
    selected_plants = np.random.choice(
        selected_plants, size=int(len(selected_plants) * 0.05), replace=False
    )
    # Filter for selected plants
    df = df[df.plant_id_eia.isin(selected_plants)]

    return df


def manually_remove_steam_units(df):
    """
    Removes any records from CEMS that we've identified as being steam only plants that need to be removed
    """

    # get the list of plant_id_eia from the static table
    units_to_remove = pd.read_csv(
        manual_folder("steam_units_to_remove.csv"),
        dtype=get_dtypes(),
    )[["plant_id_eia", "unitid"]]

    print(
        f"    Removing {len(units_to_remove)} units that only produce steam and do not report to EIA"
    )

    df = df.merge(
        units_to_remove,
        how="outer",
        on=["plant_id_eia", "unitid"],
        indicator="source",
        validate="m:1",
    )
    df = df[df["source"] == "left_only"].drop(columns=["source"])

    return df


def remove_incomplete_unit_months(cems):

    # get a count of how many hours are reported in each month for each unit
    unit_hours_in_month = (
        cems[["plant_id_eia", "report_date", "unitid", "datetime_utc"]]
        .groupby(["plant_id_eia", "report_date", "unitid"], dropna=False)
        .count()
        .reset_index()
    )

    # identify months where there is not complete data
    # The fewest number of hours in a month is 28*24 = 672
    unit_months_to_remove = unit_hours_in_month[
        unit_hours_in_month["datetime_utc"] < 600
    ].drop(columns="datetime_utc")

    print(
        f"    Removing {len(unit_months_to_remove)} unit-months with incomplete hourly data"
    )

    cems = cems.merge(
        unit_months_to_remove,
        how="outer",
        on=["plant_id_eia", "report_date", "unitid"],
        validate="m:1",
        indicator="to_remove",
    )

    cems = cems[cems["to_remove"] != "both"].drop(columns="to_remove")

    return cems


def assign_fuel_type_to_cems(cems, year, primary_fuel_table):
    "Assigns a fuel type to each observation in CEMS"

    fuel_types = get_epa_unit_fuel_types(year)

    # merge in the reported fuel type
    cems = cems.merge(
        fuel_types, how="left", on=["plant_id_epa", "unitid"], validate="m:1"
    )

    # fill missing fuel codes for plants that only have a single fuel type
    single_fuel_plants = (
        fuel_types.drop_duplicates(subset=["plant_id_epa", "energy_source_code"])
        .drop_duplicates(subset=["plant_id_epa"], keep=False)
        .drop(columns="unitid")
    )
    cems = cems.merge(
        single_fuel_plants,
        how="left",
        on=["plant_id_epa"],
        suffixes=(None, "_plant"),
        validate="m:1",
    )
    cems["energy_source_code"] = cems["energy_source_code"].fillna(
        cems["energy_source_code_plant"]
    )
    cems = cems.drop(columns="energy_source_code_plant")

    # TODO: fill fuel codes for plants that only have a single fossil type identified in EIA
    cems = fill_missing_fuel_for_single_fuel_plant_months(cems, year)

    # fill any remaining missing fuel codes with the plant primary fuel identified from EIA-923
    cems = cems.merge(
        primary_fuel_table[["plant_id_eia", "plant_primary_fuel"]].drop_duplicates(),
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )
    cems["energy_source_code"] = cems["energy_source_code"].fillna(
        cems["plant_primary_fuel"]
    )
    cems = cems.drop(columns=["plant_primary_fuel"])

    # update
    cems = update_energy_source_codes(cems)

    return cems


def get_epa_unit_fuel_types(year):
    """
    Loads the energy source code assigned to each CAMD unit in the EPA-EIA crosswalk.

    If the EIA fuel type is missing, uses the CAMD fuel type to fill.
    """
    # get a unique list of plant unit fuels
    fuel_types = load_data.load_epa_eia_crosswalk(year)[
        ["plant_id_epa", "unitid", "energy_source_code_eia"]
    ].drop_duplicates()

    fuel_types = fuel_types.dropna(subset="energy_source_code_eia")

    # remove any entries where there are multiple fuel types listed
    fuel_types = fuel_types[
        ~fuel_types[["plant_id_epa", "unitid"]].duplicated(keep=False)
    ]

    # rename the column
    fuel_types = fuel_types.rename(
        columns={"energy_source_code_eia": "energy_source_code"}
    )

    return fuel_types


def fill_missing_fuel_for_single_fuel_plant_months(df, year):
    """
    Identifies all plant-months where a single fuel was burned based on the EIA-923 generation fuel table
    Uses this to fill in the energy source code if a match was not made based on the PSDC
    """

    # identify plant-months for which there is a single fossil fuel type reported
    gf = load_data.load_pudl_table("generation_fuel_eia923", year=year)[
        ["plant_id_eia", "report_date", "energy_source_code", "fuel_consumed_mmbtu"]
    ]

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
    gf = gf.rename(columns={"energy_source_code": "energy_source_code_single"}).drop(
        columns=["num_fuels", "fuel_consumed_mmbtu"]
    )
    gf["report_date"] = pd.to_datetime(gf["report_date"])

    # merge this data into the df
    df = df.merge(gf, how="left", on=["plant_id_eia", "report_date"], validate="m:1")

    # fill missing fuel types with this data
    df["energy_source_code"] = df["energy_source_code"].fillna(
        df["energy_source_code_single"]
    )

    # remove the intermediate column
    df = df.drop(columns=["energy_source_code_single"])

    return df


def remove_cems_with_zero_monthly_data(cems):
    """
    Identifies months where zero generation or heat inputare reported.
    from each unit and removes associated hours from CEMS so that these can be filled using the eia923 data
    Inputs:
        cems: pandas dataframe of hourly cems data containing columns "plant_id_eia", "unitid" and "report_date"
    Returns:
        cems df with hourly observations for months when no emissions reported removed
    """
    # calculate the totals reported in each month
    cems_with_zero_monthly_emissions = cems.groupby(
        ["plant_id_eia", "unitid", "report_date"], dropna=False
    ).sum()[["gross_generation_mwh", "fuel_consumed_mmbtu"]]
    # identify unit-months where zero emissions reported
    cems_with_zero_monthly_emissions = cems_with_zero_monthly_emissions[
        cems_with_zero_monthly_emissions.sum(axis=1) == 0
    ]
    # add a flag to these observations
    cems_with_zero_monthly_emissions["missing_data_flag"] = "remove"

    # merge the missing data flag into the cems data
    cems = cems.merge(
        cems_with_zero_monthly_emissions.reset_index()[
            ["plant_id_eia", "unitid", "report_date", "missing_data_flag"]
        ],
        how="left",
        on=["plant_id_eia", "unitid", "report_date"],
        validate="m:1",
    )
    # remove any observations with the missing data flag
    print(
        f"    Removing {len(cems[cems['missing_data_flag'] == 'remove'])} observations from cems for unit-months where no data reported"
    )
    cems = cems[cems["missing_data_flag"] != "remove"]
    # drop the missing data flag column
    cems = cems.drop(columns="missing_data_flag")

    return cems


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
        )
        .sum()[["fuel_consumed_mmbtu", "fuel_consumed_for_electricity_mmbtu"]]
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
        eia923_allocated.groupby(["plant_id_eia", "report_date"], dropna=False)
        .sum()[["fuel_consumed_mmbtu", "fuel_consumed_for_electricity_mmbtu"]]
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

    # remove intermediate columns
    cems = cems.drop(columns=["subplant_fuel_ratio", "plant_fuel_ratio"])

    # add adjusted emissions columns
    cems = emissions.adjust_fuel_and_emissions_for_CHP(cems)

    return cems


def identify_hourly_data_source(eia923_allocated, cems, year):
    """Identifies whether there is hourly CEMS data available for each subplant-month.
    Possible categories:
        1. `cems`: For subplant-months for which we have hourly CEMS data for all CEMS units that make up that subplant,
            we will use the hourly values reported in CEMS. (Add a validation check for the net generation and fuel consumption totals)
        2. `partial_cems`: For subplant-months for which we have hourly CEMS data
            for only some of the CEMS units that make up a subplant, we will use the reported
            EIA-923 values to scale the partial hourly CEMS data from the other units to match the total value for the entire subplant. This will also calculate a partial subplant scaling factor for each data column (e.g. net generation, fuel consumption) by comparing the total monthly CEMS data to the monthly EIA-923 data.
        3. `eia`: for subplant-months for which no hourly data is reported in CEMS,
            we will attempt to use EIA-930 data to assign an hourly profile to the monthly EIA-923 data
    Inputs:
        eia923_allocated:
        cems:
        year:
    Returns:
        eia923_allocated with new column `hourly_data_source`
    """

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

    # remove any generator-months for which there is no data reported in either data source
    all_data = all_data[
        ~(
            (all_data["reported_eia923"] == 0)
            & (all_data["hourly_data_source"] == "eia")
        )
    ]

    all_data = all_data.drop(columns=["reported_eia923"])

    # identify the partial cems plants
    all_data = identify_partial_cems_plants(all_data)

    return all_data


def identify_partial_cems_subplants(year, cems, eia923_allocated):
    """Identifies subplants for which data for only some units is reported in CEMS.

    If the number of unique unitid reported in cems for a subplant is less than the
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
            ["plant_id_eia", "subplant_id", "unitid", "report_date"], dropna=False
        )
        .sum()[["fuel_consumed_mmbtu"]]
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
        )
        .sum(min_count=1)["fuel_consumed_mmbtu"]
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
    # we identify this as subplant months where the number of CEMS reported units is less
    # than the total number of units in that subplant, AND the total fuel reported in CEMS
    # is less than 95% of the fuel reported in EIA
    cems_status = cems_units_reported.assign(
        hourly_data_source=lambda x: np.where(
            (
                (x.reported_units_in_subplant < x.units_in_subplant)
                & (x.fuel_consumed_mmbtu_cems < (x.fuel_consumed_mmbtu_eia * 0.95))
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
    # load the subplant crosswalk and identify unique unitids in each subplant
    units_in_subplant = (
        pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
            dtype=get_dtypes(),
            parse_dates=["current_planned_operating_date", "retirement_date"],
        )[["plant_id_eia", "unitid", "subplant_id", "retirement_date"]]
        .drop_duplicates()
        .dropna(subset="unitid")
    )

    # remove units that retired before the current year
    units_in_subplant = units_in_subplant[
        ~(units_in_subplant["retirement_date"].dt.year < year)
    ]

    # get a count of the number of CEMS units in each subplant
    units_in_subplant = (
        units_in_subplant.groupby(["plant_id_eia", "subplant_id"], dropna=False)
        .count()["unitid"]
        .reset_index()
        .rename(columns={"unitid": "units_in_subplant"})
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
        .agg({"unitid": "count", "fuel_consumed_mmbtu": "sum"})
        .reset_index()
        .rename(columns={"unitid": "reported_units_in_subplant"})
    )

    return reported_units_in_subplant


def identify_partial_cems_plants(all_data):
    """Identifies subplants for which CEMS data is reported for at least one other subplant at the plant.

    Args:
        all_data: dataframe identifying the hourly data source for each subplant-month
    Returns:
        all_data with updated hourly_data_source column indicating partial cems plants"""

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

    # drop intermediate columns
    partial_plant = partial_plant.drop(columns=["eia_data", "cems_data"])

    # merge this data into all_data
    all_data = all_data.merge(
        partial_plant,
        how="outer",
        on=["plant_id_eia", "report_date"],
        indicator="partial_plant",
        validate="m:1",
    )

    # for plants where the data source is eia only, and it exists in partial plant, change the source to partial_cems_plant
    all_data.loc[
        (all_data["hourly_data_source"] == "eia")
        & (all_data["partial_plant"] == "both"),
        "hourly_data_source",
    ] = "partial_cems_plant"

    # remove the intermediate indicator column
    all_data = all_data.drop(columns=["partial_plant", "eia_data", "cems_data"])

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

    filtered_cems = filtered_cems[filtered_cems["source"] == "left_only"].drop(
        columns=["source"]
    )

    return filtered_cems


def aggregate_plant_data_to_ba_fuel(combined_plant_data, plant_attributes_table):

    # create a table that has data for the sythetic plant attributes
    shaped_plant_attributes = (
        plant_attributes_table[["shaped_plant_id", "ba_code", "fuel_category"]]
        .drop_duplicates()
        .dropna(subset="shaped_plant_id")
        .rename(columns={"shaped_plant_id": "plant_id_eia"})
    )

    combined_plant_attributes = pd.concat(
        [
            plant_attributes_table[["plant_id_eia", "ba_code", "fuel_category"]],
            shaped_plant_attributes,
        ],
        axis=0,
    )

    ba_fuel_data = combined_plant_data.merge(
        combined_plant_attributes, how="left", on=["plant_id_eia"], validate="m:1"
    )
    # check that there is no missing ba or fuel codes
    if (
        len(
            ba_fuel_data[
                (ba_fuel_data["ba_code"].isna())
                | (ba_fuel_data["fuel_category"].isna())
            ]
        )
        > 0
    ):
        raise UserWarning(
            "The plant attributes table is missing ba code or fuel_category data for some plants. This will result in incomplete power sector results."
        )
    ba_fuel_data = (
        ba_fuel_data.groupby(
            ["ba_code", "fuel_category", "datetime_utc", "report_date"], dropna=False
        )[DATA_COLUMNS]
        .sum()
        .reset_index()
    )
    return ba_fuel_data


def combine_plant_data(
    cems,
    partial_cems_subplant,
    partial_cems_plant,
    eia_data,
    resolution,
    validate=False,
):
    """
    Combines final hourly subplant data from each source into a single dataframe.
    Inputs:
        Pandas dataframes of shaped or original hourly data
        resolution: string, either 'monthly' or 'hourly'

    Note: returns dask dataframe (not used before this point in pipeline) because of data size

    """

    if resolution == "hourly":
        KEY_COLUMNS = [
            "plant_id_eia",
            "datetime_utc",
            "report_date",
        ]
    elif resolution == "monthly":
        KEY_COLUMNS = [
            "plant_id_eia",
            "report_date",
        ]
    else:
        raise UserWarning(
            "arg 'resolution' for `combine_plant_data` must be either 'monthly' or 'hourly'"
        )

    ALL_COLUMNS = KEY_COLUMNS + DATA_COLUMNS

    if validate:
        validation.ensure_non_overlapping_data_from_all_sources(
            cems, partial_cems_subplant, partial_cems_plant, eia_data
        )

    # group data by plant-hour and filter columns
    cems = (
        cems.groupby(
            KEY_COLUMNS,
            dropna=False,
        )
        .sum()
        .reset_index()[[col for col in cems.columns if col in ALL_COLUMNS]]
    )
    # don't group if there is no data in the dataframe
    if len(partial_cems_subplant) > 0:
        partial_cems_subplant = (
            partial_cems_subplant.groupby(
                KEY_COLUMNS,
                dropna=False,
            )
            .sum()
            .reset_index()[
                [col for col in partial_cems_subplant.columns if col in ALL_COLUMNS]
            ]
        )
    if len(partial_cems_plant) > 0:
        partial_cems_plant = (
            partial_cems_plant.groupby(
                KEY_COLUMNS,
                dropna=False,
            )
            .sum()
            .reset_index()[
                [col for col in partial_cems_plant.columns if col in ALL_COLUMNS]
            ]
        )
    eia_data = (
        eia_data.groupby(
            KEY_COLUMNS,
            dropna=False,
        )
        .sum()
        .reset_index()[[col for col in eia_data.columns if col in ALL_COLUMNS]]
    )

    # concat together
    combined_plant_data = pd.concat(
        [cems, partial_cems_subplant, partial_cems_plant, eia_data],
        axis=0,
        ignore_index=True,
        copy=False,
    )

    # groupby plant
    combined_plant_data = (
        combined_plant_data.groupby(KEY_COLUMNS, dropna=False).sum().reset_index()
    )

    combined_plant_data[DATA_COLUMNS] = combined_plant_data[DATA_COLUMNS].round(2)

    # re-order the columns
    combined_plant_data = combined_plant_data[ALL_COLUMNS]

    return combined_plant_data


def create_plant_attributes_table(cems, eia923_allocated, year, primary_fuel_table):

    # create a table with the unique plantids from both dataframes
    eia_plants = eia923_allocated[
        ["plant_id_eia", "plant_primary_fuel"]
    ].drop_duplicates()
    cems_plants = cems[["plant_id_eia"]].drop_duplicates()

    # merge primary fuel into cems
    cems_plants = cems_plants.merge(
        primary_fuel_table.drop_duplicates(subset="plant_id_eia")[
            ["plant_id_eia", "plant_primary_fuel"]
        ],
        how="left",
        on="plant_id_eia",
        validate="1:1",
    )

    # identify any CEMS-only plants that are missing a primary fuel assignment
    plants_missing_primary_fuel = list(
        cems_plants.loc[cems_plants["plant_primary_fuel"].isna(), "plant_id_eia"]
    )

    # calculate primary fuel for each of these missing plants
    cems_primary_fuel = cems.loc[
        cems["plant_id_eia"].isin(plants_missing_primary_fuel),
        ["plant_id_eia", "energy_source_code", "fuel_consumed_mmbtu"],
    ]

    # calculate the total fuel consumption by fuel type and keep the fuel code with the largest fuel consumption
    cems_primary_fuel = (
        cems_primary_fuel.groupby(["plant_id_eia", "energy_source_code"], dropna=False)
        .sum()
        .reset_index()
        .sort_values(by="fuel_consumed_mmbtu", ascending=False)
        .drop_duplicates(subset="plant_id_eia", keep="first")
        .drop(columns="fuel_consumed_mmbtu")
    )

    # merge the cems primary fuel back in and use it to fill any missing fuel codes
    cems_plants = cems_plants.merge(
        cems_primary_fuel,
        how="left",
        on="plant_id_eia",
        validate="1:1",
    )
    cems_plants["plant_primary_fuel"] = cems_plants["plant_primary_fuel"].fillna(
        cems_plants["energy_source_code"]
    )
    cems_plants = cems_plants.drop(columns="energy_source_code")

    # concat the two lists together
    plant_attributes = eia_plants.merge(
        cems_plants,
        how="outer",
        on="plant_id_eia",
        indicator="data_availability",
        suffixes=(None, "_cems"),
        validate="1:1",
    )
    plant_attributes["plant_primary_fuel"] = plant_attributes[
        "plant_primary_fuel"
    ].fillna(plant_attributes["plant_primary_fuel_cems"])
    plant_attributes = plant_attributes.drop(columns=["plant_primary_fuel_cems"])
    plant_attributes["data_availability"] = plant_attributes[
        "data_availability"
    ].replace(
        {"left_only": "eia_only", "right_only": "cems_only", "both": "cems_and_eia"}
    )

    # assign a BA code and state code to each plant
    plant_attributes = assign_ba_code_to_plant(plant_attributes, year)

    # add a flag about whether the plant is distribution connected
    plant_attributes = identify_distribution_connected_plants(
        plant_attributes, year, voltage_threshold_kv=60
    )

    # assign a fuel category to each plant based on what is most likely to match with the category used in EIA-930
    plant_attributes = assign_fuel_category_to_ESC(
        df=plant_attributes,
        esc_column="plant_primary_fuel",
    )

    # add tz info
    plant_attributes = add_plant_local_timezone(plant_attributes, year)

    plant_attributes = apply_dtypes(plant_attributes)

    return plant_attributes


def assign_ba_code_to_plant(df, year):
    """
    Assigns a balancing authority code and state to each plant based on the plant id
    Inputs:
        df: a pandas dataframe containing a 'plant_id_eia' column
        year: four digit year number for the data
    Returns:
        df with a new column for 'ba_code' and 'state'
    """

    plant_ba = create_plant_ba_table(year)[
        ["plant_id_eia", "ba_code", "ba_code_physical", "state"]
    ]

    # merge the ba code into the dataframe
    df = df.merge(plant_ba, how="left", on="plant_id_eia", validate="m:1")

    if len(df[df["ba_code"].isna()]) > 0:
        print("    WARNING: the following plants are missing ba_code:")
        print(df[df["ba_code"].isna()])

    # replace missing ba codes with NA
    df["ba_code"] = df["ba_code"].fillna("NA")
    df["ba_code_physical"] = df["ba_code_physical"].fillna("NA")

    return df


def create_plant_ba_table(year):
    """
    Creates a table assigning a ba_code and physical ba code to each plant id.
    """
    pudl_out = load_data.initialize_pudl_out(year=year)

    plant_ba = pudl_out.plants_eia860().loc[
        :,
        [
            "plant_id_eia",
            "balancing_authority_code_eia",
            "balancing_authority_name_eia",
            "utility_name_eia",
            "transmission_distribution_owner_name",
            "state",
        ],
    ]

    # convert the dtype of the balancing authority code column from string to object
    # this will allow for missing values to be filled
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].astype(object)
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(value=np.NaN)

    # specify a ba code for certain utilities
    utility_as_ba_code = pd.read_csv(manual_folder("utility_name_ba_code_map.csv"))
    utility_as_ba_code = dict(
        zip(
            utility_as_ba_code["name"],
            utility_as_ba_code["ba_code"],
        )
    )

    # fill missing BA codes first based on the BA name, then utility name, then on the transmisison owner name
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["balancing_authority_name_eia"].map(utility_as_ba_code))
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["utility_name_eia"].map(utility_as_ba_code))
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["transmission_distribution_owner_name"].map(utility_as_ba_code))

    # rename the ba column
    plant_ba = plant_ba.rename(columns={"balancing_authority_code_eia": "ba_code"})

    plant_ba["ba_code"] = plant_ba["ba_code"].replace("None", np.NaN)

    # get a list of all of the BAs that retired prior to the current year
    retired_bas = load_data.load_ba_reference()[["ba_code", "retirement_date"]]
    retired_bas = list(
        retired_bas.loc[
            retired_bas["retirement_date"].dt.year < year, "ba_code"
        ].unique()
    )
    # if there are any plants that have been assigned to a retired BA, set its BA code as missing
    plant_ba.loc[plant_ba["ba_code"].isin(retired_bas), "ba_code"] = np.NaN

    # for plants without a BA code assign the miscellaneous BA code based on the state
    plant_ba["ba_code"] = plant_ba["ba_code"].fillna(plant_ba["state"] + "MS")

    # add a physical ba code based on the owner of the transmission system
    plant_ba["ba_code_physical"] = plant_ba["ba_code"]
    plant_ba["ba_code_physical"].update(
        plant_ba["transmission_distribution_owner_name"].map(utility_as_ba_code)
    )

    # update based on mapping table when ambiguous
    physical_ba = pd.read_csv(manual_folder("physical_ba.csv"), dtype=get_dtypes())
    plant_ba = plant_ba.merge(
        physical_ba,
        how="left",
        on=["ba_code", "transmission_distribution_owner_name"],
        suffixes=("", "_map"),
        validate="m:1",
    )
    plant_ba["ba_code_physical"].update(plant_ba["ba_code_physical_map"])

    return plant_ba


def identify_distribution_connected_plants(df, year, voltage_threshold_kv=60):
    """
    Identifies which plant_id_eia are "distribution grid connected" based on a voltage threshold.

    The distribution grid is generally considered to operate at 60-69kV and below.
    Thus, any plants that have a grid voltage under this threshold will be flagged as distribution connected.
    Args:
        df: pandas dataframe with a column for plant_id_eia
        voltage_threshold_kv: the voltage (kV) under which a plant will be considered to be a distribution asset
    Returns:
        df: with additional binary column `distribution_flag`
    """
    # load the EIA-860 data
    pudl_out = load_data.initialize_pudl_out(year=year)

    plant_voltage = pudl_out.plants_eia860().loc[
        :, ["plant_id_eia", "grid_voltage_1_kv"]
    ]

    plant_voltage = plant_voltage.assign(
        distribution_flag=lambda x: np.where(
            x.grid_voltage_1_kv <= voltage_threshold_kv, True, False
        )
    )

    df = df.merge(
        plant_voltage[["plant_id_eia", "distribution_flag"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )

    return df


def assign_fuel_category_to_ESC(
    df,
    fuel_category_names=["fuel_category", "fuel_category_eia930"],
    esc_column="energy_source_code",
):
    """
    Assigns a fuel category to each energy source code in a dataframe.
    Args:
        df: pandas dataframe with column name that matches fuel_category_name and contains energy source codes
        fuel_category_name: list of the columns in energy_source_groups.csv that contains the desired category mapping
        esc_column: name of the column in df that contains the energy source codes to assign a category to
    Returns:
        df with additional column for fuel category
    """
    # load the fuel category table
    energy_source_groups = pd.read_csv(
        manual_folder("energy_source_groups.csv"), dtype=get_dtypes()
    )[["energy_source_code"] + fuel_category_names].rename(
        columns={"energy_source_code": esc_column}
    )
    # assign a fuel category to the monthly eia data
    df = df.merge(
        energy_source_groups[[esc_column] + fuel_category_names],
        how="left",
        on=esc_column,
        validate="m:1",
    )

    return df


def add_plant_local_timezone(df, year):
    plant_tz = load_data.load_pudl_table("plants_entity_eia")[
        ["plant_id_eia", "timezone"]
    ]
    df = df.merge(plant_tz, how="left", on=["plant_id_eia"], validate="m:1")

    return df


def aggregate_cems_to_subplant(cems):

    GROUPBY_COLUMNS = ["plant_id_eia", "subplant_id", "datetime_utc", "report_date"]

    cems_columns_to_aggregate = [
        "gross_generation_mwh",
        "steam_load_1000_lb",
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
        cems.groupby(GROUPBY_COLUMNS, dropna=False)
        .sum()[cems_columns_to_aggregate]
        .reset_index()
    )

    cems = apply_dtypes(cems)

    return cems
