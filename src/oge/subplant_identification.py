import pandas as pd

import pudl.analysis.epacamd_eia as epacamd_eia
from pudl.etl.glue_assets import make_subplant_ids

import oge.load_data as load_data
import oge.validation as validation
from oge.constants import latest_validated_year
from oge.logging_util import get_logger

logger = get_logger(__name__)


def generate_subplant_ids() -> pd.DataFrame:
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

    # load all unique CEMS IDs from 2001 to present
    logger.info("loading CEMS ids")
    cems_ids = load_data.load_cems_ids()

    # load the crosswalk and filter it by the data that actually exists in cems
    crosswalk = load_data.load_epa_eia_crosswalk(latest_validated_year)

    # filter the crosswalk to drop any units that don't exist in CEMS
    filtered_crosswalk = epacamd_eia.filter_crosswalk(crosswalk, cems_ids)

    # update the subplant_crosswalk to ensure completeness
    # prepare the subplant crosswalk by adding a complete list of generators and adding
    # the unit_id_pudl column
    complete_gens = load_data.load_complete_eia_generators_for_subplants()
    filtered_crosswalk = filtered_crosswalk.merge(
        complete_gens,
        how="outer",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )
    # also add a complete list of cems emissions_unit_id_epa
    filtered_crosswalk = filtered_crosswalk.merge(
        cems_ids[["plant_id_eia", "emissions_unit_id_epa"]].drop_duplicates(),
        how="outer",
        on=["plant_id_eia", "emissions_unit_id_epa"],
        validate="m:1",
    )

    # use graph analysis to identify subplants
    crosswalk_with_subplant_ids = make_subplant_ids(filtered_crosswalk)

    # change the eia plant id to int
    crosswalk_with_subplant_ids["plant_id_eia"] = crosswalk_with_subplant_ids[
        "plant_id_eia"
    ].astype(int)

    # change the order of the columns
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids[
        [
            "plant_id_epa",
            "emissions_unit_id_epa",
            "plant_id_eia",
            "generator_id",
            "unit_id_pudl",
            "unit_id_eia",
            "unit_id_eia_numeric",
            "subplant_id",
            "operational_status_code",
            "earliest_report_date",
            "generator_operating_date",
            "generator_retirement_date",
            "original_planned_generator_operating_date",
            "current_planned_generator_operating_date",
            "prime_mover_code",
        ]
    ]

    # In order to help maintain a static subplant_id order, we want to sort the generators
    # by their commercial operating date from oldest to newest (this means that the oldest
    # plant will get the first subplant_id). Sometimes the generator operating date is
    # missing, if so fill with the retirement date (this covers generators that retired
    # a long time ago), then with the original planned operating date, which covers
    # proposed generators.
    # NOTE: because of shifting schedules, it is possible that a currently proposed
    # generator will come online in a different order than originally proposed, so there
    # is a small chance that the subplant_id will change if this is run with updated
    # EIA data.
    crosswalk_with_subplant_ids["sort_date"] = (
        crosswalk_with_subplant_ids["generator_operating_date"]
        .fillna(crosswalk_with_subplant_ids["generator_retirement_date"])
        .fillna(
            crosswalk_with_subplant_ids["original_planned_generator_operating_date"]
        )
    )

    # sort values to ensure static order
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.sort_values(
        by=[
            "plant_id_eia",
            "earliest_report_date",
            "sort_date",
            "generator_id",
            "emissions_unit_id_epa",
        ],
        ascending=True,
    ).copy()

    # update the subplant ids for each plant
    subplant_crosswalk_complete = update_subplant_ids(crosswalk_with_subplant_ids)

    subplant_crosswalk_complete = manually_update_subplant_id(
        subplant_crosswalk_complete
    )

    # Update the subplant ID to use the new_subplant.
    # we will also update the subplant_id to not use zero-indexed IDs, and instead
    # start all subplant_id at 1
    subplant_crosswalk_complete["subplant_id"] = (
        subplant_crosswalk_complete["new_subplant"] + 1
    )
    subplant_crosswalk_complete = subplant_crosswalk_complete.reset_index(drop=True)[
        [
            "plant_id_epa",
            "emissions_unit_id_epa",
            "plant_id_eia",
            "generator_id",
            "subplant_id",
            "unit_id_pudl",
            "generator_operating_date",
            "generator_retirement_date",
            "current_planned_generator_operating_date",
            "prime_mover_code",
        ]
    ].drop_duplicates(
        subset=[
            "plant_id_epa",
            "emissions_unit_id_epa",
            "plant_id_eia",
            "generator_id",
            "subplant_id",
        ],
        keep="last",
    )

    validation.check_for_1_to_many_subplant_mappings(
        subplant_crosswalk_complete, "generator_id"
    )
    validation.check_for_1_to_many_subplant_mappings(
        subplant_crosswalk_complete, "emissions_unit_id_epa"
    )

    # validate that there are no orphaned combined cycle plant parts in a subplant
    validation.check_for_orphaned_cc_part_in_subplant(
        subplant_crosswalk_complete, latest_validated_year
    )

    return subplant_crosswalk_complete


def manually_update_subplant_id(subplant_crosswalk: pd.DataFrame) -> pd.DataFrame:
    """This function corrects subplant mappings not caught by update_subplant_id. This
    is temporary until the pudl subplant crosswalk includes boiler-generator id matches.
    """

    # set all generators in plant 1391 to the same subplant
    subplant_crosswalk.loc[
        subplant_crosswalk["plant_id_eia"] == 1391, "subplant_id"
    ] = 0

    return subplant_crosswalk


def update_subplant_ids(subplant_crosswalk: pd.DataFrame) -> pd.DataFrame:
    """Ensures a complete and accurate subplant_id mapping for all generators.

    NOTE:
        1. This function is a temporary placeholder until the
        `pudl.analysis.epacamd_eia` code is updated.

    Data Preparation
        Because the existing subplant_id crosswalk was only meant to map CAMD units to
        EIA generators, it is missing a large number of subplant_ids for generators
        that do not report to CEMS. Before applying this function to the subplant
        crosswalk, the crosswalk must be completed with all generators by outer merging
        in the complete list of generators from EIA-860 (specifically the gens_eia860
        table from pudl). This dataframe also contains the complete list of
        `unit_id_pudl` mappings that will be necessary.

    High-level overview of method:
        1. Use the PUDL subplant_id if available. In the case where a unit_id_pudl
        groups several subplants, we overwrite these multiple existing subplant_id with
        a single subplant_id.
        2. Where there is no PUDL subplant_id, we use the unit_id_pudl to assign a
        unique subplant_id
        3. Where there is neither a pudl subplant_id nor unit_id_pudl, we use the
        generator ID to assign a unique subplant_id
        4. All of the new unique ids are renumbered in consecutive ascending order


    Detailed explanation of steps:
        1. Because the current subplant_id code does not take boiler-generator
        associations into account, there may be instances where the code assigns
        generators to different subplants when in fact, according to the
        boiler-generator association table, these generators are grouped into a single
        unit based on their boiler associations. The first step of this function is
        thus to identify if multiple subplant_id have been assigned to a single
        unit_id_pudl. If so, we replace the existing subplant_ids with a single
        subplant_id. For example, if a generator A was assigned subplant_id 0 and
        generator B was assigned subplant_id 1, but both generators A and B are part of
        unit_id_pudl 1, we would re-assign the subplant_id to both generators to 0 (we
        always use the lowest number subplant_id in each unit_id_pudl group). This may
        result in some subplant_id being skipped, but this is okay because we will
        later renumber all subplant ids (i.e. if there were also a generator C with
        subplant_id 2, there would no be no subplant_id 1 at the plant). Likewise,
        sometimes multiple unit_id_pudl are connected to a single subplant_id, so we
        also correct the unit_id_pudl basedon these connections.
        2. The second issue is that there are many NA subplant_id that we should fill.
        To do this, we first look at unit_id_pudl. If a group of generators are
        assigned a unit_id_pudl but have NA subplant_ids, we assign a single new
        subplant_id to this group of generators. If there are still generators at a
        plant that have both NA subplant_id and NA unit_id_pudl, we for now assume that
        each of these generators consitutes its own subplant. We thus assign a unique
        subplant_id to each generator that is unique from any existing subplant_id
        already at the plant. In the case that there are multiple emissions_unit_id_epa
        at a plant that are not matched to any other identifiers (generator_id,
        unit_id_pudl, or subplant_id), as is the case when there are units that report
        to CEMS but which do not exist in the EIA data, we assign these units to a single subplant.

        Args:
            subplant_crosswalk: a dataframe containing the output of
            `pudl.etl.glue_assets.make_subplant_ids` with
    """
    # Step 1: Correct unit_id_pudl using complete values loaded from the raw EIA-860
    # update unit_id_pudl using the unit_id_eia loaded from EIA-860
    subplant_crosswalk["unit_id_eia_numeric"] = pd.to_numeric(
        subplant_crosswalk["unit_id_eia_numeric"]
    )
    subplant_crosswalk = connect_ids(
        subplant_crosswalk,
        id_to_update="unit_id_pudl",
        connecting_id="unit_id_eia_numeric",
    )
    subplant_crosswalk["unit_id_pudl"] = subplant_crosswalk[
        "unit_id_pudl_connected_by_unit_id_eia_numeric"
    ]

    # Step 2: Create corrected versions of subplant_id and unit_id_pudl
    # if multiple unit_id_pudl are connected by a single subplant_id,
    # unit_id_pudl_connected groups these unit_id_pudl together
    subplant_crosswalk = connect_ids(
        subplant_crosswalk, id_to_update="unit_id_pudl", connecting_id="subplant_id"
    )

    # if multiple subplant_id are connected by a single unit_id_pudl, group these
    # subplant_id together
    subplant_crosswalk = connect_ids(
        subplant_crosswalk,
        id_to_update="subplant_id",
        connecting_id="unit_id_pudl_connected_by_subplant_id",
    )

    # Step 3: Fill missing subplant_id

    # create a numeric version of each generator_id
    subplant_crosswalk["numeric_generator_id"] = subplant_crosswalk.groupby(
        ["plant_id_eia"], dropna=False, sort=False
    )["generator_id"].transform(lambda x: pd.factorize(x, use_na_sentinel=False)[0])

    # when filling in missing subplant_id, we don't want these unit_id or
    # numeric_generator_id to overlap each other.
    # To ensure this, we will add 1000 to the unit_ids and 1000000 to the generator_id
    subplant_crosswalk["subplant_id_filled"] = (
        subplant_crosswalk[
            "subplant_id_connected_by_unit_id_pudl_connected_by_subplant_id"
        ]
        .fillna(subplant_crosswalk["unit_id_pudl_connected_by_subplant_id"] + 1000)
        .fillna(subplant_crosswalk["numeric_generator_id"] + 1000000)
    )

    # create a new unique subplant_id based on the connected subplant ids and the
    # filled unit_id
    subplant_crosswalk["new_subplant"] = subplant_crosswalk.groupby(
        ["plant_id_eia"], dropna=False, sort=False
    )["subplant_id_filled"].transform(
        lambda x: pd.factorize(x, use_na_sentinel=False)[0]
    )

    return subplant_crosswalk


def connect_ids(df, id_to_update, connecting_id):
    """Corrects an id value if it is connected by an id value in another column.

    if multiple subplant_id are connected by a single unit_id_pudl, this groups these
    subplant_id together
    if multiple unit_id_pudl are connected by a single subplant_id, this groups these
    unit_id_pudl together

    Args:
        df: dataframe containing columns with id_to_update and connecting_id columns
        subplant_unit_pairs
    """

    # get a table with all unique subplant to unit pairs
    subplant_unit_pairs = df[
        ["plant_id_eia", id_to_update, connecting_id]
    ].drop_duplicates()

    # identify if any non-NA id_to_update are duplicated, indicated that it is
    # associated with multiple connecting_id
    duplicates = subplant_unit_pairs[
        (
            subplant_unit_pairs.duplicated(
                subset=["plant_id_eia", connecting_id], keep=False
            )
        )
        & (~subplant_unit_pairs[connecting_id].isna())
    ].copy()

    # if there are any duplicate units, indicating an incorrect id_to_update,
    # fix the id_to_update
    df[f"{id_to_update}_connected_by_{connecting_id}"] = df[id_to_update]
    if len(duplicates) > 0:
        # find the lowest number subplant id associated with each duplicated
        # unit_id_pudl
        ids_to_replace = (
            duplicates.groupby(["plant_id_eia", connecting_id])[id_to_update]
            .min()
            .reset_index()
            .rename(columns={id_to_update: f"{id_to_update}_to_replace"})
        )

        # merge this replacement subplant_id into the dataframe and use it to update
        # the existing subplant id
        df = df.merge(
            ids_to_replace,
            how="left",
            on=["plant_id_eia", connecting_id],
            validate="m:1",
        )
        df.update(
            {
                f"{id_to_update}_connected_by_{connecting_id}": df[
                    f"{id_to_update}_to_replace"
                ]
            }
        )
        df = df.drop(columns=f"{id_to_update}_to_replace")

    return df
