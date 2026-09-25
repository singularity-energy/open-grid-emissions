import numpy as np
import pandas as pd

import oge.load_data as load_data
from oge.constants import ENERGY_STORAGE_PRIME_MOVERS
from oge.logging_util import get_logger

logger = get_logger(__name__)


# prime movers used to identify the storage type (standalone, co-located, or hybrid) of
# energy storage resources. CP is included here but not in ENERGY_STORAGE_PRIME_MOVERS,
# so that concentrated solar power resources keep their solar fuel category
STORAGE_TYPE_PRIME_MOVERS = ENERGY_STORAGE_PRIME_MOVERS + [
    "CP",  # Energy Storage, Concentrated Solar Power
]

# prime movers where energy storage can be an integral part of a generator that uses
# another energy source (e.g. compressed air storage that burns natural gas, or
# concentrated solar power with thermal storage). These are only considered hybrid if
# the generator reports an energy source code other than MWH
HYBRID_STORAGE_PRIME_MOVERS = ["CE", "CP"]

# the flags from the EIA-860 energy storage table used to identify storage types
STORAGE_FLAG_COLUMNS = [
    "is_dc_coupled_tightly",
    "is_direct_support",
    "served_co_located_renewable_firming",
    "is_independent",
]
DIRECT_SUPPORT_PLANT_COLUMNS = [
    "plant_id_eia_direct_support_1",
    "plant_id_eia_direct_support_2",
    "plant_id_eia_direct_support_3",
]

# the method used to assign each storage type, in the order in which the methods are
# applied (see assign_storage_type_to_generators())
STORAGE_TYPE_BY_METHOD = {
    "hybrid_prime_mover": "hybrid",
    "dc_coupled_tightly": "hybrid",
    "same_plant": "co_located",
    "direct_support_other_plant": "co_located",
    "same_location": "co_located",
    "eia860_flag": "co_located",
    "is_independent": "standalone",
    "no_evidence": "standalone",
}
STORAGE_TYPE_METHODS = list(STORAGE_TYPE_BY_METHOD)


def identify_energy_storage_types(
    primary_fuel_table: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Coordinating function for identifying the storage type of storage resources.

    Each energy storage generator in `primary_fuel_table` is assigned one of three
    storage types:

    - "standalone": the storage is an independent power plant with no non-storage
      generators at the same plant or location.
    - "co_located": the storage is at the same plant or location as a non-storage
      generator.
    - "hybrid": the storage is an integral part of a generator and is metered together
      with it.

    The storage type of each generator is then used to assign a `subplant_storage_type`
    to each subplant that contains a storage generator, and a `plant_storage_type` to
    each plant that contains a storage generator, regardless of the primary fuel of the
    subplant or plant.

    For storage that is co-located with a generator at a different plant (identified by
    the "direct_support_other_plant" or "same_location" rules), the `plant_id_eia` of
    the other plant(s) are recorded in `subplant_co_located_plant_ids` and
    `plant_co_located_plant_ids`, as a comma-separated string.

    Args:
        primary_fuel_table (pd.DataFrame): table of primary fuels by generator, with
            `plant_id_eia`, `subplant_id`, and `generator_id` columns.
        year (int): the data year.

    Returns:
        pd.DataFrame: `primary_fuel_table` with `subplant_storage_type`,
            `subplant_storage_type_method`, `subplant_co_located_plant_ids`,
            `plant_storage_type`, and `plant_co_located_plant_ids` columns added.

    Raises:
        ValueError: if storage generators at the same plant are assigned different
            storage types.
    """
    logger.info("Identifying energy storage types")

    # load EIA-860 attributes for each generator
    generators = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year=year,
        columns=[
            "plant_id_eia",
            "generator_id",
            "prime_mover_code",
            "energy_source_code_1",
            "operational_status",
        ],
    )
    plant_locations = load_data.load_pudl_table(
        "core_eia__entity_plants", columns=["plant_id_eia", "latitude", "longitude"]
    )
    generators = generators.merge(
        plant_locations, how="left", on="plant_id_eia", validate="m:1"
    )

    # a generator is considered to be operating if it is reported as existing in
    # EIA-860, or if it is in the primary fuel table (meaning that it reported data to
    # EIA-923 or CEMS). This includes generators that are producing energy while
    # testing before commercial operation (operational_status_code "TS")
    generators = generators.merge(
        primary_fuel_table[["plant_id_eia", "generator_id"]]
        .dropna(subset="generator_id")
        .drop_duplicates(),
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="1:1",
        indicator="in_primary_fuel_table",
    )
    generators["is_operating"] = (generators["operational_status"] == "existing") | (
        generators["in_primary_fuel_table"] == "both"
    )
    generators = generators[
        [
            "plant_id_eia",
            "generator_id",
            "prime_mover_code",
            "energy_source_code_1",
            "is_operating",
            "latitude",
            "longitude",
        ]
    ]

    # identify the storage generators in the primary fuel table
    storage_generators = (
        primary_fuel_table[["plant_id_eia", "subplant_id", "generator_id"]]
        .dropna(subset="generator_id")
        .drop_duplicates()
        .merge(
            generators[
                [
                    "plant_id_eia",
                    "generator_id",
                    "prime_mover_code",
                    "energy_source_code_1",
                ]
            ],
            how="left",
            on=["plant_id_eia", "generator_id"],
            validate="m:1",
        )
    )
    storage_generators = storage_generators[
        storage_generators["prime_mover_code"].isin(STORAGE_TYPE_PRIME_MOVERS)
    ]

    # the hybrid prime mover rule assumes that these generators use another energy
    # source (e.g. compressed air storage that burns natural gas). Warn if any report
    # MWH instead, since these may be newer technologies that do not fit this assumption
    hybrid_pm_reporting_mwh = storage_generators[
        storage_generators["prime_mover_code"].isin(HYBRID_STORAGE_PRIME_MOVERS)
        & (storage_generators["energy_source_code_1"] == "MWH")
    ]
    if len(hybrid_pm_reporting_mwh) > 0:
        logger.warning(
            "The following generators have a prime mover in "
            f"{HYBRID_STORAGE_PRIME_MOVERS} but report an energy source code of MWH, so "
            "they are not assumed to be hybrid storage. Check whether the storage type "
            "and fuel category assumptions for these technologies are still valid:\n"
            f"{hybrid_pm_reporting_mwh.to_string()}"
        )

    # add the storage flags reported in EIA-860. Pumped storage is not reported in the
    # energy storage table, so it will not have any flags
    storage_generators = storage_generators.merge(
        load_energy_storage_flags(year),
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )
    storage_generators[STORAGE_FLAG_COLUMNS] = (
        storage_generators[STORAGE_FLAG_COLUMNS]
        .astype("boolean")
        .fillna(False)
        .astype(bool)
    )

    storage_generators = assign_storage_type_to_generators(
        storage_generators, generators
    )
    validate_one_storage_type_per_plant(storage_generators)

    # the storage type method of each subplant is the method of the generator that was
    # assigned a storage type by the highest-priority rule
    storage_generators["rule_priority"] = storage_generators["storage_type_method"].map(
        {method: i for i, method in enumerate(STORAGE_TYPE_METHODS)}
    )
    subplant_storage_types = (
        storage_generators.sort_values("rule_priority")
        .drop_duplicates(subset=["plant_id_eia", "subplant_id"], keep="first")[
            ["plant_id_eia", "subplant_id", "storage_type", "storage_type_method"]
        ]
        .rename(
            columns={
                "storage_type": "subplant_storage_type",
                "storage_type_method": "subplant_storage_type_method",
            }
        )
    )
    plant_storage_types = (
        storage_generators[["plant_id_eia", "storage_type"]]
        .drop_duplicates()
        .rename(columns={"storage_type": "plant_storage_type"})
    )

    # identify the other plants that each storage subplant is co-located with, based on
    # the generators that were assigned the subplant's storage type method
    subplant_co_located_plants = (
        storage_generators.merge(
            subplant_storage_types[
                ["plant_id_eia", "subplant_id", "subplant_storage_type_method"]
            ],
            how="inner",
            left_on=["plant_id_eia", "subplant_id", "storage_type_method"],
            right_on=["plant_id_eia", "subplant_id", "subplant_storage_type_method"],
            validate="m:1",
        )
        .groupby(["plant_id_eia", "subplant_id"], dropna=False)["co_located_plant_ids"]
        .agg(lambda plant_lists: sorted(set().union(*plant_lists)))
        .reset_index()
    )
    plant_co_located_plants = (
        subplant_co_located_plants.groupby("plant_id_eia", dropna=False)[
            "co_located_plant_ids"
        ]
        .agg(lambda plant_lists: sorted(set().union(*plant_lists)))
        .reset_index()
    )
    subplant_storage_types = subplant_storage_types.merge(
        subplant_co_located_plants.assign(
            subplant_co_located_plant_ids=lambda df: df["co_located_plant_ids"].map(
                format_plant_ids
            )
        )[["plant_id_eia", "subplant_id", "subplant_co_located_plant_ids"]],
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="1:1",
    )
    plant_storage_types = plant_storage_types.merge(
        plant_co_located_plants.assign(
            plant_co_located_plant_ids=lambda df: df["co_located_plant_ids"].map(
                format_plant_ids
            )
        )[["plant_id_eia", "plant_co_located_plant_ids"]],
        how="left",
        on="plant_id_eia",
        validate="1:1",
    )

    primary_fuel_table = primary_fuel_table.merge(
        subplant_storage_types,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )
    primary_fuel_table = primary_fuel_table.merge(
        plant_storage_types, how="left", on="plant_id_eia", validate="m:1"
    )

    logger.info(
        "Storage types assigned to subplants:\n"
        + subplant_storage_types.groupby(
            ["subplant_storage_type", "subplant_storage_type_method"], dropna=False
        )
        .size()
        .to_string()
    )

    return primary_fuel_table


def load_energy_storage_flags(year: int) -> pd.DataFrame:
    """Loads the storage flags reported for each energy storage generator in EIA-860.

    Args:
        year (int): the data year.

    Returns:
        pd.DataFrame: table with one row per storage generator, with `plant_id_eia`,
            `generator_id`, the `STORAGE_FLAG_COLUMNS`, and the
            `DIRECT_SUPPORT_PLANT_COLUMNS`. The table is empty for years in which the
            energy storage table is not available.
    """
    storage_flags = load_data.load_pudl_table(
        "core_eia860__scd_generators_energy_storage",
        year=year,
        columns=["plant_id_eia", "generator_id"]
        + STORAGE_FLAG_COLUMNS
        + DIRECT_SUPPORT_PLANT_COLUMNS,
    )
    # flags are reported as 1 if true, and are otherwise either 0 or missing
    for column in STORAGE_FLAG_COLUMNS:
        storage_flags[column] = storage_flags[column].eq(1).fillna(False).astype(bool)

    return storage_flags


def assign_storage_type_to_generators(
    storage_generators: pd.DataFrame, generators: pd.DataFrame
) -> pd.DataFrame:
    """Assigns a storage type to each energy storage generator.

    Each generator is assigned the storage type of the first of the following rules
    that applies, and the rule is recorded in the `storage_type_method` column:

    1. "hybrid_prime_mover" (hybrid): the prime mover is in
       `HYBRID_STORAGE_PRIME_MOVERS` and the energy source code is not MWH.
    2. "dc_coupled_tightly" (hybrid): the storage is reported as tightly DC-coupled.
    3. "same_plant" (co_located): the plant has an operating non-storage generator.
    4. "direct_support_other_plant" (co_located): the storage is reported as directly
       supporting a generator at a different plant.
    5. "same_location" (co_located): a different plant with an operating non-storage
       generator has exactly the same latitude and longitude.
    6. "eia860_flag" (co_located): the storage is reported as directly supporting
       another generator or as firming co-located renewables.
    7. "is_independent" (standalone): the storage is reported as independent.
    8. "no_evidence" (standalone): none of the above rules apply.

    Args:
        storage_generators (pd.DataFrame): storage generators with `plant_id_eia`,
            `generator_id`, `prime_mover_code`, `energy_source_code_1`,
            `STORAGE_FLAG_COLUMNS`, and `DIRECT_SUPPORT_PLANT_COLUMNS` columns.
        generators (pd.DataFrame): EIA-860 attributes of all generators in the data
            year, with `plant_id_eia`, `generator_id`, `prime_mover_code`,
            `energy_source_code_1`, `is_operating`, `latitude`, and `longitude`
            columns.

    Returns:
        pd.DataFrame: `storage_generators` with `storage_type`, `storage_type_method`,
            and `co_located_plant_ids` columns added. `co_located_plant_ids` is a list
            of the `plant_id_eia` of other plants that the storage is co-located with,
            and is only filled for generators assigned by the
            "direct_support_other_plant" or "same_location" rules.
    """
    # identify plants and locations with operating non-storage generators
    operating_non_storage = generators[
        generators["is_operating"]
        & ~generators["prime_mover_code"].isin(STORAGE_TYPE_PRIME_MOVERS)
    ]
    non_storage_plants = operating_non_storage["plant_id_eia"].unique()
    non_storage_locations = (
        operating_non_storage.dropna(subset=["latitude", "longitude"])
        .groupby(["latitude", "longitude"], dropna=False)["plant_id_eia"]
        .unique()
    )

    # add the location of each storage generator
    storage_generators = storage_generators.merge(
        generators[["plant_id_eia", "generator_id", "latitude", "longitude"]],
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )

    # identify the plants other than its own that each storage generator directly
    # supports. Only consider the supported plants if the storage is reported as
    # providing direct support, since some generators report a supported plant without
    # doing so
    def get_directly_supported_plants(generator: dict) -> list[int]:
        if not generator["is_direct_support"]:
            return []
        return sorted(
            {
                int(generator[column])
                for column in DIRECT_SUPPORT_PLANT_COLUMNS
                if pd.notna(generator[column])
                and generator[column] != generator["plant_id_eia"]
            }
        )

    # identify other plants with operating non-storage generators that are at exactly
    # the same location as each storage generator
    def get_plants_at_same_location(generator: dict) -> list[int]:
        if pd.isna(generator["latitude"]) or pd.isna(generator["longitude"]):
            return []
        plants_at_location = non_storage_locations.get(
            (generator["latitude"], generator["longitude"]), []
        )
        return sorted(
            {
                int(plant)
                for plant in plants_at_location
                if plant != generator["plant_id_eia"]
            }
        )

    storage_generator_records = storage_generators.to_dict("records")
    directly_supported_plants = [
        get_directly_supported_plants(generator)
        for generator in storage_generator_records
    ]
    plants_at_same_location = [
        get_plants_at_same_location(generator)
        for generator in storage_generator_records
    ]
    supports_other_plant = pd.Series(
        [len(plants) > 0 for plants in directly_supported_plants],
        index=storage_generators.index,
        dtype=bool,
    )
    same_location = pd.Series(
        [len(plants) > 0 for plants in plants_at_same_location],
        index=storage_generators.index,
        dtype=bool,
    )

    # the conditions for each storage type method, in the order they are applied
    conditions = [
        storage_generators["prime_mover_code"].isin(HYBRID_STORAGE_PRIME_MOVERS)
        & (storage_generators["energy_source_code_1"] != "MWH"),
        storage_generators["is_dc_coupled_tightly"],
        storage_generators["plant_id_eia"].isin(non_storage_plants),
        supports_other_plant,
        same_location,
        storage_generators["is_direct_support"]
        | storage_generators["served_co_located_renewable_firming"],
        storage_generators["is_independent"],
    ]
    conditions = [
        condition.fillna(False).astype(bool).to_numpy() for condition in conditions
    ]
    storage_generators["storage_type_method"] = np.select(
        conditions, STORAGE_TYPE_METHODS[:-1], default=STORAGE_TYPE_METHODS[-1]
    )
    storage_generators["storage_type"] = storage_generators["storage_type_method"].map(
        STORAGE_TYPE_BY_METHOD
    )

    # for storage that is co-located with a generator at a different plant, record the
    # plant(s) that it is co-located with
    storage_generators["co_located_plant_ids"] = [
        supported
        if method == "direct_support_other_plant"
        else same_location_plants
        if method == "same_location"
        else []
        for method, supported, same_location_plants in zip(
            storage_generators["storage_type_method"],
            directly_supported_plants,
            plants_at_same_location,
        )
    ]

    return storage_generators.drop(columns=["latitude", "longitude"])


def format_plant_ids(plant_ids: list[int]) -> str | None:
    """Formats a list of plant IDs as a comma-separated string.

    Args:
        plant_ids (list[int]): a list of `plant_id_eia` values.

    Returns:
        str | None: the plant IDs separated by commas, or None if the list is empty.
    """
    if len(plant_ids) == 0:
        return None
    return ", ".join(str(plant_id) for plant_id in plant_ids)


def validate_one_storage_type_per_plant(storage_generators: pd.DataFrame) -> None:
    """Checks that all storage generators at each plant have the same storage type.

    Args:
        storage_generators (pd.DataFrame): storage generators with `plant_id_eia`,
            `generator_id`, `prime_mover_code`, `storage_type`, and
            `storage_type_method` columns.

    Raises:
        ValueError: if storage generators at the same plant have different storage
            types.
    """
    logger.info("Checking that each plant has a single storage type...  ")
    storage_types_per_plant = storage_generators.groupby("plant_id_eia", dropna=False)[
        "storage_type"
    ].transform("nunique")
    mismatched_plants = storage_generators[storage_types_per_plant > 1]
    if len(mismatched_plants) > 0:
        raise ValueError(
            "Storage generators at the following plants were assigned different "
            "storage types:\n"
            + mismatched_plants[
                [
                    "plant_id_eia",
                    "generator_id",
                    "prime_mover_code",
                    "storage_type",
                    "storage_type_method",
                ]
            ].to_string()
        )
    logger.info("OK")
