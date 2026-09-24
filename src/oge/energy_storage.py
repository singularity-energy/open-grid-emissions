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

# prime movers where energy storage is an integral part of a generator that is metered
# together with the storage
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

    Args:
        primary_fuel_table (pd.DataFrame): table of primary fuels by generator, with
            `plant_id_eia`, `subplant_id`, and `generator_id` columns.
        year (int): the data year.

    Returns:
        pd.DataFrame: `primary_fuel_table` with `subplant_storage_type`,
            `subplant_storage_type_method`, and `plant_storage_type` columns added.

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
            "operational_status",
        ],
    )
    plant_locations = load_data.load_pudl_table(
        "core_eia__entity_plants", columns=["plant_id_eia", "latitude", "longitude"]
    )
    generators = generators.merge(
        plant_locations, how="left", on="plant_id_eia", validate="m:1"
    )
    generators = generators[
        [
            "plant_id_eia",
            "generator_id",
            "prime_mover_code",
            "operational_status",
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
            generators[["plant_id_eia", "generator_id", "prime_mover_code"]],
            how="left",
            on=["plant_id_eia", "generator_id"],
            validate="m:1",
        )
    )
    storage_generators = storage_generators[
        storage_generators["prime_mover_code"].isin(STORAGE_TYPE_PRIME_MOVERS)
    ]

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
       `HYBRID_STORAGE_PRIME_MOVERS`.
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
            `generator_id`, `prime_mover_code`, `STORAGE_FLAG_COLUMNS`, and
            `DIRECT_SUPPORT_PLANT_COLUMNS` columns.
        generators (pd.DataFrame): EIA-860 attributes of all generators in the data
            year, with `plant_id_eia`, `generator_id`, `prime_mover_code`,
            `operational_status`, `latitude`, and `longitude` columns.

    Returns:
        pd.DataFrame: `storage_generators` with `storage_type` and
            `storage_type_method` columns added.
    """
    # identify plants and locations with operating non-storage generators
    operating_non_storage = generators[
        (generators["operational_status"] == "existing")
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

    # identify storage generators that directly support a generator at a different
    # plant. Only consider the supported plants if the storage is reported as providing
    # direct support, since some generators report a supported plant without doing so
    supports_other_plant = pd.Series(False, index=storage_generators.index)
    for column in DIRECT_SUPPORT_PLANT_COLUMNS:
        supports_other_plant = supports_other_plant | (
            storage_generators[column].notna()
            & (storage_generators[column] != storage_generators["plant_id_eia"])
        ).fillna(False)
    supports_other_plant = (
        supports_other_plant & storage_generators["is_direct_support"]
    )

    # identify storage generators at the same location as a different plant with
    # operating non-storage generators
    def shares_location_with_non_storage_plant(row: pd.Series) -> bool:
        if pd.isna(row["latitude"]) or pd.isna(row["longitude"]):
            return False
        plants_at_location = non_storage_locations.get(
            (row["latitude"], row["longitude"]), []
        )
        return any(plant != row["plant_id_eia"] for plant in plants_at_location)

    same_location = storage_generators.apply(
        shares_location_with_non_storage_plant, axis=1
    ).astype(bool)

    # the conditions for each storage type method, in the order they are applied
    conditions = [
        storage_generators["prime_mover_code"].isin(HYBRID_STORAGE_PRIME_MOVERS),
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

    return storage_generators.drop(columns=["latitude", "longitude"])


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
