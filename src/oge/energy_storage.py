import numpy as np
import pandas as pd

import oge.load_data as load_data
import oge.validation as validation
from oge.column_checks import STORAGE_DATA_COLUMNS
from oge.constants import ENERGY_STORAGE_PRIME_MOVERS
from oge.logging_util import get_logger

logger = get_logger(__name__)


# prime movers where energy storage can be an integral part of a generator that uses
# another energy source (e.g. compressed air storage that burns natural gas). These are
# only considered hybrid if the generator reports an energy source code other than MWH.
# NOTE: concentrated solar power with thermal storage (CP) is not considered energy
# storage, since it stores solar thermal energy rather than using electricity as an
# input, so it is not included in ENERGY_STORAGE_PRIME_MOVERS or here
HYBRID_STORAGE_PRIME_MOVERS = ["CE"]

# the flags from the EIA-860 energy storage table used to identify storage categories
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

# the method used to assign each storage category, in the order in which the methods are
# applied (see assign_storage_category_to_generators())
STORAGE_CATEGORY_BY_METHOD = {
    "hybrid_prime_mover": "hybrid",
    "pumped_storage_with_inflow": "hybrid",
    "dc_coupled_tightly": "hybrid",
    "same_plant": "co_located",
    "direct_support_other_plant": "co_located",
    "same_location": "co_located",
    "eia860_flag": "co_located",
    "is_independent": "standalone",
    "no_evidence": "standalone",
}
STORAGE_CATEGORY_METHODS = list(STORAGE_CATEGORY_BY_METHOD)


def identify_energy_storage_categories(
    primary_fuel_table: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Coordinating function for identifying the storage category of storage resources.

    Each energy storage generator in `primary_fuel_table` is assigned one of three
    storage categories:

    - "standalone": the storage is an independent power plant with no non-storage
      generators at the same plant or location.
    - "co_located": the storage is at the same plant or location as a non-storage
      generator.
    - "hybrid": the storage is an integral part of a generator and is metered together
      with it.

    The storage category of each generator is then used to assign a `subplant_storage_category`
    to each subplant that contains a storage generator, and a `plant_storage_category` to
    each plant that contains a storage generator, regardless of the primary fuel of the
    subplant or plant.

    For storage that is co-located with a generator at a different plant (identified by
    the "direct_support_other_plant" or "same_location" rules), the `plant_id_eia` of
    the other plant(s) are recorded for each storage subplant in `co_located_plant_ids`,
    as a comma-separated string.

    Args:
        primary_fuel_table (pd.DataFrame): table of primary fuels by generator, with
            `plant_id_eia`, `subplant_id`, and `generator_id` columns.
        year (int): the data year.

    Returns:
        pd.DataFrame: `primary_fuel_table` with `subplant_storage_category`,
            `subplant_storage_category_method`, `co_located_plant_ids`, and
            `plant_storage_category` columns added.

    Raises:
        ValueError: if storage generators at the same plant are assigned different
            storage categories.
    """
    logger.info("Identifying energy storage categories")

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
        storage_generators["prime_mover_code"].isin(ENERGY_STORAGE_PRIME_MOVERS)
    ]

    # the hybrid prime mover rule assumes that these generators burn a fuel to supplement
    # the stored energy (e.g. compressed air storage that burns natural gas). Warn if any
    # report MWH instead, since these may be newer technologies that do not fit this
    # assumption
    hybrid_pm_reporting_mwh = storage_generators[
        storage_generators["prime_mover_code"].isin(HYBRID_STORAGE_PRIME_MOVERS)
        & (storage_generators["energy_source_code_1"] == "MWH")
    ]
    if len(hybrid_pm_reporting_mwh) > 0:
        logger.warning(
            "The following generators have a prime mover in "
            f"{HYBRID_STORAGE_PRIME_MOVERS} but report an energy source code of MWH, so "
            "they are not assumed to be hybrid storage. Check whether the storage category "
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

    # identify pumped storage that also generates electricity from natural inflow,
    # based on the charging and discharging reported in EIA-923
    storage_generators["is_pumped_storage_with_inflow"] = (
        storage_generators["prime_mover_code"] == "PS"
    ) & storage_generators["plant_id_eia"].isin(
        identify_pumped_storage_with_inflow(load_energy_storage_dispatch(year))
    )

    storage_generators = assign_storage_category_to_generators(
        storage_generators, generators
    )
    validate_one_storage_category_per_plant(storage_generators)

    # the storage category method of each subplant is the method of the generator that was
    # assigned a storage category by the highest-priority rule
    storage_generators["rule_priority"] = storage_generators[
        "storage_category_method"
    ].map({method: i for i, method in enumerate(STORAGE_CATEGORY_METHODS)})
    subplant_storage_categories = (
        storage_generators.sort_values("rule_priority")
        .drop_duplicates(subset=["plant_id_eia", "subplant_id"], keep="first")[
            [
                "plant_id_eia",
                "subplant_id",
                "storage_category",
                "storage_category_method",
            ]
        ]
        .rename(
            columns={
                "storage_category": "subplant_storage_category",
                "storage_category_method": "subplant_storage_category_method",
            }
        )
    )
    plant_storage_categories = (
        storage_generators[["plant_id_eia", "storage_category"]]
        .drop_duplicates()
        .rename(columns={"storage_category": "plant_storage_category"})
    )

    # identify the other plants that each storage subplant is co-located with, based on
    # the generators that were assigned the subplant's storage category method
    subplant_co_located_plants = (
        storage_generators.merge(
            subplant_storage_categories[
                ["plant_id_eia", "subplant_id", "subplant_storage_category_method"]
            ],
            how="inner",
            left_on=["plant_id_eia", "subplant_id", "storage_category_method"],
            right_on=[
                "plant_id_eia",
                "subplant_id",
                "subplant_storage_category_method",
            ],
            validate="m:1",
        )
        .groupby(["plant_id_eia", "subplant_id"], dropna=False)[
            "co_located_plant_id_list"
        ]
        .agg(lambda plant_lists: sorted(set().union(*plant_lists)))
        .reset_index()
    )
    subplant_storage_categories = subplant_storage_categories.merge(
        subplant_co_located_plants.assign(
            co_located_plant_ids=lambda df: df["co_located_plant_id_list"].map(
                format_plant_ids
            )
        )[["plant_id_eia", "subplant_id", "co_located_plant_ids"]],
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="1:1",
    )

    primary_fuel_table = primary_fuel_table.merge(
        subplant_storage_categories,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )
    primary_fuel_table = primary_fuel_table.merge(
        plant_storage_categories, how="left", on="plant_id_eia", validate="m:1"
    )

    logger.info(
        "Storage categories assigned to subplants:\n"
        + subplant_storage_categories.groupby(
            ["subplant_storage_category", "subplant_storage_category_method"],
            dropna=False,
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


def assign_storage_category_to_generators(
    storage_generators: pd.DataFrame, generators: pd.DataFrame
) -> pd.DataFrame:
    """Assigns a storage category to each energy storage generator.

    Each generator is assigned the storage category of the first of the following rules
    that applies, and the rule is recorded in the `storage_category_method` column:

    1. "hybrid_prime_mover" (hybrid): the prime mover is in
       `HYBRID_STORAGE_PRIME_MOVERS` and the energy source code is not MWH.
    2. "pumped_storage_with_inflow" (hybrid): the generator is pumped storage that
       also generates electricity from natural inflow (see
       `identify_pumped_storage_with_inflow()`).
    3. "dc_coupled_tightly" (hybrid): the storage is reported as tightly DC-coupled.
    4. "same_plant" (co_located): the plant has an operating non-storage generator.
    5. "direct_support_other_plant" (co_located): the storage is reported as directly
       supporting a generator at a different plant.
    6. "same_location" (co_located): a different plant with an operating non-storage
       generator has exactly the same latitude and longitude.
    7. "eia860_flag" (co_located): the storage is reported as directly supporting
       another generator or as firming co-located renewables.
    8. "is_independent" (standalone): the storage is reported as independent.
    9. "no_evidence" (standalone): none of the above rules apply.

    Args:
        storage_generators (pd.DataFrame): storage generators with `plant_id_eia`,
            `generator_id`, `prime_mover_code`, `energy_source_code_1`,
            `is_pumped_storage_with_inflow`, `STORAGE_FLAG_COLUMNS`, and
            `DIRECT_SUPPORT_PLANT_COLUMNS` columns.
        generators (pd.DataFrame): EIA-860 attributes of all generators in the data
            year, with `plant_id_eia`, `generator_id`, `prime_mover_code`,
            `energy_source_code_1`, `is_operating`, `latitude`, and `longitude`
            columns.

    Returns:
        pd.DataFrame: `storage_generators` with `storage_category`, `storage_category_method`,
            and `co_located_plant_id_list` columns added. `co_located_plant_id_list` is a list
            of the `plant_id_eia` of other plants that the storage is co-located with,
            and is only filled for generators assigned by the
            "direct_support_other_plant" or "same_location" rules.
    """
    # identify plants and locations with operating non-storage generators
    operating_non_storage = generators[
        generators["is_operating"]
        & ~generators["prime_mover_code"].isin(ENERGY_STORAGE_PRIME_MOVERS)
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

    # the conditions for each storage category method, in the order they are applied
    conditions = [
        storage_generators["prime_mover_code"].isin(HYBRID_STORAGE_PRIME_MOVERS)
        & (storage_generators["energy_source_code_1"] != "MWH"),
        storage_generators["is_pumped_storage_with_inflow"],
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
    storage_generators["storage_category_method"] = np.select(
        conditions, STORAGE_CATEGORY_METHODS[:-1], default=STORAGE_CATEGORY_METHODS[-1]
    )
    storage_generators["storage_category"] = storage_generators[
        "storage_category_method"
    ].map(STORAGE_CATEGORY_BY_METHOD)

    # for storage that is co-located with a generator at a different plant, record the
    # plant(s) that it is co-located with
    storage_generators["co_located_plant_id_list"] = [
        supported
        if method == "direct_support_other_plant"
        else same_location_plants
        if method == "same_location"
        else []
        for method, supported, same_location_plants in zip(
            storage_generators["storage_category_method"],
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


def validate_one_storage_category_per_plant(storage_generators: pd.DataFrame) -> None:
    """Checks that all storage generators at each plant have the same storage category.

    Args:
        storage_generators (pd.DataFrame): storage generators with `plant_id_eia`,
            `generator_id`, `prime_mover_code`, `storage_category`, and
            `storage_category_method` columns.

    Raises:
        ValueError: if storage generators at the same plant have different storage
            types.
    """
    logger.info("Checking that each plant has a single storage category...  ")
    storage_categories_per_plant = storage_generators.groupby(
        "plant_id_eia", dropna=False
    )["storage_category"].transform("nunique")
    mismatched_plants = storage_generators[storage_categories_per_plant > 1]
    if len(mismatched_plants) > 0:
        raise ValueError(
            "Storage generators at the following plants were assigned different "
            "storage categories:\n"
            + mismatched_plants[
                [
                    "plant_id_eia",
                    "generator_id",
                    "prime_mover_code",
                    "storage_category",
                    "storage_category_method",
                ]
            ].to_string()
        )
    logger.info("OK")


def create_monthly_energy_storage_data(
    primary_fuel_table: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Coordinating function for creating monthly energy storage charging and discharging data.

    Monthly charging and discharging data for each energy storage resource is loaded
    from EIA-923 (see `load_energy_storage_dispatch()`), checked for data quality
    issues, and allocated to each storage subplant. The discharge of pumped storage
    that also generates electricity from natural inflow (subplants with a
    `subplant_storage_category_method` of "pumped_storage_with_inflow") is left blank,
    since it cannot be separated from the electricity generated from natural inflow.

    Args:
        primary_fuel_table (pd.DataFrame): table of primary fuels by generator, with
            `plant_id_eia`, `subplant_id`, `generator_id`, and
            `subplant_storage_category_method` columns.
        year (int): the data year.

    Returns:
        pd.DataFrame: table with one row per storage subplant-month, with
            `plant_id_eia`, `subplant_id`, `report_date`, and `STORAGE_DATA_COLUMNS`.
    """
    logger.info("Creating monthly energy storage charging and discharging data")
    storage_dispatch = load_energy_storage_dispatch(year)
    log_energy_storage_data_quality(storage_dispatch)

    monthly_storage_data = allocate_energy_storage_dispatch_to_subplants(
        storage_dispatch, primary_fuel_table, year
    )

    # the discharge of pumped storage with natural inflow cannot be separated from the
    # electricity generated from natural inflow, so leave it blank
    pumped_storage_with_inflow = primary_fuel_table.loc[
        primary_fuel_table["subplant_storage_category_method"]
        == "pumped_storage_with_inflow",
        ["plant_id_eia", "subplant_id"],
    ].drop_duplicates()
    monthly_storage_data = monthly_storage_data.merge(
        pumped_storage_with_inflow,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
        indicator="pumped_storage_with_inflow",
    )
    monthly_storage_data.loc[
        monthly_storage_data["pumped_storage_with_inflow"] == "both",
        "storage_discharge_mwh",
    ] = np.nan
    monthly_storage_data = monthly_storage_data.drop(
        columns="pumped_storage_with_inflow"
    )

    return monthly_storage_data


def load_energy_storage_dispatch(year: int) -> pd.DataFrame:
    """Loads the monthly charging and discharging of each energy storage resource.

    Charging and discharging data is reported in EIA-923 for each plant, prime mover,
    and energy source code. Only energy storage prime movers are kept. For storage
    reported in MWh, charging is the reported fuel consumed for electricity, and
    discharging is the reported gross generation. For storage that supplements the
    stored energy with a combustion fuel (prime movers in
    `HYBRID_STORAGE_PRIME_MOVERS` that report a fuel other than MWh, e.g. compressed air
    storage that burns natural gas), discharging is the reported gross generation and
    charging is the gross generation minus the net generation.

    If a month reports zero discharge, but filling the discharge with the net
    generation plus the charge would make that month consistent, and all other months
    for that resource are consistent, the discharge is filled with this value.

    Args:
        year (int): the data year.

    Returns:
        pd.DataFrame: table with one row per plant, prime mover, energy source code, and
            month, with `STORAGE_DATA_COLUMNS`, `net_generation_mwh`, and
            `discharge_filled` columns. The table is empty for years in which the
            energy storage table is not available.
    """
    storage_dispatch = load_data.load_pudl_table(
        "core_eia923__monthly_energy_storage",
        year=year,
        columns=[
            "plant_id_eia",
            "report_date",
            "prime_mover_code",
            "energy_source_code",
            "fuel_units",
            "fuel_consumed_for_electricity_units",
            "gross_generation_mwh",
            "net_generation_mwh",
        ],
    )
    storage_dispatch = storage_dispatch[
        storage_dispatch["prime_mover_code"].isin(ENERGY_STORAGE_PRIME_MOVERS)
    ].copy()

    # identify how the charging of each resource is reported
    reported_in_mwh = storage_dispatch["fuel_units"] == "mwh"
    supplemented_with_fuel = (
        storage_dispatch["prime_mover_code"].isin(HYBRID_STORAGE_PRIME_MOVERS)
        & ~reported_in_mwh
    )
    unexpected_units = storage_dispatch[~reported_in_mwh & ~supplemented_with_fuel]
    if len(unexpected_units) > 0:
        logger.warning(
            "The following energy storage data is not reported in MWh, so charging and "
            "discharging data will not be included for these resources:\n"
            + validation.limit_error_output_df(unexpected_units).to_string()
        )
    storage_dispatch = storage_dispatch[reported_in_mwh | supplemented_with_fuel]
    supplemented_with_fuel = supplemented_with_fuel[storage_dispatch.index]

    storage_dispatch["storage_discharge_mwh"] = storage_dispatch["gross_generation_mwh"]
    storage_dispatch["storage_charge_mwh"] = storage_dispatch[
        "fuel_consumed_for_electricity_units"
    ].where(
        ~supplemented_with_fuel,
        storage_dispatch["gross_generation_mwh"]
        - storage_dispatch["net_generation_mwh"],
    )

    # fill months that report zero discharge if doing so makes the month consistent
    # and all other months for that resource are consistent
    is_inconsistent = identify_inconsistent_net_generation(storage_dispatch)
    is_fillable = (
        is_inconsistent
        & (storage_dispatch["storage_discharge_mwh"] == 0)
        & (
            storage_dispatch["net_generation_mwh"]
            + storage_dispatch["storage_charge_mwh"]
            > 0
        )
    )
    resource_keys = [
        storage_dispatch[column]
        for column in ["plant_id_eia", "prime_mover_code", "energy_source_code"]
    ]
    inconsistent_months = is_inconsistent.groupby(
        resource_keys, dropna=False
    ).transform("sum")
    fillable_months = is_fillable.groupby(resource_keys, dropna=False).transform("sum")
    storage_dispatch["discharge_filled"] = is_fillable & (
        inconsistent_months == fillable_months
    )
    storage_dispatch.loc[
        storage_dispatch["discharge_filled"], "storage_discharge_mwh"
    ] = storage_dispatch["net_generation_mwh"] + storage_dispatch["storage_charge_mwh"]

    return storage_dispatch[
        [
            "plant_id_eia",
            "report_date",
            "prime_mover_code",
            "energy_source_code",
            "net_generation_mwh",
            "discharge_filled",
        ]
        + STORAGE_DATA_COLUMNS
    ].reset_index(drop=True)


def identify_inconsistent_net_generation(storage_dispatch: pd.DataFrame) -> pd.Series:
    """Identifies records where net generation does not equal discharge minus charge.

    Args:
        storage_dispatch (pd.DataFrame): table with `net_generation_mwh` and
            `STORAGE_DATA_COLUMNS`.

    Returns:
        pd.Series: boolean series that is True where the net generation differs from the
            discharge minus the charge by more than 1 MWh, or where any of these values
            are missing.
    """
    is_consistent = np.isclose(
        storage_dispatch["storage_discharge_mwh"]
        - storage_dispatch["storage_charge_mwh"],
        storage_dispatch["net_generation_mwh"],
        rtol=0,
        atol=1.0,
    )
    return pd.Series(~is_consistent, index=storage_dispatch.index)


def identify_pumped_storage_with_inflow(storage_dispatch: pd.DataFrame) -> list[int]:
    """Identifies pumped storage plants that also generate electricity from inflow.

    Pumped storage that only discharges energy that it previously pumped will always
    discharge less energy than it charges, due to round-trip efficiency losses. Pumped
    storage plants that discharge at least as much energy as they charge over the year
    (and discharge some energy) are therefore assumed to also generate electricity from
    natural inflow to the upper reservoir.

    Args:
        storage_dispatch (pd.DataFrame): monthly charging and discharging data, from
            `load_energy_storage_dispatch()`.

    Returns:
        list[int]: the `plant_id_eia` of pumped storage plants with natural inflow.
    """
    annual_pumped_storage = (
        storage_dispatch[storage_dispatch["prime_mover_code"] == "PS"]
        .groupby("plant_id_eia", dropna=False)[STORAGE_DATA_COLUMNS]
        .sum(min_count=1)
        .reset_index()
    )
    return annual_pumped_storage.loc[
        (
            annual_pumped_storage["storage_discharge_mwh"]
            >= annual_pumped_storage["storage_charge_mwh"]
        )
        & (annual_pumped_storage["storage_discharge_mwh"] > 0),
        "plant_id_eia",
    ].tolist()


def log_energy_storage_data_quality(storage_dispatch: pd.DataFrame) -> None:
    """Logs data quality issues in the monthly energy storage data.

    Args:
        storage_dispatch (pd.DataFrame): monthly charging and discharging data, from
            `load_energy_storage_dispatch()`.
    """
    resource_keys = ["plant_id_eia", "prime_mover_code", "energy_source_code"]

    filled_months = storage_dispatch[storage_dispatch["discharge_filled"]]
    if len(filled_months) > 0:
        logger.info(
            "Filled missing energy storage discharge using net generation plus charge "
            "for the following months:\n"
            + validation.limit_error_output_df(filled_months).to_string()
        )

    inconsistent_months = storage_dispatch[
        identify_inconsistent_net_generation(storage_dispatch)
    ]
    if len(inconsistent_months) > 0:
        inconsistent_summary = (
            inconsistent_months.assign(
                difference_mwh=lambda df: (
                    df["storage_discharge_mwh"]
                    - df["storage_charge_mwh"]
                    - df["net_generation_mwh"]
                )
            )
            .groupby(resource_keys, dropna=False)
            .agg(
                inconsistent_months=("report_date", "count"),
                difference_mwh=("difference_mwh", "sum"),
            )
            .reset_index()
        )
        logger.warning(
            "Net generation does not equal energy storage discharge minus charge for "
            "the following resources:\n"
            + validation.limit_error_output_df(inconsistent_summary).to_string()
        )

    annual_dispatch = (
        storage_dispatch.groupby(resource_keys, dropna=False)[STORAGE_DATA_COLUMNS]
        .sum(min_count=1)
        .reset_index()
    )
    pumped_storage_with_inflow = annual_dispatch[
        annual_dispatch["plant_id_eia"].isin(
            identify_pumped_storage_with_inflow(storage_dispatch)
        )
        & (annual_dispatch["prime_mover_code"] == "PS")
    ]
    if len(pumped_storage_with_inflow) > 0:
        logger.info(
            "The following pumped storage plants discharge at least as much energy as "
            "they charge, so are assumed to generate electricity from natural inflow. "
            "Their discharge data will be left blank:\n"
            + validation.limit_error_output_df(pumped_storage_with_inflow).to_string()
        )
    discharge_exceeds_charge = annual_dispatch[
        (annual_dispatch["prime_mover_code"] != "PS")
        & (
            annual_dispatch["storage_discharge_mwh"]
            > annual_dispatch["storage_charge_mwh"]
        )
    ]
    if len(discharge_exceeds_charge) > 0:
        logger.warning(
            "The following energy storage resources report more annual discharge than "
            "charge:\n"
            + validation.limit_error_output_df(discharge_exceeds_charge).to_string()
        )


def allocate_energy_storage_dispatch_to_subplants(
    storage_dispatch: pd.DataFrame, primary_fuel_table: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Allocates plant-level energy storage charging and discharging to subplants.

    The charging and discharging reported for each plant and prime mover is allocated
    to each storage generator with that prime mover at the plant based on its share of
    nameplate capacity. Generator-level data is not reported for energy storage in
    EIA-923, so this matches how net generation is allocated to these generators. If
    the generators at a plant do not report any nameplate capacity, the data is
    allocated equally.

    Args:
        storage_dispatch (pd.DataFrame): monthly charging and discharging data, from
            `load_energy_storage_dispatch()`.
        primary_fuel_table (pd.DataFrame): table of primary fuels by generator, used to
            identify the generators and subplant of each plant.
        year (int): the data year.

    Returns:
        pd.DataFrame: table with one row per storage subplant-month, with
            `plant_id_eia`, `subplant_id`, `report_date`, and `STORAGE_DATA_COLUMNS`.
    """
    # identify the storage generators at each plant, and their share of capacity
    storage_generators = (
        primary_fuel_table[["plant_id_eia", "subplant_id", "generator_id"]]
        .dropna(subset="generator_id")
        .drop_duplicates()
        .merge(
            load_data.load_pudl_table(
                "core_eia860__scd_generators",
                year=year,
                columns=[
                    "plant_id_eia",
                    "generator_id",
                    "prime_mover_code",
                    "capacity_mw",
                ],
            ),
            how="left",
            on=["plant_id_eia", "generator_id"],
            validate="m:1",
        )
    )
    storage_generators = storage_generators[
        storage_generators["prime_mover_code"].isin(ENERGY_STORAGE_PRIME_MOVERS)
    ].copy()
    storage_generators["capacity_mw"] = storage_generators["capacity_mw"].fillna(0)
    plant_capacity = storage_generators.groupby(
        ["plant_id_eia", "prime_mover_code"], dropna=False
    )["capacity_mw"].transform("sum")
    generators_at_plant = storage_generators.groupby(
        ["plant_id_eia", "prime_mover_code"], dropna=False
    )["generator_id"].transform("count")
    storage_generators["allocation_share"] = np.where(
        plant_capacity > 0,
        storage_generators["capacity_mw"] / plant_capacity,
        1 / generators_at_plant,
    )

    # allocate the data to each generator. A plant and prime mover can have multiple
    # generators, and can report data for multiple energy source codes
    allocated_dispatch = storage_dispatch.merge(
        storage_generators[
            [
                "plant_id_eia",
                "subplant_id",
                "generator_id",
                "prime_mover_code",
                "allocation_share",
            ]
        ],
        how="left",
        on=["plant_id_eia", "prime_mover_code"],
        validate="m:m",
    )
    # warn about any non-zero data that could not be allocated to a generator
    unallocated_dispatch = (
        allocated_dispatch[allocated_dispatch["generator_id"].isna()]
        .groupby(["plant_id_eia", "prime_mover_code"], dropna=False)[
            STORAGE_DATA_COLUMNS
        ]
        .sum(min_count=1)
        .reset_index()
    )
    unallocated_dispatch = unallocated_dispatch[
        (unallocated_dispatch[STORAGE_DATA_COLUMNS].fillna(0) != 0).any(axis=1)
    ]
    if len(unallocated_dispatch) > 0:
        logger.warning(
            "Energy storage data for the following plants could not be allocated to a "
            "storage generator, and will not be included in the results:\n"
            + validation.limit_error_output_df(unallocated_dispatch).to_string()
        )
    allocated_dispatch = allocated_dispatch.dropna(subset="generator_id")
    for column in STORAGE_DATA_COLUMNS:
        allocated_dispatch[column] = (
            allocated_dispatch[column] * allocated_dispatch["allocation_share"]
        )

    # aggregate the data to each subplant
    monthly_storage_data = (
        allocated_dispatch.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )[STORAGE_DATA_COLUMNS]
        .sum(min_count=1)
        .round(1)
        .reset_index()
    )

    return monthly_storage_data


def add_monthly_energy_storage_data(
    monthly_subplant_data: pd.DataFrame, monthly_storage_data: pd.DataFrame
) -> pd.DataFrame:
    """Adds monthly energy storage charging and discharging data to subplant data.

    Storage data for a subplant-month that is not in `monthly_subplant_data` is added
    as a new row if the subplant is in `monthly_subplant_data` for other months. This
    can happen when there is no other data for the subplant in that month (for example,
    when CEMS data for the month is removed because it is all zero). The other data
    columns for these rows are left blank. Storage data for subplants that are not in
    `monthly_subplant_data` at all is not added.

    Args:
        monthly_subplant_data (pd.DataFrame): combined monthly data for all subplants,
            with one row per subplant-month.
        monthly_storage_data (pd.DataFrame): monthly energy storage charging and
            discharging data, from `create_monthly_energy_storage_data()`.

    Returns:
        pd.DataFrame: `monthly_subplant_data` with `STORAGE_DATA_COLUMNS` added. These
            columns are blank for subplants that are not energy storage.
    """
    subplant_month_keys = ["plant_id_eia", "subplant_id", "report_date"]

    # identify storage data for subplant-months that are not in the subplant data
    monthly_storage_data = monthly_storage_data.merge(
        monthly_subplant_data[subplant_month_keys],
        how="left",
        on=subplant_month_keys,
        validate="1:1",
        indicator="subplant_month_in_data",
    )
    missing_months = monthly_storage_data[
        monthly_storage_data["subplant_month_in_data"] == "left_only"
    ].merge(
        monthly_subplant_data[["plant_id_eia", "subplant_id"]].drop_duplicates(),
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
        indicator="subplant_in_data",
    )

    # add storage data for missing months of subplants that are in the subplant data
    rows_to_add = missing_months.loc[
        missing_months["subplant_in_data"] == "both",
        subplant_month_keys + STORAGE_DATA_COLUMNS,
    ]
    if len(rows_to_add) > 0:
        logger.info(
            f"Adding {len(rows_to_add)} subplant-months that only contain energy "
            "storage data to the subplant data"
        )

    # warn about any non-zero storage data for subplants that are not in the data
    unmatched_storage_data = (
        missing_months[missing_months["subplant_in_data"] == "left_only"]
        .groupby(["plant_id_eia", "subplant_id"], dropna=False)[STORAGE_DATA_COLUMNS]
        .sum(min_count=1)
        .reset_index()
    )
    unmatched_storage_data = unmatched_storage_data[
        (unmatched_storage_data[STORAGE_DATA_COLUMNS].fillna(0) != 0).any(axis=1)
    ]
    if len(unmatched_storage_data) > 0:
        logger.warning(
            "Energy storage data for the following subplants could not be matched to "
            "subplant data, and will not be included in the results:\n"
            + validation.limit_error_output_df(unmatched_storage_data).to_string()
        )

    monthly_subplant_data = monthly_subplant_data.merge(
        monthly_storage_data[subplant_month_keys + STORAGE_DATA_COLUMNS],
        how="left",
        on=subplant_month_keys,
        validate="1:1",
    )
    monthly_subplant_data = pd.concat(
        [monthly_subplant_data, rows_to_add], axis=0, ignore_index=True
    )

    return monthly_subplant_data
