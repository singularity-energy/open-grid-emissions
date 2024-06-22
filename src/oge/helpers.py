import numpy as np
import pandas as pd

from oge.column_checks import get_dtypes, apply_dtypes
from oge.constants import earliest_data_year, latest_validated_year
from oge.filepaths import reference_table_folder, outputs_folder
import oge.load_data as load_data
from oge.logging_util import get_logger
import oge.validation as validation

logger = get_logger(__name__)


def create_plant_attributes_table(
    cems: pd.DataFrame,
    eia923_allocated: pd.DataFrame,
    year: int,
    primary_fuel_table: pd.DataFrame,
) -> pd.DataFrame:
    """Creates the plant attributes table.

    Args:
        cems (pd.DataFrame): CEMS table.
        eia923_allocated (pd.DataFrame): allocated EIA-923 data.
        year (int): a four-digit year
        primary_fuel_table (pd.DataFrame): primary fuel table.

    Returns:
        pd.DataFrame: the plants attributes table. Timezone, geographical and fuel
            information can be found in this table.
    """
    # create a table with the unique plantids from both dataframes
    eia_plants = eia923_allocated.copy()[
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

    # calculate the total fuel consumption by fuel type and keep the fuel code with the
    # largest fuel consumption
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
    ].cat.rename_categories(
        {
            "left_only": "eia_only",
            "right_only": "cems_only",
            "both": "cems_and_eia",
        }
    )

    # assign a BA code to each plant
    plant_attributes = assign_ba_code_to_plant(plant_attributes, year)

    # add a flag about whether the plant is distribution connected
    plant_attributes = identify_distribution_connected_plants(
        plant_attributes, year, voltage_threshold_kv=60
    )

    # assign a fuel category to each plant based on what is most likely to match with
    # the category used in EIA-930
    plant_attributes = assign_fuel_category_to_ESC(
        df=plant_attributes,
        esc_column="plant_primary_fuel",
    )

    # add geographical info
    plant_attributes = add_plant_entity(plant_attributes)

    # add nameplate capacity
    plant_attributes = add_plant_nameplate_capacity(year, plant_attributes)

    # add operating and retirement dates
    plant_attributes = add_plant_operating_and_retirement_dates(plant_attributes)

    # convert types
    plant_attributes = apply_dtypes(plant_attributes)

    # change order of columns
    new_column_ordering = [
        "plant_id_eia",
        "plant_name_eia",
        "capacity_mw",
        "plant_primary_fuel",
        "fuel_category",
        "fuel_category_eia930",
        "state",
        "county",
        "city",
        "ba_code",
        "ba_code_physical",
        "latitude",
        "longitude",
        "plant_operating_date",
        "plant_retirement_date",
        "distribution_flag",
        "timezone",
        "data_availability",
    ]
    plant_attributes = plant_attributes[new_column_ordering]

    return plant_attributes


def assign_ba_code_to_plant(df: pd.DataFrame, year: int) -> pd.DataFrame:
    """Assigns a balancing authority code to each plant based on the plant id.

    Args:
         df (pd.DataFrame): data frame containing a 'plant_id_eia' column.
         year (int): a four-digit year.

     Returns:
         pd.DataFrame: original data frame with additional 'ba_code' column.
    """
    plant_ba = create_plant_ba_table(year)[
        ["plant_id_eia", "ba_code", "ba_code_physical"]
    ]

    # merge the ba code into the dataframe
    df = df.merge(plant_ba, how="left", on="plant_id_eia", validate="m:1")

    if len(df[df["ba_code"].isna()]) > 0:
        logger.warning("the following plants are missing ba_code:")
        logger.warning("\n" + df[df["ba_code"].isna()].to_string())

    # replace missing ba codes with NA
    df["ba_code"] = df["ba_code"].fillna("NA")
    df["ba_code_physical"] = df["ba_code_physical"].fillna("NA")

    return df


def create_plant_ba_table(year: int) -> pd.DataFrame:
    """Creates a table assigning a BA code and physical BA code to each plant ID.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table relating BA codes to plant ID.
    """

    plant_ba = load_data.load_pudl_table(
        "out_eia__yearly_plants",
        columns=[
            "plant_id_eia",
            "report_date",
            "balancing_authority_code_eia",
            "balancing_authority_name_eia",
            "utility_id_eia",
            "transmission_distribution_owner_name",
        ],
    )

    # for some earlier years, the plants data is missing BA codes.
    # backfill and forwardfill to make sure that we have complete data for all years, if
    # data is available for any year
    for col in [
        "balancing_authority_code_eia",
        "utility_id_eia",
        "balancing_authority_name_eia",
        "transmission_distribution_owner_name",
    ]:
        plant_ba[col] = plant_ba.groupby(["plant_id_eia"])[col].bfill()
        plant_ba[col] = plant_ba.groupby(["plant_id_eia"])[col].ffill()

    # some plants only have a record for years after the current year. To help ensure
    # that we have complete BA codes, create a dataframe containing only those plants
    # whose first record is after the current year, so that we can add these plants back
    # to plant_ba after filtering
    plant_ba_only_data_after_year = plant_ba[
        plant_ba.groupby(["plant_id_eia"])["report_date"].transform("min").dt.year
        > year
    ]
    # only keep the oldest record
    plant_ba_only_data_after_year = plant_ba_only_data_after_year[
        plant_ba_only_data_after_year["report_date"]
        == plant_ba_only_data_after_year.groupby(["plant_id_eia"])[
            "report_date"
        ].transform("min")
    ]

    # remove report dates newer than the current year
    plant_ba = plant_ba[plant_ba["report_date"].dt.year <= year]

    # sort the data from newest to oldest
    plant_ba = plant_ba.sort_values(by=["plant_id_eia", "report_date"], ascending=False)

    # add back plants that only have records after the current year
    # if for some reason this adds a duplicate plant, this will be dropped in the next
    # step since these records will be added to the end of the dataframe
    plant_ba = pd.concat([plant_ba, plant_ba_only_data_after_year], axis=0)

    # only keep the most recent row of data
    plant_ba = plant_ba.drop_duplicates(subset=["plant_id_eia"], keep="first")

    # merge utility name
    utilities_eia = load_data.load_pudl_table(
        "core_eia__entity_utilities", columns=["utility_id_eia", "utility_name_eia"]
    )
    plant_ba = plant_ba.merge(
        utilities_eia, how="left", on="utility_id_eia", validate="m:1"
    )
    # merge plant state
    plant_states = load_data.load_pudl_table(
        "core_eia__entity_plants", columns=["plant_id_eia", "state"]
    )
    plant_ba = plant_ba.merge(
        plant_states, how="left", on="plant_id_eia", validate="m:1"
    )

    # load the ba name reference
    ba_name_to_ba_code = pd.read_csv(
        reference_table_folder("ba_reference.csv"),
        dtype={"ba_name": "string", "ba_code": "string"},
    )
    ba_name_to_ba_code = dict(
        zip(
            ba_name_to_ba_code["ba_name"],
            ba_name_to_ba_code["ba_code"],
        )
    )

    # specify a ba code for certain utilities
    utility_as_ba_code = pd.read_csv(
        reference_table_folder("utility_name_ba_code_map.csv"),
        dtype={"name": "string", "ba_code": "string"},
    )
    utility_as_ba_code = dict(
        zip(
            utility_as_ba_code["name"],
            utility_as_ba_code["ba_code"],
        )
    )

    # fill missing BA codes first based on the BA name, then utility name, then on
    # the transmisison owner name
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["balancing_authority_name_eia"].map(ba_name_to_ba_code))
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
    # if there are any plants that have been assigned to a retired BA, set its BA code
    # as missing
    plant_ba.loc[plant_ba["ba_code"].isin(retired_bas), "ba_code"] = pd.NA

    # for plants without a BA code assign the miscellaneous BA code based on the state
    plant_ba["ba_code"] = plant_ba["ba_code"].fillna(plant_ba["state"] + "MS")

    # add a physical ba code based on the owner of the transmission system
    plant_ba["ba_code_physical"] = plant_ba["ba_code"]
    plant_ba.update(
        {
            "ba_code_physical": plant_ba["transmission_distribution_owner_name"].map(
                utility_as_ba_code
            )
        }
    )

    # update based on mapping table when ambiguous
    physical_ba = pd.read_csv(
        reference_table_folder("physical_ba.csv"),
        dtype=get_dtypes(),
    )
    plant_ba = plant_ba.merge(
        physical_ba,
        how="left",
        on=["ba_code", "transmission_distribution_owner_name"],
        suffixes=("", "_map"),
        validate="m:1",
    )
    plant_ba.update({"ba_code_physical": plant_ba["ba_code_physical_map"]})
    plant_ba.drop(columns="state")

    return plant_ba


def add_plant_operating_and_retirement_dates(df: pd.DataFrame) -> pd.DataFrame:
    """Adds the operating and retirement dates of a plant to input data frame. The
    operating date of a plant is taken as the earliest date among all generators'
    operating date over all report dates. Likewise, the retirement date of a plant is
    taken as the latest date among all generators' retirement date over all report
    dates.

    Note that the operating date is the date the generator began commercial operation.
    The retirement date is the date of the scheduled or effected retirement of the
    generator.

    Args:
        df (pd.DataFrame): table with a 'plant_id_eia' column.

    Returns:
        pd.DataFrame: original data frame with additional 'plant_operating_date' and
            'plant_retirement_date' column.
    """
    generator_dates = load_data.load_pudl_table(
        "out_eia__yearly_generators",
        year=earliest_data_year,
        end_year=latest_validated_year,
        columns=[
            "plant_id_eia",
            "generator_id",
            "report_date",
            "generator_operating_date",
            "generator_retirement_date",
        ],
    ).sort_values(by=["plant_id_eia", "generator_id", "report_date"], ascending=True)

    # fill missing dates
    date_columns = ["generator_operating_date", "generator_retirement_date"]

    for col in date_columns:
        generator_dates[col] = generator_dates.groupby(
            ["plant_id_eia", "generator_id"]
        )[col].bfill()
        generator_dates[col] = generator_dates.groupby(
            ["plant_id_eia", "generator_id"]
        )[col].ffill()

    # keep only the most recent year of data
    generator_dates = generator_dates.drop_duplicates(
        subset=["plant_id_eia", "generator_id"], keep="last"
    )

    plant_dates = (
        generator_dates.groupby("plant_id_eia")[
            ["generator_operating_date", "generator_retirement_date"]
        ]
        .agg(
            {
                "generator_operating_date": "min",
                "generator_retirement_date": lambda x: x.max(skipna=False),
            }
        )
        .rename(
            columns={
                "generator_operating_date": "plant_operating_date",
                "generator_retirement_date": "plant_retirement_date",
            }
        )
    )

    df = df.merge(plant_dates, how="left", on=["plant_id_eia"], validate="1:1")

    return df


def add_plant_nameplate_capacity(year: int, df: pd.DataFrame) -> pd.DataFrame:
    """Adds nameplate capacity to input data frame.

    Args:
        year (int): a four-digit year.
        df (pd.DataFrame): table with a 'plant_id_eia' column.

    Returns:
        pd.DataFrame: original data frame with additional 'capacity_mw' column.
    """
    generator_capacity = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year=earliest_data_year,
        end_year=latest_validated_year,
        columns=["plant_id_eia", "generator_id", "report_date", "capacity_mw"],
    ).sort_values(by=["plant_id_eia", "generator_id", "report_date"], ascending=True)

    generator_capacity["capacity_mw"] = generator_capacity.groupby(
        ["plant_id_eia", "generator_id"]
    )["capacity_mw"].bfill()
    generator_capacity["capacity_mw"] = generator_capacity.groupby(
        ["plant_id_eia", "generator_id"]
    )["capacity_mw"].ffill()

    # keep only the specified year of data
    generator_capacity = generator_capacity[
        generator_capacity["report_date"].dt.year == year
    ]

    plant_capacity = (
        generator_capacity.groupby(["plant_id_eia"])["capacity_mw"]
        .sum()
        .round(2)
        .reset_index()
    )

    df = df.merge(plant_capacity, how="left", on=["plant_id_eia"], validate="1:1")

    return df


def identify_distribution_connected_plants(
    df: pd.DataFrame, year: int, voltage_threshold_kv: int = 60
) -> pd.DataFrame:
    """Identifies which plant are "distribution grid connected" based on a voltage
    threshold.

    Args:
        df (pd.DataFrame): table with a 'plant_id_eia' column.
        year (int): a four-digit year.
        voltage_threshold_kv (int, optional): the voltage (kV) under which a plant
            will be considered to be a distribution asset. Defaults to 60.

    Returns:
        pd.DataFrame: original column with additional 'distribution_flag' binary column.
    """

    # load the EIA-860 data
    plant_voltage = load_data.load_pudl_table(
        "out_eia__yearly_plants", year, columns=["plant_id_eia", "grid_voltage_1_kv"]
    )

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
    df: pd.DataFrame,
    fuel_category_names: list = ["fuel_category", "fuel_category_eia930"],
    esc_column: str = "energy_source_code",
) -> pd.DataFrame:
    """Assigns a fuel category to each energy source code in a dataframe.

    Args:
        df (pd.DataFrame): table with column name that matches `fuel_category_names`
            and 'energy_source_codes
        fuel_category_names (list, optional): columns in
            reference_tables/energy_source_groups.csv that contains the desired
            category mapping. Defaults to ["fuel_category", "fuel_category_eia930"].
        esc_column (str, optional): name of the column in `df` that contains the energy
            source codes to assign a category to. Defaults to "energy_source_code".

    Returns:
        pd.DataFrame: original data frame with additional 'fuel_category_names'
            column(s).
    """
    # load the fuel category table
    energy_source_groups = pd.read_csv(
        reference_table_folder("energy_source_groups.csv"),
        dtype=get_dtypes(),
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


def add_plant_entity(df: pd.DataFrame) -> pd.DataFrame:
    """Adds timezone and geographical information to input data frame.

    Args:
        df (pd.DataFrame): table with a 'plant_id_eia' column.

    Returns:
        pd.DataFrame: original data frame with additional 'timezone', 'lat', 'lon',
            'state', 'county', 'city', 'zip_code' and 'street_address' columns
    """
    eia860_info = [
        "latitude",
        "longitude",
        "state",
        "county",
        "city",
        "plant_name_eia",
    ]
    plants_entity = load_data.load_pudl_table(
        "core_eia__entity_plants",
        columns=["plant_id_eia", "timezone"] + eia860_info,
    )
    plants_entity_from_eia860 = load_data.load_raw_eia860_plant_geographical_info(
        latest_validated_year
    )
    complete_plants_entity = plants_entity.merge(
        plants_entity_from_eia860,
        how="left",
        on=["plant_id_eia"],
        validate="1:1",
        suffixes=["", "_eia"],
    )

    for c in eia860_info:
        if complete_plants_entity[c].isna().sum() > 0:
            complete_plants_entity[c] = complete_plants_entity[c].fillna(
                complete_plants_entity[f"{c}_eia"]
            )
        complete_plants_entity = complete_plants_entity.drop(columns=f"{c}_eia")

    df = df.merge(
        complete_plants_entity, how="left", on=["plant_id_eia"], validate="m:1"
    )

    return df


def add_subplant_ids_to_df(
    df: pd.DataFrame,
    year: int,
    plant_part_to_map: str,
    how_merge: str,
    validate_merge: str,
) -> pd.DataFrame:
    """Helper function for adding subplant_id's to a dataframe based on the `plant_part_to_map`.

    Validates that all rows have non-missing subplant_id's assigned.

    Args:
        df (pd.DataFrame): Dataframe to add subplant_id's to
        year (int): the data year
        plant_part_to_map (str): one of the plant parts in the subplant crosswalk that
            exists in df to use to map the subplant ids. Acceptable options include:
            "generator_id", "boiler_id", or "emissions_unit_id_epa"
        how_merge (str): used for `how=` argument of pd.merge(): e.g. "left", "right",
            "inner", "outer"
        validate_merge (str): used for `validate=` argument of pd.merge(): eg. "1:1",
            "1:m","m:1","m:m". Depends on the merge key specified in plant_part_to_map
            and the structure of `df`

    Returns:
        pd.DataFrame: df with subplant_id column added
    """
    # add subplant id, dropping duplicate entries for the relevant plant part
    subplant_crosswalk = (
        pd.read_csv(
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv.zip"),
            dtype=get_dtypes(),
        )[["plant_id_eia", plant_part_to_map, "subplant_id"]]
        .drop_duplicates()
        .dropna(subset=plant_part_to_map)
    )
    df = df.merge(
        subplant_crosswalk,
        how=how_merge,
        on=["plant_id_eia", plant_part_to_map],
        validate=validate_merge,
    )
    validation.test_for_missing_subplant_id(df, plant_part_to_map)

    return df
