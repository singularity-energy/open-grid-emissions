import numpy as np
import pandas as pd

from oge.column_checks import get_dtypes, apply_dtypes
from oge.filepaths import reference_table_folder, outputs_folder
import oge.load_data as load_data
from oge.logging_util import get_logger
import oge.validation as validation

logger = get_logger(__name__)


def create_plant_attributes_table(cems, eia923_allocated, year, primary_fuel_table):
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
    ].cat.rename_categories(
        {
            "left_only": "eia_only",
            "right_only": "cems_only",
            "both": "cems_and_eia",
        }
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
        logger.warning("the following plants are missing ba_code:")
        logger.warning("\n" + df[df["ba_code"].isna()].to_string())

    # replace missing ba codes with NA
    df["ba_code"] = df["ba_code"].fillna("NA")
    df["ba_code_physical"] = df["ba_code_physical"].fillna("NA")

    return df


def create_plant_ba_table(year):
    """
    Creates a table assigning a ba_code and physical ba code to each plant id.
    """

    plant_ba = load_data.load_pudl_table(
        "plants_eia860",
        columns=[
            "plant_id_eia",
            "report_date",
            "balancing_authority_code_eia",
            "balancing_authority_name_eia",
            "utility_id_eia",
            "transmission_distribution_owner_name",
        ],
    )

    # remove report dates newer than the current year
    plant_ba = plant_ba[plant_ba["report_date"].dt.year <= year]

    # sort the data from newest to oldest
    plant_ba = plant_ba.sort_values(by=["plant_id_eia", "report_date"], ascending=False)

    # only keep the most recent row of data
    plant_ba = plant_ba.drop_duplicates(subset=["plant_id_eia"], keep="first")

    # merge utility name
    utilities_eia = load_data.load_pudl_table(
        "utilities_eia", columns=["utility_id_eia", "utility_name_eia"]
    )
    plant_ba = plant_ba.merge(
        utilities_eia, how="left", on="utility_id_eia", validate="m:1"
    )
    # merge plant state
    plant_states = load_data.load_pudl_table(
        "plants_entity_eia", columns=["plant_id_eia", "state"]
    )
    plant_ba = plant_ba.merge(
        plant_states, how="left", on="plant_id_eia", validate="m:1"
    )

    # convert the dtype of the balancing authority code column from string to object
    # this will allow for missing values to be filled
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].astype(object)
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(value=np.NaN)

    # load the ba name reference
    ba_name_to_ba_code = pd.read_csv(reference_table_folder("ba_reference.csv"))
    ba_name_to_ba_code = dict(
        zip(
            ba_name_to_ba_code["ba_name"],
            ba_name_to_ba_code["ba_code"],
        )
    )

    # specify a ba code for certain utilities
    utility_as_ba_code = pd.read_csv(
        reference_table_folder("utility_name_ba_code_map.csv")
    )
    utility_as_ba_code = dict(
        zip(
            utility_as_ba_code["name"],
            utility_as_ba_code["ba_code"],
        )
    )

    # fill missing BA codes first based on the BA name, then utility name, then on the transmisison owner name
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
    # if there are any plants that have been assigned to a retired BA, set its BA code as missing
    plant_ba.loc[plant_ba["ba_code"].isin(retired_bas), "ba_code"] = np.NaN

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
        reference_table_folder("physical_ba.csv"), dtype=get_dtypes()
    )
    plant_ba = plant_ba.merge(
        physical_ba,
        how="left",
        on=["ba_code", "transmission_distribution_owner_name"],
        suffixes=("", "_map"),
        validate="m:1",
    )
    plant_ba.update({"ba_code_physical": plant_ba["ba_code_physical_map"]})

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
    plant_voltage = load_data.load_pudl_table(
        "plants_eia860", year, columns=["plant_id_eia", "grid_voltage_1_kv"]
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
        reference_table_folder("energy_source_groups.csv"), dtype=get_dtypes()
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
    plant_tz = load_data.load_pudl_table(
        "plants_entity_eia", columns=["plant_id_eia", "timezone"]
    )
    df = df.merge(plant_tz, how="left", on=["plant_id_eia"], validate="m:1")

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
            outputs_folder(f"{year}/subplant_crosswalk_{year}.csv"),
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
