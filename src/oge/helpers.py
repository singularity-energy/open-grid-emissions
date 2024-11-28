import numpy as np
import pandas as pd

from geopy.geocoders import Nominatim
from geopy.exc import GeocoderUnavailable
from timezonefinder import TimezoneFinder
from urllib3.exceptions import ReadTimeoutError

from oge.column_checks import get_dtypes, apply_dtypes, DATA_COLUMNS
from oge.constants import (
    earliest_data_year,
    latest_validated_year,
    current_early_release_year,
)
from oge.filepaths import reference_table_folder, outputs_folder

import oge.load_data as load_data
from oge.logging_util import get_logger
import oge.validation as validation

logger = get_logger(__name__)

tf = TimezoneFinder()
geolocator = Nominatim(user_agent="oge")


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
    plant_attributes = assign_fuel_category_to_esc(
        df=plant_attributes,
        esc_column="plant_primary_fuel",
    )

    # add geographical info
    plant_attributes = add_plant_entity(plant_attributes)

    # fill out missing location/coordinates
    plant_attributes = add_missing_location(plant_attributes)

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

    # test for missing values
    validation.test_for_missing_values(
        plant_attributes, skip_cols=["plant_retirement_date"]
    )

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


def assign_fleet_to_subplant_data(
    subplant_data: pd.DataFrame,
    plant_attributes_table: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    ba_col: str = "ba_code",
    primary_fuel_col: str = "subplant_primary_fuel",
    fuel_category_col: str = "fuel_category",
    other_attribute_cols: list[str] = [],
    drop_primary_fuel_col: bool = True,
) -> pd.DataFrame:
    """Assigns a BA code and fuel category to each subplant in order to facilitate
    aggregating the data to the fleet level.

    When assigning a primary fuel/fuel category, the general options we should follow
    are:
        - For applying to CEMS data for the residual hourly profile calculation:
            - primary_fuel_col = "subplant_primary_fuel_from_capacity_mw"
            - fuel_category_col = "fuel_category_eia930"
            - Notes: this matches how generators would likely be classified for 930
            reporting
        - For applying to monthly EIA-923 data to shape:
            - primary_fuel_col = "subplant_primary_fuel_from_capacity_mw"
            - fuel_category_col = "fuel_category"
            - Notes: This means that we are assigning our "flat" profiles for all the
            fuels that are categorized as "other" in 930, rather than using the "other"
            profile for all of these individual categories
        - For aggregating subplant data to fleet-level results:
            - primary_fuel_col = "subplant_primary_fuel"
            - fuel_category_col = "fuel_category"


    Args:
        subplant_data (pd.DataFrame): dataframe to assign fleet to
        plant_attributes_table (pd.DataFrame): static plant attributes
        primary_fuel_table (pd.DataFrame): table of subplant-level primary fuels
        ba_col (str, optional): Whether to use commercial "ba_code" balancing area
            definition or physical BA "ba_code_physical". Defaults to "ba_code".
        primary_fuel_col (str, optional): Name of column from primary_fuel_table to use
            to assign a fuel type to the subplant. Defaults to "subplant_primary_fuel".
        fuel_category_col (str, optional): name of fuel category column to map to the
            energy source code specified by the primary_fuel_col in primary_fuel_table.
            Defaults to "fuel_category".
        other_attribute_cols (list[str], optional): a list of additional columns from
            plant_attributes_table to add to subplant_data. Defaults to [].
        drop_primary_fuel_col (bool): Whether to drop the ESC-level primary_fuel_col
            before returning the table. Can be set to False for use of this function
            in validation.identify_percent_of_data_by_input_source() Defaults to True.

    Raises:
        UserWarning: If a BA code or fuel type cannot be assigned to a subplant

    Returns:
        pd.DataFrame: subplant_data with ba_code and fuel_category columns added
    """

    # check to make sure the ba_col and primary_fuel_col are not already in the dataframe
    # if so, drop them before merging
    cols_to_add = [ba_col, primary_fuel_col] + other_attribute_cols
    fleet_cols_already_in_subplant_data = [
        col for col in subplant_data.columns if col in cols_to_add
    ]
    if len(fleet_cols_already_in_subplant_data) > 0:
        subplant_data = subplant_data.drop(columns=fleet_cols_already_in_subplant_data)

    # Assign a BA to the data
    subplant_data = subplant_data.merge(
        plant_attributes_table[["plant_id_eia", ba_col] + other_attribute_cols].rename(
            columns={ba_col: "ba_code"}
        ),
        how="left",
        on=["plant_id_eia"],
        validate="m:1",
    )

    # assign a fuel category to the subplant
    # ensure no missing primary fuel
    if "subplant" in primary_fuel_col:
        default_col = "subplant_primary_fuel"
    else:
        default_col = "plant_primary_fuel"
    primary_fuel_table[primary_fuel_col] = primary_fuel_table[primary_fuel_col].fillna(
        primary_fuel_table[default_col]
    )
    subplant_primary_fuel = primary_fuel_table[
        ["plant_id_eia", "subplant_id", primary_fuel_col]
    ].drop_duplicates()
    subplant_primary_fuel = assign_fuel_category_to_esc(
        subplant_primary_fuel,
        fuel_category_names=[fuel_category_col],
        esc_column=primary_fuel_col,
    )
    if drop_primary_fuel_col:
        subplant_primary_fuel = subplant_primary_fuel.drop(columns=[primary_fuel_col])
    # merge in the fuel data
    subplant_data = subplant_data.merge(
        subplant_primary_fuel,
        how="left",
        on=["plant_id_eia", "subplant_id"],
        validate="m:1",
    )

    # check that there is no missing ba or fuel codes for subplants with nonzero gen
    # for CEMS data, check only units that report positive gross generaiton
    if "gross_generation_mwh" in subplant_data.columns:
        missing_fleet_keys = subplant_data[
            (
                (subplant_data["ba_code"].isna())
                | (subplant_data[fuel_category_col].isna())
                & (subplant_data["gross_generation_mwh"] > 0)
            )
        ]
    # otherwise, check units that report non-zero net generation
    else:
        missing_fleet_keys = subplant_data[
            (
                (subplant_data["ba_code"].isna())
                | (subplant_data[fuel_category_col].isna())
                & (subplant_data["net_generation_mwh"] != 0)
            )
        ]
    if len(missing_fleet_keys) > 0:
        logger.warning(
            missing_fleet_keys.groupby(
                [
                    "plant_id_eia",
                    "subplant_id",
                    "ba_code",
                    fuel_category_col,
                ],
                dropna=False,
            )[
                [
                    "net_generation_mwh",
                ]
            ]
            .sum()
            .to_string()
        )
        raise UserWarning(
            "The plant attributes table is missing ba_code or fuel_category data for some plants. This will result in incomplete power sector results."
        )

    return subplant_data


def combine_subplant_data(
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
    eia_data: pd.DataFrame,
    resolution: str,
    validate: bool = True,
) -> pd.DataFrame:
    """Combines subplant-level data from multiple sources into a single dataframe.
    Data can be returned at either the monthly or hourly resolution.

    When passing hourly data in later in the pipeline, we only want to combine cems data
    so we can pass an empty dataframe for eia_data.

    Args:
        cems (pd.DataFrame): One of the dfs to combine
        partial_cems_subplant (pd.DataFrame): One of the dfs to combine
        partial_cems_plant (pd.DataFrame): One of the dfs to combine
        eia_data (pd.DataFrame): One of the dfs to combine
        resolution (str): Whether to combine "hourly" data or "monthly" data. All input
            dataframes should have "datetime_utc" columns for the former, and
            "report_date" columns for the latter
        validate (bool, optional): Whether to ensure non-overlapping data from all
            sources. Sometimes not necessary based on whether this has already been
            checked. Defaults to True.

    Raises:
        UserWarning: If acceptable option for resolution arg not passed

    Returns:
        pd.DataFrame: combined data from the four input dfs at the resolution
    """

    KEY_COLUMNS = [
        "plant_id_eia",
        "subplant_id",
        "report_date",
    ]
    if resolution == "hourly":
        KEY_COLUMNS += ["datetime_utc"]
    elif resolution == "monthly":
        pass
    else:
        raise UserWarning(
            f"`resolution` must be 'monthly' or 'hourly'. '{resolution}' specified"
        )

    ALL_COLUMNS = KEY_COLUMNS + DATA_COLUMNS

    if validate:
        validation.ensure_non_overlapping_data_from_all_sources(
            cems, partial_cems_subplant, partial_cems_plant, eia_data
        )

    # group data by subplant-month or subplant-hour and filter columns
    cems = (
        cems.groupby(
            KEY_COLUMNS,
            dropna=False,
        )
        .sum(numeric_only=True)
        .reset_index()[[col for col in cems.columns if col in ALL_COLUMNS]]
    )
    # don't group if there is no data in the dataframe
    if len(partial_cems_subplant) > 0:
        partial_cems_subplant = (
            partial_cems_subplant.groupby(
                KEY_COLUMNS,
                dropna=False,
            )
            .sum(numeric_only=True)
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
            .sum(numeric_only=True)
            .reset_index()[
                [col for col in partial_cems_plant.columns if col in ALL_COLUMNS]
            ]
        )
    eia_data = (
        eia_data.groupby(
            KEY_COLUMNS,
            dropna=False,
        )
        .sum(numeric_only=True)
        .reset_index()[[col for col in eia_data.columns if col in ALL_COLUMNS]]
    )

    # concat together
    combined_subplant_data = pd.concat(
        [cems, partial_cems_subplant, partial_cems_plant, eia_data],
        axis=0,
        ignore_index=True,
        copy=False,
    )

    # groupby subplant after combining in case subplant reported multiple places
    combined_subplant_data = (
        combined_subplant_data.groupby(KEY_COLUMNS, dropna=False)[DATA_COLUMNS]
        .sum(numeric_only=True)
        .reset_index()
    )

    # re-order the columns
    combined_subplant_data = combined_subplant_data[ALL_COLUMNS]

    return combined_subplant_data


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
        plant_states, how="outer", on="plant_id_eia", validate="m:1"
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

    plant_ba["ba_code"] = plant_ba["ba_code"].replace("None", pd.NA)

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
        end_year=max(latest_validated_year, current_early_release_year),
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

    Includes multiple steps for ensuring complete data for capacity. These values will
    only be used if there is no "existing" capcity reported for a year.
    1. Some plants report data before reporting to EIA-860. In this case, find the
        earliest reported capacity for the plant and fill missing capacity values.
    2. Sometimes plants report data before they are operational. Find the capacity of
        any proposed generators and use this to fill missing capacity values.
    3. Sometimes plants report data after they have retired. In this case, find the
        capacity of all generators in the latest year in which the plant was operating.

    Args:
        year (int): a four-digit year.
        df (pd.DataFrame): table with a 'plant_id_eia' column.

    Returns:
        pd.DataFrame: original data frame with additional 'capacity_mw' column.
    """
    generator_capacity = load_data.load_pudl_table(
        "core_eia860__scd_generators",
        year=earliest_data_year,
        end_year=max(latest_validated_year, current_early_release_year),
        columns=[
            "plant_id_eia",
            "generator_id",
            "report_date",
            "capacity_mw",
            "operational_status",
            "generator_retirement_date",
        ],
    ).sort_values(by=["plant_id_eia", "generator_id", "report_date"], ascending=True)

    generator_capacity["capacity_mw"] = generator_capacity.groupby(
        ["plant_id_eia", "generator_id"]
    )["capacity_mw"].bfill()
    generator_capacity["capacity_mw"] = generator_capacity.groupby(
        ["plant_id_eia", "generator_id"]
    )["capacity_mw"].ffill()

    # fill missing operational status with the last known status
    generator_capacity["operational_status"] = generator_capacity.groupby(
        ["plant_id_eia", "generator_id"]
    )["operational_status"].ffill()

    # if a generator retires in a year, change the status to existing for that year
    generator_capacity.loc[
        generator_capacity["generator_retirement_date"]
        >= generator_capacity["report_date"],
        "operational_status",
    ] = "existing"

    # find the earliest reported plant capacity - used for filling missing values later
    # only keep the earliest year of reported data
    earliest_plant_capacity = generator_capacity.drop_duplicates(
        subset=["plant_id_eia", "generator_id"], keep="first"
    )
    # only keep capacities for generators reported in the earliest year
    earliest_plant_capacity = earliest_plant_capacity[
        earliest_plant_capacity["report_date"]
        == earliest_plant_capacity.groupby(["plant_id_eia", "generator_id"])[
            "report_date"
        ].transform("min")
    ].reset_index()
    # calculate the plant capacity
    earliest_plant_capacity = (
        earliest_plant_capacity.groupby("plant_id_eia")["capacity_mw"]
        .sum()
        .reset_index()
    )

    # remove any years after the current year
    generator_capacity = generator_capacity[
        generator_capacity["report_date"].dt.year <= year
    ]

    # get all proposed generators in the current year
    proposed_plant_capacity = (
        generator_capacity[
            (generator_capacity["operational_status"] == "proposed")
            & (generator_capacity["report_date"].dt.year == year)
        ]
        .groupby(["plant_id_eia"])["capacity_mw"]
        .sum()
        .reset_index()
    )

    # Only consider generators that are existing (operating, standby, etc.)
    generator_capacity = generator_capacity[
        generator_capacity["operational_status"] == "existing"
    ]

    # get the total plant capacity in the year when the plant was most recently operating
    latest_plant_capacity = (
        generator_capacity[
            generator_capacity["report_date"]
            == generator_capacity.groupby(["plant_id_eia", "generator_id"])[
                "report_date"
            ].transform("max")
        ]
        .groupby(["plant_id_eia"])["capacity_mw"]
        .sum()
        .reset_index()
    )

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

    # if there are any plants that are missing from plant_capacity because they've
    # already retired, add the most recent capacity for that plant
    plant_capacity = pd.concat(
        [
            plant_capacity,
            latest_plant_capacity[
                ~latest_plant_capacity["plant_id_eia"].isin(
                    list(plant_capacity.plant_id_eia.unique())
                )
            ],
            proposed_plant_capacity[
                ~proposed_plant_capacity["plant_id_eia"].isin(
                    list(plant_capacity.plant_id_eia.unique())
                )
                & ~proposed_plant_capacity["plant_id_eia"].isin(
                    list(latest_plant_capacity.plant_id_eia.unique())
                )
            ],
            earliest_plant_capacity[
                ~earliest_plant_capacity["plant_id_eia"].isin(
                    list(plant_capacity.plant_id_eia.unique())
                )
                & ~earliest_plant_capacity["plant_id_eia"].isin(
                    list(latest_plant_capacity.plant_id_eia.unique())
                )
                & ~earliest_plant_capacity["plant_id_eia"].isin(
                    list(proposed_plant_capacity.plant_id_eia.unique())
                )
            ],
        ],
        axis=0,
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


def assign_fuel_category_to_esc(
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
        max(latest_validated_year, current_early_release_year)
    )
    complete_plants_entity = plants_entity.merge(
        plants_entity_from_eia860,
        how="left",
        on=["plant_id_eia"],
        validate="1:1",
        suffixes=["", "_eia"],
    )

    for c in eia860_info:
        # Handle NAs
        if complete_plants_entity[c].isna().sum() > 0:
            complete_plants_entity[c] = complete_plants_entity[c].fillna(
                complete_plants_entity[f"{c}_eia"]
            )
        # Handle positive longitude
        if c == "longitude" and (complete_plants_entity[c] > 0).any():
            # Replace if EIA-860 longitude is negative, otherwise flip the sign.
            for i in complete_plants_entity[complete_plants_entity[c] > 0].index:
                lat_eia = complete_plants_entity.loc[i, "latitude_eia"]
                lon_eia = complete_plants_entity.loc[i, "longitude_eia"]
                if lon_eia < 0:
                    complete_plants_entity.loc[i, "latitude"] = lat_eia
                    complete_plants_entity.loc[i, "longitude"] = lon_eia
                # Otherwise flip the sign of longitude and keep PUDL latitude
                else:
                    complete_plants_entity.loc[
                        i, "longitude"
                    ] = -complete_plants_entity.loc[i, "longitude"]
                # Get new timezone
                complete_plants_entity.loc[i, "timezone"] = tf.timezone_at(
                    lng=complete_plants_entity.loc[i, "longitude"],
                    lat=complete_plants_entity.loc[i, "latitude"],
                )

    # Clean data frame
    complete_plants_entity = complete_plants_entity.drop(
        columns=[f"{c}_eia" for c in eia860_info]
    )

    df = df.merge(
        complete_plants_entity, how="left", on=["plant_id_eia"], validate="m:1"
    )

    return df


def add_missing_location(df: pd.DataFrame) -> pd.DataFrame:
    """Add missing latitude, longitude, state, county and city when possible.

    Args:
        df (pd.DataFrame): table with 'latitude', 'longitude', 'state', 'county' and
            'city' columns

    Returns:
        pd.DataFrame: original data frame with missing 'latitude', 'longitude',
            'state', 'county' and 'city' filled out when possible.
    """
    # get lat/lon
    missing_coord = df[df["longitude"].isna() | df["latitude"].isna()]
    if len(missing_coord) > 0:
        # only get coordinates when state, county and city are available
        for i in missing_coord.index:
            state = df.loc[i, "state"]
            county = df.loc[i, "county"]
            city = df.loc[i, "city"]

            lat, lon = get_coordinates_of_location(state, county, city)
            df.loc[i, "latitude"] = lat
            df.loc[i, "longitude"] = lon

    # get missing state, county and city from coordinates
    missing_location = df[df["state"].isna() | df["county"].isna() | df["city"].isna()]
    if len(missing_location) > 0:
        for i in missing_location.index:
            if df.loc[i, ["latitude", "longitude"]].isna().sum() == 0:
                state, county, city = search_location_from_coordinates(
                    df.loc[i, "latitude"],
                    df.loc[i, "longitude"],
                )
                if pd.isna(df.loc[i, "state"]):
                    df.loc[i, "state"] = state
                if pd.isna(df.loc[i, "county"]):
                    df.loc[i, "county"] = county
                if pd.isna(df.loc[i, "city"]):
                    df.loc[i, "city"] = city

    return df


def search_location_from_coordinates(latitude: float, longitude: float) -> tuple[str]:
    """Get state, county, city at latitude/longitude.

    Example:
        >>> latitude = 33.458665
        >>> longitude = -87.35682
        >>> location = geolocator.reverse(f"{latitude}, {longitude}").raw
        >>> location
        {'place_id': 149439, 'licence': 'Data © OpenStreetMap contributors,
        ODbL 1.0. http://osm.org/copyright', 'osm_type': 'way', 'osm_id': 8885591,
        'lat': '33.460586', 'lon': '-87.359444', 'class': 'highway',
        'type': 'unclassified', 'place_rank': 26, 'importance': 0.10000999999999993,
        'addresstype': 'road', 'name': 'County Road 38', 'display_name':
        'County Road 38, Tuscaloosa County, Alabama, United States', 'address':
        {'road': 'County Road 38', 'county': 'Tuscaloosa County', 'state': 'Alabama',
        'ISO3166-2-lvl4': 'US-AL', 'country': 'United States', 'country_code': 'us'},
        'boundingbox': ['33.4552300', '33.4758370', '-87.3873120', '-87.3589650']}
    Args:
        latitude (float): latitude of the location.
        longitude (float): longitude of the location.

    Returns:
        tuple[str]: state, county and city of the location.
    """
    try:
        address = geolocator.reverse(f"{latitude}, {longitude}").raw["address"]
        if address["country_code"] != "us":
            return pd.NA, pd.NA, pd.NA
    except ReadTimeoutError:
        return pd.NA, pd.NA, pd.NA

    # Check for State
    state = (
        address["ISO3166-2-lvl4"].split("-")[1]
        if "ISO3166-2-lvl4" in address.keys()
        else pd.NA
    )
    county = address["county"].split(" ")[0] if "county" in address.keys() else pd.NA
    city = address["city"] if "city" in address.keys() else pd.NA
    return state, county, city


def get_coordinates_of_location(state: str, county: str, city: str) -> tuple[float]:
    """Use state, county and city information to get coordinates.

    Example:
            >>> location = geolocator.geocode("Bucks, Mobile county, AL").raw
            >>> location
            {'place_id': 379392554, 'licence': 'Data © OpenStreetMap contributors,
            ODbL 1.0. http://osm.org/copyright', 'osm_type': 'relation',
            'osm_id': 17019769, 'lat': '31.01630685', 'lon': '-88.02448016014876',
            'class': 'boundary', 'type': 'census', 'place_rank': 25,
            'importance': 0.4106648553911631, 'addresstype': 'census', 'name': 'Bucks',
            'display_name': 'Bucks, Mobile County, Alabama, United States',
            'boundingbox': ['31.0072629', '31.0244919', '-88.0286929', '-88.0198599']}

        Args:
            state (str): state of the location.
            county (str): county of the location.
            city (str): city of the location.


        Returns:
            tuple[float]: the latitude and longitude.
    """
    if pd.isna(state):
        return np.NaN, np.NaN
    if pd.isna(city) | (city == "unsited"):
        if not (pd.isna(county) | (county == "NOT IN FILE")):
            query = f"{county} county, {state}, USA"
        else:
            query = f"{state}, USA"
    else:
        query = f"{city}, {state}, USA"

    try:
        location = geolocator.geocode(query, country_codes="us")
        if location is None:
            logger.warning(f"No location returned for {query}")
            return np.NaN, np.NaN
        else:
            return float(location.raw["lat"]), float(location.raw["lon"])
    except ReadTimeoutError:
        logger.warning(f"ReadTimeoutError for {query}")
        return np.NaN, np.NaN
    except GeocoderUnavailable:
        logger.warning(f"GeocoderUnavailable for {query}")
        return np.NaN, np.NaN


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
