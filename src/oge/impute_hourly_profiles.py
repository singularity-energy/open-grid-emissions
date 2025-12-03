import gc
import pandas as pd
import numpy as np

# import open-grid-emissions modules
from oge.column_checks import apply_dtypes, DATA_COLUMNS
import oge.load_data as load_data
from oge.filepaths import reference_table_folder
import oge.validation as validation
import oge.output_data as output_data
from oge.logging_util import get_logger
from oge.helpers import assign_fleet_to_subplant_data

logger = get_logger(__name__)


def calculate_hourly_profiles(
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
    eia930_data: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    monthly_eia_data_to_shape: pd.DataFrame,
    year: int,
    transmission_only: bool = False,
    ba_column_name: str = "ba_code",
    use_flat: bool = False,
) -> pd.DataFrame:
    """Coordinating function for calculating the hourly residual profiles that will be
    used to shape monthly EIA data.

    First, we calculate residual profiles by subtracting known, hourly CEMS data from
    the hourly EIA-930, fleet-level profiles.

    Where these residual profiles cannot be calculated, we then estimate the profiles in
    impute_missing_hourly_profiles, either by using regional wind and solar profiles, or
    assuming flat generation for certain resources.

    Next, we add the average hourly profile of all subplants in a fleet that report CEMS
    data.

    Finally, of the various profile options created, we select the best available
    profile based on a hierarchy. The available profile types to select from include:
        - "residual_profile"
        - "shifted_residual_profile"
        - "eia930_profile"
        - "cems_profile"
        - "DIBA_average"
        - "national_average"
        - "assumed_flat"

    Args:
        cems (pd.DataFrame): Hourly, unit- or subplant-level CEMS data
        partial_cems_subplant (pd.DataFrame): Hourly data for partial-CEMS subplants
        partial_cems_plant (pd.DataFrame): Hourly data for partial-CEMS plants
        eia930_data (pd.DataFrame): Hourly EIA-930 net generation data for each fleet
        plant_attributes (pd.DataFrame): table of plant-level static attributes
        monthly_eia_data_to_shape (pd.DataFrame): Monthly, subplant-level EIA data for
            subplants that only report to EIA-923
        year (int): The data year
        transmission_only (bool, optional): Whether or not to only consider subplants
            that are connected to the transmission grid when calculating residual
            profiles. Defaults to False.
        ba_column_name (str, optional): Whether to use "ba_code" or "ba_code_physical"
            to define fleet composition. Defaults to "ba_code".
        use_flat (bool, optional): Whether or not to use a flat profile for all fleets.
            Defaults to False.

    Returns:
        pd.DataFrame: dataframe of hourly profiles (in MW) for each fleet. The table
            includes a profile for each method, as well as a selected best "profile".
            Profiles must still be converted to a percent of the monthly total using
            convert_profile_to_percent()
    """
    residual_profiles = calculate_residual(
        cems,
        partial_cems_subplant,
        partial_cems_plant,
        eia930_data,
        plant_attributes,
        primary_fuel_table,
        year,
        transmission_only=transmission_only,
        ba_column_name=ba_column_name,
    )

    # add backstop profiles where residuals were not available
    hourly_profiles = impute_missing_hourly_profiles(
        monthly_eia_data_to_shape,
        residual_profiles,
        plant_attributes,
        primary_fuel_table,
        year,
    )
    hourly_profiles = add_cems_backstop_profile(
        hourly_profiles, cems, plant_attributes, primary_fuel_table, year
    )

    # if there are any months that have incomplete cems data, replace the cems profile with na
    incomplete_cems = hourly_profiles.loc[
        hourly_profiles["cems_profile"].isna(),
        ["ba_code", "fuel_category", "report_date"],
    ].drop_duplicates()
    hourly_profiles = hourly_profiles.merge(
        incomplete_cems,
        how="outer",
        on=["ba_code", "fuel_category", "report_date"],
        indicator="source",
        validate="m:1",
    )
    hourly_profiles.loc[(hourly_profiles["source"] == "both"), "cems_profile"] = np.NaN
    hourly_profiles.drop(columns="source", inplace=True)

    hourly_profiles = select_best_available_profile(hourly_profiles)

    # round the data to the nearest tenth
    hourly_profiles["profile"] = hourly_profiles["profile"].round(1)

    # add a flat profile for negative generation
    hourly_profiles["flat_profile"] = 1.0

    # use flat profile?
    if use_flat:
        hourly_profiles["profile"] = hourly_profiles["flat_profile"]
        hourly_profiles["profile_method"] = "flat_profile"

    logger.info(
        "Summary of methods used to estimate missing hourly profiles (count of ba-months):"
    )
    summary_table = (
        hourly_profiles[["ba_code", "fuel_category", "report_date", "profile_method"]]
        .drop_duplicates()
        .drop(columns=["ba_code"])
        .pivot_table(index="fuel_category", columns="profile_method", aggfunc="count")
        .fillna(0)
        .astype(int)
        .droplevel(level=0, axis=1)
    )

    # Format outputs
    profile_methods = [
        "residual_profile",
        "shifted_residual_profile",
        "eia930_profile",
        "cems_profile",
        "DIBA_average",
        "national_average",
        "assumed_flat",
    ]
    # Add default value of zero: let users know this method exists but was not used
    for pm in profile_methods:
        if pm not in summary_table.columns:
            summary_table[pm] = 0
    # re-order columns
    summary_table = summary_table.loc[
        :,
        profile_methods,
    ]
    logger.info("\n" + summary_table.to_string())

    return hourly_profiles


def select_best_available_profile(hourly_profiles: pd.DataFrame) -> pd.DataFrame:
    """Selects the best available hourly profile from the options available.
    The order of preference is:
        1. If the residual profile does not have a negative total for a month, use that
        2. If the eia930 profile doesn't have missing data, use that next
        3. If there are at least 3 CEMS plants available, use the combined CEMS profile
        4. Use the imputed profile

    We could create two different profiles - one for positive and one for negative values

    Args:
        hourly_profiles (pd.DataFrame): dataframe of hourly profiles

    Returns:
        pd.DataFrame: hourly_profiles with a new "profile" column added
    """

    # create a filtered version of the residual profile, removing months where the 
    # residual contains negative values
    residual_filter = (
        hourly_profiles.groupby(["ba_code", "fuel_category", "report_date"])[
            "residual_profile"
        ]
        .min()
        .reset_index()
    )
    residual_filter = residual_filter[residual_filter["residual_profile"] < 0]

    hourly_profiles = hourly_profiles.merge(
        residual_filter[["ba_code", "fuel_category", "report_date"]],
        how="outer",
        on=["ba_code", "fuel_category", "report_date"],
        indicator="negative_filter",
        validate="m:1",
    )
    hourly_profiles["residual_profile_filtered"] = hourly_profiles["residual_profile"]
    hourly_profiles.loc[
        hourly_profiles["negative_filter"] == "both", "residual_profile_filtered"
    ] = np.NaN
    hourly_profiles.drop(columns=["negative_filter"], inplace=True)

    # implement a filter on the shifted residual profile so that we don't use it if 
    # greater than the eia930 data
    shifted_filter = (
        hourly_profiles.groupby(["ba_code", "fuel_category", "report_date"])
        .sum(numeric_only=True)
        .reset_index()
    )
    shifted_filter = shifted_filter.loc[
        shifted_filter["shifted_residual_profile"] > shifted_filter["eia930_profile"],
        ["ba_code", "fuel_category", "report_date"],
    ]

    hourly_profiles = hourly_profiles.merge(
        shifted_filter,
        how="outer",
        on=["ba_code", "fuel_category", "report_date"],
        indicator="shifted_filter",
        validate="m:1",
    )
    # create a new filtered column, replacing filtered values with nan
    hourly_profiles["shifted_residual_profile_filtered"] = hourly_profiles[
        "shifted_residual_profile"
    ]
    hourly_profiles.loc[
        hourly_profiles["shifted_filter"] == "both", "shifted_residual_profile_filtered"
    ] = np.NaN
    hourly_profiles.drop(columns=["shifted_filter"], inplace=True)

    # pick the profile

    hourly_profiles["profile"] = np.NaN
    hourly_profiles["profile_method"] = pd.NA
    # specify the profile as the best available data
    profile_hierarchy = [
        "residual_profile_filtered",
        "shifted_residual_profile_filtered",
        "eia930_profile",
        "cems_profile",
        "imputed_profile",
    ]
    for profile in profile_hierarchy:
        # fill in the name of the method
        hourly_profiles.loc[
            hourly_profiles["profile"].isna() & ~hourly_profiles[profile].isna(),
            "profile_method",
        ] = profile
        # fill missing values with non-missing values from the filtered profile data
        hourly_profiles["profile"] = hourly_profiles["profile"].fillna(
            hourly_profiles[profile]
        )

    # for imputed profiles, identify the specific imputation method
    hourly_profiles.loc[
        hourly_profiles["profile_method"] == "imputed_profile", "profile_method"
    ] = hourly_profiles.loc[
        hourly_profiles["profile_method"] == "imputed_profile", "imputation_method"
    ]
    # replace the filtered method names with the regular method name
    hourly_profiles["profile_method"] = hourly_profiles["profile_method"].replace(
        {
            "shifted_residual_profile_filtered": "shifted_residual_profile",
            "residual_profile_filtered": "residual_profile",
        }
    )
    hourly_profiles.drop(
        columns=[
            "imputation_method",
            "residual_profile_filtered",
            "shifted_residual_profile_filtered",
        ],
        inplace=True,
    )

    return hourly_profiles


def aggregate_cems_to_fleet_for_residual_calc(
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year: int,
    datetime_col: str = "datetime_utc",
    ba_column_name: str = "ba_code",
    transmission_only: bool = False,
) -> pd.DataFrame:
    """Utility function for trying different BA aggregations in 930 and 923 data

    Args:
         cems (pd.DataFrame): Hourly, unit-level or subplant-level CEMS data
        partial_cems_subplant (pd.DataFrame): Hourly data for partial-subplant CEMS
            subplants
        partial_cems_plant (pd.DataFrame): Hourly data for partial-plant CEMS subplants
        plant_attributes (pd.DataFrame): table of plant-level static attributes
        year (int): The data year
        datetime_col (str, optional): Name of datetime column to use. Defaults to
            "datetime_utc".
        ba_column_name (str, optional): Whether to use "ba_code" or "ba_code_physical"
            to define fleet composition.. Defaults to "ba_code".
        transmission_only (bool, optional): Whether or not to only consider subplants
            that are connected to the transmission grid when calculating residual
            profiles. Defaults to False.

    Raises:
        UserWarning: If a fuel category cannot be assigned to a CEMS subplant

    Returns:
        pd.DataFrame: dataframe of all hourly CEMS data aggregated to the fleet level
    """

    # add the partial cems data
    cems_agg = pd.concat(
        [cems, partial_cems_subplant, partial_cems_plant], axis=0, copy=False
    )
    validation.validate_unique_datetimes(
        year, cems_agg, "cems_for_residual", ["plant_id_eia", "subplant_id"]
    )

    # Assign CEMS-subplant level data to a fleet based on the subplant-specific,
    # capacity-based primary fuel, and using the EIA-930 fuel categories. This should
    # best match the way data is aggregated to fleet in 930.
    cems_agg = assign_fleet_to_subplant_data(
        subplant_data=cems_agg,
        plant_attributes_table=plant_attributes,
        primary_fuel_table=primary_fuel_table,
        year=year,
        ba_col=ba_column_name,
        primary_fuel_col="subplant_primary_fuel_from_capacity_mw",
        fuel_category_col="fuel_category_eia930",
        other_attribute_cols=["distribution_flag"],
    )

    if transmission_only:
        cems_agg = cems_agg[cems_agg["distribution_flag"] is False]

    cems_agg = (
        cems_agg.groupby(
            [ba_column_name, "fuel_category_eia930", datetime_col], dropna=False
        )["net_generation_mwh"]
        .sum()
        .reset_index()
    )

    if ba_column_name == "ba_code_physical":
        cems_agg = cems_agg.rename(columns={"ba_code_physical": "ba_code"})

    # clean up the eia930 data before merging
    cems_agg = cems_agg.rename(columns={"fuel_category_eia930": "fuel_category"})

    # rename the net generation column to cems profile
    cems_agg = cems_agg.rename(columns={"net_generation_mwh": "cems_profile"})

    return cems_agg


def aggregate_non_930_fuel_categories(cems, plant_attributes):
    # get a list of the fuel categories not in EIA-930
    fuel_categories_not_in_eia930 = list(
        set(plant_attributes.fuel_category.unique())
        - set(plant_attributes.fuel_category_eia930.unique())
    )

    # create a new dataframe with only these fuels grouped together
    cems_profiles_for_non_930_fuels = cems.loc[
        cems["fuel_category"].isin(fuel_categories_not_in_eia930),
        ["ba_code", "fuel_category", "datetime_utc", "net_generation_mwh"],
    ]
    cems_profiles_for_non_930_fuels = (
        cems_profiles_for_non_930_fuels.groupby(
            ["ba_code", "fuel_category", "datetime_utc"], dropna=False
        )
        .sum()
        .reset_index()
    )
    return cems_profiles_for_non_930_fuels


def calculate_residual(
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
    eia930_data: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year: int,
    transmission_only: bool = False,
    ba_column_name: str = "ba_code",
) -> pd.DataFrame:
    """Creates a dataframe of residual hourly profiles for each fleet.

    NOTE: We do not return the CEMS profile used for the residual calculation. The back-
    stop CEMS value will be calculated in add_cems_backstop_profile()

    Args:
        cems (pd.DataFrame): Hourly, unit- or subplant-level CEMS data
        partial_cems_subplant (pd.DataFrame): Hourly data for partial-subplant CEMS
            subplants
        partial_cems_plant (pd.DataFrame): Hourly data for partial-plant CEMS subplants
        eia930_data (pd.DataFrame): Hourly EIA-930 net generation data for each fleet
        plant_attributes (pd.DataFrame): table of plant-level static attributes
        plant_attributes (pd.DataFrame): Static plant attributes table
        year (int): The data year
        transmission_only (bool, optional): Whether or not to only consider subplants
            that are connected to the transmission grid when calculating residual
            profiles. Defaults to False.
        ba_column_name (str, optional): Whether to use "ba_code" or "ba_code_physical"
            to define fleet composition. Defaults to "ba_code".

    Returns:
        pd.DataFrame: Dataframe of hourly profiles, containing columns
    """

    cems_agg = aggregate_cems_to_fleet_for_residual_calc(
        cems,
        partial_cems_subplant,
        partial_cems_plant,
        plant_attributes,
        primary_fuel_table,
        year,
        "datetime_utc",
        ba_column_name,
        transmission_only,
    )

    # clean up the eia930 data before merging
    eia930_data = eia930_data.rename(
        columns={
            "fuel_category_eia930": "fuel_category",
            "net_generation_mwh_930": "eia930_profile",
        }
    )
    # only keep rows where local datetime is in the current year
    eia930_data = eia930_data[eia930_data["datetime_local"].str[:4].astype(int) == year]

    # combine the data from both sources together
    combined_data = eia930_data.merge(
        cems_agg,
        how="left",
        on=["ba_code", "fuel_category", "datetime_utc"],
        validate="1:1",
    )

    # if there is no cems data for a ba-fuel, and there is eia profile data replace
    # missing values with zero. This means that the profile will match the EIA profile
    combined_data.loc[~combined_data["eia930_profile"].isna(), "cems_profile"] = (
        combined_data.loc[
            ~combined_data["eia930_profile"].isna(), "cems_profile"
        ].fillna(0)
    )

    combined_data = calculate_scaled_residual(combined_data)
    combined_data = calculate_shifted_residual(combined_data)

    # calculate the residual
    combined_data["residual_profile"] = (
        combined_data["eia930_profile"] - combined_data["cems_profile"]
    )

    return combined_data[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "eia930_profile",
            "residual_profile",
            "scaled_residual_profile",
            "shifted_residual_profile",
        ]
    ]


def calculate_scaled_residual(combined_data: pd.DataFrame) -> pd.DataFrame:
    """Scaling:
    If the residual is ever negative, we want to scale the cems net generation data to
    always be less than or equal to the 930 net generation.
    To do this, we'll try scaling the data as a percentage:
        1. For each hour, calculate the ratio between 930 NG and CEMS NG.
        2. For each BA-fuel, find the minimum ratio. If the minimum ratio is >= 1,
            it means that 930 is always greater than CEMS and doesn't need to be
            scaled. For any BA-fuels where the ratio is < 1, we will use this as a
            scaling factor to scale the CEMS data such that the scaled data is
            always <= the 930 data
        3. Multiply all hourly CEMS values by the scaling factor

    Args:
        combined_data (pd.DataFrame): dataframe containing combined CEMS and EIA-930
            data from previous steps

    Returns:
        pd.DataFrame: combined_data with an added "scaled_residual_profile" column
    """
    # Find scaling factor
    # only keep data where the cems data is greater than zero
    scaling_factors = combined_data[combined_data["cems_profile"] > 0].copy()
    # calculate the ratio of 930 net generation to cems net generation
    # if correct, ratio should be >=1
    scaling_factors["scaling_factor"] = (
        scaling_factors["eia930_profile"] / scaling_factors["cems_profile"]
    )
    # find the minimum ratio for each ba-fuel
    scaling_factors = (
        scaling_factors.groupby(
            ["ba_code", "fuel_category", "report_date"], dropna=False
        )["scaling_factor"]
        .min()
        .reset_index()
    )

    # only keep scaling factors < 1, which means the data needs to be scaled
    scaling_factors = scaling_factors[
        (scaling_factors["scaling_factor"] < 1)
        & (scaling_factors["scaling_factor"] > 0)
    ]

    # merge the scaling factor into the combined data
    # for any BA-fuels without a scaling factor, fill with 1 (scale to 100% of the origina data)
    combined_data = combined_data.merge(
        scaling_factors,
        how="left",
        on=["ba_code", "fuel_category", "report_date"],
        validate="m:1",
    )
    combined_data["scaling_factor"] = combined_data["scaling_factor"].fillna(1)

    # calculate the scaled cems data
    combined_data["cems_profile_scaled"] = (
        combined_data["cems_profile"] * combined_data["scaling_factor"]
    )

    # calculate the residual
    combined_data["scaled_residual_profile"] = (
        combined_data["eia930_profile"] - combined_data["cems_profile_scaled"]
    )

    return combined_data


def calculate_shifted_residual(combined_data):
    # Find scaling factor
    # only keep data where the cems data is greater than zero
    scaling_factors = combined_data[combined_data["cems_profile"] != 0].copy()
    # calculate the ratio of 930 net generation to cems net generation
    # if correct, ratio should be >=1
    scaling_factors["shift_factor"] = (
        scaling_factors["eia930_profile"] - scaling_factors["cems_profile"]
    )
    # find the minimum factor for each ba-fuel
    scaling_factors = (
        scaling_factors.groupby(
            ["ba_code", "fuel_category", "report_date"], dropna=False
        )["shift_factor"]
        .min()
        .reset_index()
    )

    # only keep scaling factors < 0, which means the data needs to be scaled
    scaling_factors = scaling_factors[scaling_factors["shift_factor"] < 0]

    # merge the scaling factor into the combined data
    # for any BA-fuels without a scaling factor, fill with 1 (scale to 100% of the origina data)
    combined_data = combined_data.merge(
        scaling_factors,
        how="left",
        on=["ba_code", "fuel_category", "report_date"],
        validate="m:1",
    )
    combined_data["shift_factor"] = combined_data["shift_factor"].fillna(0)

    # calculate the scaled cems data
    combined_data["cems_profile_shifted"] = (
        combined_data["cems_profile"] + combined_data["shift_factor"]
    )

    # calculate the residual
    combined_data["shifted_residual_profile"] = (
        combined_data["eia930_profile"] - combined_data["cems_profile_shifted"]
    )

    return combined_data


def create_flat_profile(report_date, ba, fuel):
    year = report_date.year

    df_temporary = pd.DataFrame(
        index=pd.date_range(
            start=f"{year - 1}-12-31 00:00:00",
            end=f"{year + 1}-01-01 23:00:00",
            freq="h",
            tz="UTC",
            name="datetime_utc",
        ),
        columns=["ba_code", "fuel_category"],
    ).reset_index()

    df_temporary["imputed_profile"] = 1.0
    df_temporary["datetime_local"] = ""
    df_temporary["datetime_local"] = (
        df_temporary["datetime_utc"]
        .dt.tz_convert(load_data.ba_timezone(ba=ba, type="local"))
        .astype(str)
    )
    # only keep data for which the local datetime is in the current year
    df_temporary = df_temporary[
        df_temporary["datetime_local"].str[:4].astype(int) == year
    ]

    # create a report date column
    df_temporary["report_date"] = df_temporary["datetime_local"].str[:7]
    df_temporary["report_date"] = pd.to_datetime(df_temporary["report_date"]).astype(
        "datetime64[s]"
    )

    # only keep the report dates that match
    df_temporary = df_temporary[df_temporary["report_date"] == report_date]

    df_temporary["ba_code"] = ba
    df_temporary["fuel_category"] = fuel

    return df_temporary


def impute_missing_hourly_profiles(
    monthly_eia_data_to_shape: pd.DataFrame,
    residual_profiles: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Identifies which hourly profiles are missing after calculating residual_profiles,
    and imputes profiles for wind and solar, as well as flat profiles for "baseload" gen


    Args:
        monthly_eia_data_to_shape (pd.DataFrame): used to identify the complete set of
            fleet-months for which we have data
        residual_profiles (pd.DataFrame): Calculated from calculate_residual(), the df
            we evaluate for missing fleet-months
        plant_attributes (pd.DataFrame): Used to assign fleet keys
        primary_fuel_table (pd.DataFrame): used to assign fleet keys
        year (int): the data year

    Raises:
        UserWarning: if hourly_profiles contains an unrecognized fuel category

    Returns:
        pd.DataFrame: new hourly_profiles dataframe that includes complete data for all
            fleet-months, including calculated residual profiles and imputed backstop
            profiles
    """

    # round the data to the nearest tenth
    # residual_profiles["profile"] = residual_profiles["profile"].round(1)

    missing_profiles = identify_missing_profiles(
        monthly_eia_data_to_shape,
        residual_profiles,
        plant_attributes,
        primary_fuel_table,
        year,
    )

    # load information about directly interconnected balancing authorities (DIBAs)
    # this will help us fill profiles using data from nearby BAs
    dibas = load_data.load_diba_data(year)

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    for index, row in missing_profiles.iterrows():
        ba = row["ba_code"]
        fuel = row["fuel_category"]
        report_date = row["report_date"]

        # for wind and solar, average the wind and solar generation profiles from
        # nearby interconnected BAs
        if fuel in ["wind", "solar"]:
            # get a list of diba located in the same region and located in the same time zone
            ba_dibas = list(
                dibas.loc[
                    (dibas["ba_code"] == ba)
                    & (dibas["ba_region"] == dibas["diba_region"])
                    & (dibas["timezone_local"] == dibas["timezone_local_diba"]),
                    "diba_code",
                ].unique()
            )
            if len(ba_dibas) > 0:
                df_temporary = average_diba_wind_solar_profiles(
                    residual_profiles, ba, fuel, report_date, ba_dibas
                )
            # if there are no neighboring DIBAs, calculate a national average profile
            else:
                df_temporary = average_national_wind_solar_profiles(
                    residual_profiles, ba, fuel, report_date
                )

        # certain fuels we assume would be operated as baseload
        elif fuel in ["geothermal", "biomass", "waste", "nuclear"]:
            # assign a flat profile
            df_temporary = create_flat_profile(report_date, ba, fuel)
            df_temporary["imputation_method"] = "assumed_flat"

        # For now assume hydro is dispatched with a flat profile
        # TODO improve this assumption see: https://github.com/singularity-energy/open-grid-emissions/issues/37
        elif fuel in ["hydro"]:
            df_temporary = create_flat_profile(report_date, ba, fuel)
            df_temporary["imputation_method"] = "assumed_flat"
        # for any other fossil resources, use a flat profile
        # TODO: we need to improve this method
        # see: https://github.com/singularity-energy/open-grid-emissions/issues/96
        elif fuel in ["natural_gas", "coal", "petroleum"]:
            df_temporary = create_flat_profile(report_date, ba, fuel)
            df_temporary["imputation_method"] = "assumed_flat"
        elif fuel in ["other"]:
            # assign a flat profile
            df_temporary = create_flat_profile(report_date, ba, fuel)
            df_temporary["imputation_method"] = "assumed_flat"
        else:
            raise UserWarning(f"Fuel category {fuel} not recognized.")

        hourly_profiles_to_add.append(df_temporary)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True, copy=False
    )

    hourly_profiles = pd.concat(
        [residual_profiles, hourly_profiles_to_add], axis=0, copy=False
    )

    # convert to datetime64[s] format
    hourly_profiles["datetime_utc"] = (
        pd.to_datetime(hourly_profiles["datetime_utc"])
        .dt.tz_localize(None)
        .astype("datetime64[s]")
        .dt.tz_localize("UTC")
    )
    validation.validate_unique_datetimes(
        year, hourly_profiles, "hourly_profiles", ["ba_code", "fuel_category"]
    )

    return hourly_profiles


def identify_missing_profiles(
    monthly_eia_data_to_shape: pd.DataFrame,
    residual_profiles: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year,
) -> pd.DataFrame:
    """Identifies the fleet-months for which there is no profile data in residual_profiles

    Args:
        monthly_eia_data_to_shape (pd.DataFrame): used to determine the complete set of
            fleet-months for which we should have a profile.
        residual_profiles (pd.DataFrame): used to determine which profiles are available
        plant_attributes (pd.DataFrame): used to identify fleet keys
        primary_fuel_table (pd.DataFrame): used to identify fleet keys
        year (int): the data year

    Returns:
        pd.DataFrame: table with three columns (ba_code, fuel_category, report_date)
            that list the combinations of these keys for which a profile is missing in
            residual_profiles
    """
    # drop ba fuel months where the reported profile is all zeros
    MONTHLY_GROUP_COLUMNS = [
        "ba_code",
        "fuel_category",
        "report_date",
    ]

    # determine for which BA-fuels we are missing residual profiles
    available_profiles = residual_profiles[MONTHLY_GROUP_COLUMNS].drop_duplicates()
    # assign fleet to the subplant-level monthly data
    monthly_eia_data_to_shape = assign_fleet_to_subplant_data(
        monthly_eia_data_to_shape,
        plant_attributes,
        primary_fuel_table,
        year,
        ba_col="ba_code",
        primary_fuel_col="subplant_primary_fuel_from_capacity_mw",
        fuel_category_col="fuel_category",
    )
    ba_fuel_to_distribute = monthly_eia_data_to_shape[
        MONTHLY_GROUP_COLUMNS
    ].drop_duplicates()
    missing_profiles = ba_fuel_to_distribute.merge(
        available_profiles,
        how="outer",
        on=MONTHLY_GROUP_COLUMNS,
        indicator="source",
        validate="1:1",
    )
    # identify ba fuel months where there is no data in the available residual profiles
    missing_profiles = missing_profiles[missing_profiles.source == "left_only"]
    missing_profiles.drop(columns="source", inplace=True)
    missing_profiles = missing_profiles.sort_values(by=MONTHLY_GROUP_COLUMNS)

    return missing_profiles


def average_diba_wind_solar_profiles(
    residual_profiles, ba, fuel, report_date, ba_dibas, validation_run=False
):
    # calculate the average generation profile for the fuel in all neighboring DIBAs
    df_temporary = residual_profiles[
        (residual_profiles["ba_code"].isin(ba_dibas))
        & (residual_profiles["fuel_category"] == fuel)
        & (residual_profiles["report_date"] == report_date)
    ]
    if len(df_temporary) == 0 and not validation_run:
        # if this error is raised, we might have to implement an approach that uses average values for the wider region
        logger.warning(f"There is no {fuel} data in the DIBAs for {ba}: {ba_dibas}")
        df_temporary = average_national_wind_solar_profiles(
            residual_profiles, ba, fuel, report_date
        )
    else:
        df_temporary = (
            df_temporary.groupby(
                ["fuel_category", "datetime_utc", "datetime_local", "report_date"],
                dropna=False,
            )["eia930_profile"]
            .mean()
            .reset_index()
        )
        df_temporary["ba_code"] = ba
        df_temporary["imputation_method"] = "DIBA_average"
        df_temporary = df_temporary.rename(
            columns={"eia930_profile": "imputed_profile"}
        )

    return df_temporary


def average_national_wind_solar_profiles(residual_profiles, ba, fuel, report_date):
    df_temporary = residual_profiles[
        (residual_profiles["fuel_category"] == fuel)
        & (residual_profiles["report_date"] == report_date)
    ].copy()
    # strip the time zone information so we can group by local time
    df_temporary["datetime_local"] = df_temporary["datetime_local"].str[:-6]
    df_temporary = (
        df_temporary.groupby(
            ["fuel_category", "datetime_local", "report_date"],
            dropna=False,
        )["eia930_profile"]
        .mean()
        .reset_index()
    )
    df_temporary["ba_code"] = ba
    df_temporary["imputation_method"] = "national_average"
    df_temporary = df_temporary.rename(columns={"eia930_profile": "imputed_profile"})

    # re-localize the datetime_local
    local_tz = load_data.ba_timezone(ba, "local")
    df_temporary["datetime_local"] = pd.to_datetime(
        df_temporary["datetime_local"]
    ).astype("datetime64[s]")
    df_temporary["datetime_local"] = (
        df_temporary["datetime_local"]
        .dt.tz_localize(local_tz, nonexistent="NaT", ambiguous="NaT")
        .ffill()
    )
    df_temporary["datetime_utc"] = df_temporary["datetime_local"].dt.tz_convert("UTC")

    # drop duplicate datetimes around DST
    df_temporary = df_temporary.drop_duplicates(subset=["datetime_utc"], keep="first")

    return df_temporary


def validate_wind_solar_imputation(hourly_profiles, year):
    """Creates a table showing cross-validaton results of the wind and solar profile imputation method"""

    # calculate the results and merge together
    diba_results = validate_diba_imputation_method(hourly_profiles, year)
    nationaal_results = validate_national_imputation_method(hourly_profiles)

    imputation_results = diba_results.merge(
        nationaal_results, how="outer", on=["fuel_category", "ba_code"], validate="1:1"
    )

    return imputation_results


def validate_diba_imputation_method(hourly_profiles, year):
    """Validates the method for imputing missing wind and solar profiles.

    Calculates an imputed profile for regions where we have actual wind and solar profiles,
    then calculates how well each imputed profile is correlated with the actual profile.
    Calculates the correlation for each month, then calculates an annual average correlation coefficient.
    """

    # only keep wind and solar data
    data_to_validate = hourly_profiles[
        (hourly_profiles["fuel_category"].isin(["wind", "solar"]))
        & (~hourly_profiles["eia930_profile"].isna())
    ]
    data_to_validate = data_to_validate[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "eia930_profile",
        ]
    ]

    profiles_to_impute = data_to_validate[
        ["ba_code", "fuel_category", "report_date"]
    ].drop_duplicates()

    profiles_to_impute = profiles_to_impute[
        profiles_to_impute["report_date"].dt.year == year
    ]

    dibas = load_data.load_diba_data(year)

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    for index, row in profiles_to_impute.iterrows():
        ba = row["ba_code"]
        fuel = row["fuel_category"]
        report_date = row["report_date"]

        # for wind and solar, average the wind and solar generation profiles from
        # nearby interconnected BAs
        if fuel in ["wind", "solar"]:
            # get a list of diba located in the same region and located in the same time zone
            ba_dibas = list(
                dibas.loc[
                    (dibas["ba_code"] == ba)
                    & (dibas["ba_region"] == dibas["diba_region"])
                    & (dibas["timezone_local"] == dibas["timezone_local_diba"]),
                    "diba_code",
                ].unique()
            )
            if len(ba_dibas) > 0:
                df_temporary = average_diba_wind_solar_profiles(
                    data_to_validate, ba, fuel, report_date, ba_dibas, True
                )

                hourly_profiles_to_add.append(df_temporary)
            # if there are no neighboring DIBAs, calculate a national average profile
            else:
                pass

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True, copy=False
    )

    # merge the imputed data with the actual data
    compare_method = data_to_validate.merge(
        hourly_profiles_to_add,
        how="left",
        on=[
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "ba_code",
        ],
        validate="1:1",
    )

    # calculate the correlation coefficient for each fleet-month
    compare_method = (
        compare_method.groupby(["fuel_category", "report_date", "ba_code"])
        .corr(numeric_only=True)
        .reset_index()
    )
    compare_method = compare_method[compare_method["level_3"] == "eia930_profile"]

    # calculate the annual average correlation coefficent for each month
    compare_method = (
        compare_method.groupby(["fuel_category", "ba_code"])["imputed_profile"]
        .mean()
        .reset_index()
    )

    compare_method = compare_method.rename(
        columns={"imputed_profile": "diba_method_correlation_coefficient"}
    )

    return compare_method


def validate_national_imputation_method(hourly_profiles):
    # only keep wind and solar data
    data_to_validate = hourly_profiles[
        (hourly_profiles["fuel_category"].isin(["wind", "solar"]))
        & (~hourly_profiles["eia930_profile"].isna())
    ]
    data_to_validate = data_to_validate[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "eia930_profile",
        ]
    ]

    profiles_to_impute = data_to_validate[
        ["ba_code", "fuel_category", "report_date"]
    ].drop_duplicates()

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    for index, row in profiles_to_impute.iterrows():
        ba = row["ba_code"]
        fuel = row["fuel_category"]
        report_date = row["report_date"]

        # for wind and solar, average the wind and solar generation profiles from
        # nearby interconnected BAs
        if fuel in ["wind", "solar"]:
            # get a list of diba located in the same region and located in the same time zone
            df_temporary = average_national_wind_solar_profiles(
                data_to_validate, ba, fuel, report_date
            )

        hourly_profiles_to_add.append(df_temporary)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True, copy=False
    )

    # merge the imputed data with the actual data
    compare_method = data_to_validate.merge(
        hourly_profiles_to_add,
        how="left",
        on=["fuel_category", "datetime_utc", "report_date", "ba_code"],
        validate="1:1",
    )

    # calculate the correlation coefficient for each fleet-month
    compare_method = (
        compare_method.groupby(["fuel_category", "report_date", "ba_code"])
        .corr(numeric_only=True)
        .reset_index()
    )
    compare_method = compare_method[compare_method["level_3"] == "eia930_profile"]

    # calculate the annual average correlation coefficent for each month
    compare_method = (
        compare_method.groupby(["fuel_category", "ba_code"])["imputed_profile"]
        .mean()
        .reset_index()
    )

    compare_method = compare_method.rename(
        columns={"imputed_profile": "national_method_correlation_coefficient"}
    )

    return compare_method


def add_cems_backstop_profile(
    hourly_profiles: pd.DataFrame,
    cems: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year,
) -> pd.DataFrame:
    """Adds a "cems_profile" column to hourly_profiles to use as a backstop where
    residuals are not available.

    This function only adds a profile if there are at least 4 plants in the fleet so
    that the profile is not overly-fit to a specific plant.

    We also use OGE fuel categories to aggregate to the fleet so that if data for
    fuel categories that do not exist in EIA-930 (ie waste, biomass) are in CEMS, we
    have profiles to use rather than using the flat backstop.

    Args:
        hourly_profiles (pd.DataFrame): dataframe calculated by
            impute_missing_hourly_profiles(). "cems_profile" will be added to this.
        cems (pd.DataFrame): Used to create cems profiles
        plant_attributes (pd.DataFrame): Used to identify fleet keys
        primary_fuel_table (pd.DataFrame): Used to identify fleet keys

    Returns:
        pd.DataFrame: hourly_profiles with a "cems_profile" column added
    """
    cems_ba_fuel = assign_fleet_to_subplant_data(
        cems,
        plant_attributes,
        primary_fuel_table,
        year,
        ba_col="ba_code",
        primary_fuel_col="subplant_primary_fuel_from_capacity_mw",
        fuel_category_col="fuel_category",
    ).rename(columns={"fuel_category": "fuel_category"})

    # Count unique plants: after grouping by BA we will remove where n_unique_plants < 3
    cems_count = (
        cems_ba_fuel.groupby(["ba_code", "fuel_category", "report_date"], dropna=False)[
            "plant_id_eia"
        ]
        .nunique()
        .reset_index()
        .rename(columns={"plant_id_eia": "n_unique_plants"})
    )

    # grouping by datetime_utc and report_date will lead to some duplicate datetime
    # values since report_date is based on datetime_local, and some plants in a ba may
    # be located in different timezones. We will remove these duplicates through a later
    # groupby operation once we no longer need report_date as a merge key
    cems_ba_fuel = (
        cems_ba_fuel.groupby(
            ["ba_code", "fuel_category", "datetime_utc", "report_date"], dropna=False
        )["net_generation_mwh"]
        .sum()
        .reset_index()
    )

    # Remove data where too few plants
    cems_ba_fuel = cems_ba_fuel.merge(
        cems_count,
        how="left",
        on=["ba_code", "fuel_category", "report_date"],
        validate="m:1",
    )
    cems_ba_fuel = cems_ba_fuel[cems_ba_fuel["n_unique_plants"] > 3]
    cems_ba_fuel.drop(columns=["n_unique_plants"], inplace=True)

    # remove months where there is zero generation reported
    months_with_zero_data = (
        cems_ba_fuel.groupby(["ba_code", "fuel_category", "report_date"], dropna=False)[
            "net_generation_mwh"
        ]
        .sum()
        .reset_index()
    )
    months_with_zero_data = months_with_zero_data.loc[
        months_with_zero_data["net_generation_mwh"] == 0,
        ["ba_code", "fuel_category", "report_date"],
    ]
    cems_ba_fuel = cems_ba_fuel.merge(
        months_with_zero_data,
        how="outer",
        on=["ba_code", "fuel_category", "report_date"],
        indicator="zero_filter",
        validate="m:1",
    )
    cems_ba_fuel = cems_ba_fuel[cems_ba_fuel["zero_filter"] != "both"]
    cems_ba_fuel.drop(columns=["zero_filter"], inplace=True)

    # remove months where the cems profile has negative values
    months_with_negative_data = (
        cems_ba_fuel.groupby(["ba_code", "fuel_category", "report_date"], dropna=False)[
            "net_generation_mwh"
        ]
        .min()
        .reset_index()
    )
    months_with_negative_data = months_with_negative_data.loc[
        months_with_negative_data["net_generation_mwh"] < 0,
        ["ba_code", "fuel_category", "report_date"],
    ]
    cems_ba_fuel = cems_ba_fuel.merge(
        months_with_negative_data,
        how="outer",
        on=["ba_code", "fuel_category", "report_date"],
        indicator="neg_filter",
        validate="m:1",
    )
    cems_ba_fuel = cems_ba_fuel[cems_ba_fuel["neg_filter"] != "both"]
    cems_ba_fuel.drop(columns=["neg_filter", "report_date"], inplace=True)

    # remove duplicate datetime values
    cems_ba_fuel = (
        cems_ba_fuel.groupby(
            ["ba_code", "fuel_category", "datetime_utc"], dropna=False
        )["net_generation_mwh"]
        .sum()
        .reset_index()
    ).rename(columns={"net_generation_mwh": "cems_profile"})

    # fill missing cems profile data
    hourly_profiles = hourly_profiles.merge(
        cems_ba_fuel[["ba_code", "fuel_category", "datetime_utc", "cems_profile"]],
        how="left",
        on=["ba_code", "fuel_category", "datetime_utc"],
        validate="1:1",
    )

    return hourly_profiles


def convert_profile_to_percent(hourly_profiles, group_keys, columns_to_convert):
    """converts hourly timeseries profiles from absolute mwh to percentage of monthly total mwh."""
    # convert the profile so that each hour is a percent of the monthly total
    monthly_group_columns = group_keys + ["report_date"]

    monthly_total = (
        hourly_profiles.groupby(monthly_group_columns, dropna=False)[columns_to_convert]
        .sum()
        .reset_index()
    )

    # merge the hourly totals back into the hourly data
    hourly_profiles = hourly_profiles.merge(
        monthly_total,
        how="left",
        on=monthly_group_columns,
        suffixes=(None, "_monthly_total"),
        validate="m:1",
    )
    for col in columns_to_convert:
        hourly_profiles[col] = (
            hourly_profiles[col] / hourly_profiles[f"{col}_monthly_total"]
        )
        hourly_profiles.drop(columns=[f"{col}_monthly_total"], inplace=True)

    return hourly_profiles


def combine_and_export_hourly_plant_data(
    year: int,
    cems: pd.DataFrame,
    partial_cems_subplant: pd.DataFrame,
    partial_cems_plant: pd.DataFrame,
    monthly_eia_data_to_shape: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    hourly_profiles: pd.DataFrame,
    path_prefix: str,
    skip_outputs: bool,
    region_to_group: str = "ba_code",
):
    """
    Exports files with hourly data for each individual plant, split up by region.

    Creating hourly records for all EIA plants in the US at once will cause memory errors
    on most computers, so we need to only shape data for one subset of plants at a time.
    This function shapes the EIA monthly data for one region at a time, combines it with
    the hourly CEMS data for that region, and exports the data as a csv file.

    Process:
        1. Combine all hourly CEMS data together and aggregate to plant level
        2. Add a BA code to the CEMS data to help with filtering
        3. Add fleet identifiers to the EIA data to help with filtering and shaping
        4. For each BA:
            a. Filter the EIA data to the region and shape the data to hourly
            b. Aggregate EIA data to the plant level
            c. Combine the hourly EIA plant data with the hourly CEMS plant data
            d. Export the data

    All of the inputs are dataframes containing data from the data pipeline except for `region_to_group`
    `region_to_group` identifying whether "ba_code" or "state" should be used to group the data. "ba_code" is the default.
    """

    # specify the names of the columns that we want to use to group the data
    key_columns = [
        "plant_id_eia",
        "datetime_utc",
    ]
    # specify the data columns that should be included in the output files
    data_columns_for_plant_export = [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "ch4_mass_lb",
        "n2o_mass_lb",
        "nox_mass_lb",
        "so2_mass_lb",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_for_electricity_adjusted",
    ]
    all_columns = key_columns + data_columns_for_plant_export

    # ensure that we will not be duplicating data when we combine
    validation.ensure_non_overlapping_data_from_all_sources(
        cems, partial_cems_subplant, partial_cems_plant, monthly_eia_data_to_shape
    )

    # combine all the CEMS data together and only keep relevant columns
    combined_cems_data = pd.concat(
        [cems, partial_cems_subplant, partial_cems_plant],
        axis=0,
        ignore_index=True,
        copy=False,
    )[[col for col in cems.columns if col in all_columns]]

    # Free memory after combining
    gc.collect()

    # aggregate combined CEMS data to the plant level
    combined_cems_data = (
        combined_cems_data.groupby(
            key_columns,
            dropna=False,
        )
        .sum(numeric_only=True)
        .reset_index()
    )

    # add ba_codes to the input data to help with regional filtering
    combined_cems_data = combined_cems_data.merge(
        plant_attributes[["plant_id_eia", region_to_group]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )

    # because we need to shape the EIA based on its fleet category, we need to assign
    # fleet based on the subplant fuel, not the plant fuel
    monthly_eia_data_to_shape = assign_fleet_to_subplant_data(
        monthly_eia_data_to_shape,
        plant_attributes,
        primary_fuel_table,
        year,
        ba_col="ba_code",
        primary_fuel_col="subplant_primary_fuel_from_capacity_mw",
        fuel_category_col="fuel_category",
    )

    # for each region, shape the EIA-only data, combine with CEMS data, and export
    for region in list(
        plant_attributes.sort_values(by=region_to_group)[region_to_group].unique()
    ):
        logger.info(f"Shaping hourly plant data for {region}")
        # filter each of the data sources to the region
        # Note: copy needed because we'll be modifying this subset
        eia_region = monthly_eia_data_to_shape[
            monthly_eia_data_to_shape[region_to_group] == region
        ].copy()

        # shape the EIA subplant-level data
        shaped_eia_region_data = shape_monthly_eia_data_as_hourly(
            eia_region, hourly_profiles, year
        )
        # validate that the shaped data contains no duplicate datetimes
        validation.validate_unique_datetimes(
            year,
            df=shaped_eia_region_data,
            df_name="shaped_eia_region_data",
            keys=["plant_id_eia", "subplant_id"],
        )

        # group the subplant data to plant
        shaped_eia_region_data = (
            shaped_eia_region_data.groupby(
                key_columns,
                dropna=False,
            )[
                [
                    col
                    for col in shaped_eia_region_data.columns
                    if col in data_columns_for_plant_export
                ]
            ]
            .sum(numeric_only=True)
            .reset_index()
        )

        # filter the CEMS data to the region and combine it with the EIA data
        # Note: copy needed because we'll be modifying this subset
        cems_region = combined_cems_data[
            combined_cems_data[region_to_group] == region
        ].copy()
        combined_plant_data = pd.concat(
            [
                cems_region,
                shaped_eia_region_data,
            ],
            axis=0,
            ignore_index=True,
            copy=False,
        )

        del (
            cems_region,
            shaped_eia_region_data,
        )

        # groupby plant in case some plant data was split between multiple dfs
        combined_plant_data = (
            combined_plant_data.groupby(key_columns, dropna=False)
            .sum(numeric_only=True)
            .reset_index()
        )

        # round the data columns to two decimal places
        combined_plant_data[data_columns_for_plant_export] = combined_plant_data[
            data_columns_for_plant_export
        ].round(2)

        # re-order columns
        combined_plant_data = combined_plant_data[all_columns]

        # validate that there are complete timeseries
        validation.check_for_complete_hourly_timeseries(
            df=combined_plant_data,
            df_name=f"plant_data/hourly/{region}",
            keys=["plant_id_eia"],
            period="year",
        )

        # write data
        output_data.output_to_results(
            combined_plant_data,
            year,
            region,
            "plant_data/hourly/",
            path_prefix,
            skip_outputs,
            include_metric=False,
        )


def get_shaped_plant_id_from_ba_fuel(df: pd.DataFrame) -> pd.DataFrame:
    """Return artificial plant code. Max real plant is 64663
    Our codes look like 9BBBFF where BBB is the three digit BA number and FF is the
    two-digit fuel number

    Args:
        df (pd.DataFrame): df to add shaped plant IDs to. must contain `ba_code` and
            `fuel_category`

    Returns:
        pd.DataFrame: df with added "shaped_plant_id" column
    """
    # specify the ba numbers with leading zeros
    FUEL_NUMBERS = {
        "biomass": "01",
        "coal": "02",
        "geothermal": "03",
        "hydro": "04",
        "natural_gas": "05",
        "nuclear": "06",
        "other": "07",
        "petroleum": "08",
        "solar": "09",
        "storage": "10",
        "waste": "11",
        "wind": "12",
    }

    # load the ba reference table with all of the ba number ids
    ba_numbers = pd.read_csv(reference_table_folder("ba_reference.csv"))[
        ["ba_code", "ba_number"]
    ]
    # reformat the number with leading zeros
    ba_numbers["ba_number"] = ba_numbers["ba_number"].astype(str).str.zfill(3)
    # convert to a dictionary
    ba_numbers = dict(zip(ba_numbers["ba_code"], ba_numbers["ba_number"]))

    # make sure the ba codes are strings
    df["ba_code"] = df["ba_code"].astype(str)
    # create a new column with the shaped plant ids
    df["shaped_plant_id"] = df.apply(
        lambda row: f"9{ba_numbers[row['ba_code']]}{FUEL_NUMBERS[row['fuel_category_for_shaping']]}",
        axis=1,
    )
    # convert to an int32 column
    df["shaped_plant_id"] = df["shaped_plant_id"].astype("Int32")

    return df


def aggregate_eia_data_to_fleet(
    monthly_eia_data_to_shape: pd.DataFrame,
    plant_attributes: pd.DataFrame,
    primary_fuel_table: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Aggregates monthly EIA-923 data to the fleet level prior to assigning an hourly
    profile

    Given cleaned monthly EIA-923 data and plant attributes, aggregate to BA-fuel
    using artificial plant IDs 9XXXYYY where XXX=BA code (see `ba_reference.csv`)
    and YY=fuel (see `impute_hourly_profiles.get_shaped_plant_id_from_ba_fuel`)

    Add new artificial plants to plant_attributes frame.

    Args:
        monthly_eia_data_to_shape (pd.DataFrame): data that will be aggregated and
            shaped
        plant_attributes (pd.DataFrame): for assigning fleet keys
        primary_fuel_table (pd.DataFrame): for assigning fleet keys
        year (int): the data year

    Returns:
        pd.DataFrame: returns eia_agg, which is aggregated to the fleet-month level
            using shaped_plant_id
    """

    # NOTE: currently using ba_code, could alternatively use ba_code_physical
    # First, add capacity-based fuel category for shaping the data
    eia_agg = assign_fleet_to_subplant_data(
        monthly_eia_data_to_shape,
        plant_attributes,
        primary_fuel_table,
        year,
        ba_col="ba_code",
        primary_fuel_col="subplant_primary_fuel_from_capacity_mw",
        fuel_category_col="fuel_category",
    )
    eia_agg = eia_agg.rename(columns={"fuel_category": "fuel_category_for_shaping"})

    # now, add fuel-based fuel category for fleet aggregation
    eia_agg = assign_fleet_to_subplant_data(
        eia_agg,
        plant_attributes,
        primary_fuel_table,
        year,
        ba_col="ba_code",
        primary_fuel_col="subplant_primary_fuel",
        fuel_category_col="fuel_category",
    )

    # create a column with shaped plant ids
    eia_agg = get_shaped_plant_id_from_ba_fuel(eia_agg)

    # Group
    eia_agg = (
        eia_agg.groupby(
            [
                "shaped_plant_id",
                "ba_code",
                "report_date",
                "fuel_category",
                "fuel_category_for_shaping",
            ],
            dropna=False,
        )
        .sum(numeric_only=True)
        .reset_index()
        .drop(columns=["plant_id_eia", "subplant_id"])
        .rename(columns={"shaped_plant_id": "plant_id_eia"})
    )

    return eia_agg


def shape_monthly_eia_data_as_hourly(
    monthly_eia_data_to_shape: pd.DataFrame,
    hourly_profiles: pd.DataFrame,
    year: int,
    fuel_category_col_for_shaping: str = "fuel_category",
) -> pd.DataFrame:
    """Assigns an hourly profile to monthly-level EIA data.

    Can be used to shape subplant-level data (if monthly_eia_data_to_shape contains a
    subplant_id column), or shape fleet (BA-fuel) level data.

    The fuel_category column from both monthly_eia_data_to_shape and hourly_profiles
    is used to identify the correct profile. This column can be based on any fuel
    assignment method before passing the dataframe to this function.

    Args:
        monthly_eia_data_to_shape (pd.DataFrame): The monthly data to shape
        hourly_profiles (pd.DataFrame): The profiles used to shape the other df
        year (int): the data year, only used in the validation function
        fuel_category_col_for_shaping (str): the name of the fuel column for shaping
            the data. When shaping the data in step 16 we assign a
            fuel_category_for_shaping column

    Returns:
        pd.DataFrame: plant- or fleet-level hourly activity data
    """

    # merge the hourly profiles into each plant-month
    # this will create a row for each hour in a month
    shaped_monthly_data = monthly_eia_data_to_shape.merge(
        hourly_profiles[
            [
                "ba_code",
                "fuel_category",
                "datetime_utc",
                "report_date",
                "profile",
                "flat_profile",
                "profile_method",
            ]
        ].rename(columns={"fuel_category": fuel_category_col_for_shaping}),
        how="left",
        on=["report_date", fuel_category_col_for_shaping, "ba_code"],
        validate="m:m",
    )

    # plant-months where there is negative net generation, assign a flat profile
    shaped_monthly_data.loc[
        shaped_monthly_data["net_generation_mwh"] < 0, "profile"
    ] = shaped_monthly_data.loc[
        shaped_monthly_data["net_generation_mwh"] < 0, "flat_profile"
    ]
    # update the method column
    shaped_monthly_data.loc[
        shaped_monthly_data["net_generation_mwh"] < 0, "profile_method"
    ] = "flat_negative_generation"

    # shape the data by multiplying the monthly activity data by the fraction of total
    # monthly activity occuring in that hour ("profile")
    for column in DATA_COLUMNS:
        if column in shaped_monthly_data.columns:
            shaped_monthly_data[column] = (
                shaped_monthly_data[column] * shaped_monthly_data["profile"]
            )
        else:
            pass

    # re order the columns
    column_order = [
        "plant_id_eia",
        "subplant_id",
        "ba_code",
        "fuel_category",
        "datetime_utc",
        "report_date",
        "profile_method",
    ] + DATA_COLUMNS

    # re-order and drop intermediate columns
    shaped_monthly_data = shaped_monthly_data[
        [col for col in column_order if col in shaped_monthly_data.columns]
    ]

    # validate that the shaping did not alter data at the monthly level
    validation.validate_shaped_totals(
        shaped_monthly_data,
        monthly_eia_data_to_shape,
        year,
        group_keys=["ba_code", "fuel_category"],
    )

    return shaped_monthly_data


def shape_partial_cems_plants(cems, eia923_allocated, year):
    """Shapes the monthly data for subplants where partial plant data is available in
    CEMS."""

    SUBPLANT_KEYS = ["report_date", "plant_id_eia", "subplant_id"]

    # identify all of the partial cems plants and group by subplant-month
    eia_data_to_shape = eia923_allocated.loc[
        eia923_allocated.hourly_data_source == "partial_cems_plant"
    ].copy()

    # if there is no data in the partial cems dataframe, skip.
    if len(eia_data_to_shape) > 0:
        # group the eia data by subplant
        eia_data_to_shape = (
            eia_data_to_shape.groupby(SUBPLANT_KEYS, dropna=False)[DATA_COLUMNS]
            .sum()
            .reset_index()
        )

        # get a list of plant ids for plants with partial plant data in CEMS
        partial_cems_plant_ids = list(eia_data_to_shape.plant_id_eia.unique())

        # get the hourly cems data for the partial plants and aggregate by plant-hour to use for shaping
        partial_cems_profiles = (
            cems[cems["plant_id_eia"].isin(partial_cems_plant_ids)]
            .groupby(["plant_id_eia", "report_date", "datetime_utc"], dropna=False)[
                ["gross_generation_mwh", "fuel_consumed_mmbtu"]
            ]
            .sum()
            .reset_index()
        )
        # add a column for flat profiles
        partial_cems_profiles["flat_profile"] = 1
        # convert the profiles to a percent
        partial_cems_profiles = convert_profile_to_percent(
            partial_cems_profiles,
            group_keys=["plant_id_eia"],
            columns_to_convert=[
                "gross_generation_mwh",
                "fuel_consumed_mmbtu",
                "flat_profile",
            ],
        )
        partial_cems_profiles = partial_cems_profiles.rename(
            columns={
                "gross_generation_mwh": "generation_profile",
                "fuel_consumed_mmbtu": "fuel_profile",
            }
        )

        # prepare the profiles
        # generation or fuel profiles will be missing for a month if there was zero generation or fuel reported for that month
        # if we are missing a profile, try filling it with the other profile (e.g. fill generation with fuel profile)
        partial_cems_profiles["generation_profile"] = partial_cems_profiles[
            "generation_profile"
        ].fillna(partial_cems_profiles["fuel_profile"])
        partial_cems_profiles["fuel_profile"] = partial_cems_profiles[
            "fuel_profile"
        ].fillna(partial_cems_profiles["generation_profile"])
        # if the profile is still missing, fill it using a flat profile
        partial_cems_profiles["generation_profile"] = partial_cems_profiles[
            "generation_profile"
        ].fillna(partial_cems_profiles["flat_profile"])
        partial_cems_profiles["fuel_profile"] = partial_cems_profiles[
            "fuel_profile"
        ].fillna(partial_cems_profiles["flat_profile"])

        # merge the profiles into the monthly data
        shaped_partial_plants = eia_data_to_shape.merge(
            partial_cems_profiles,
            how="left",
            on=["plant_id_eia", "report_date"],
            validate="m:m",
        )

        # where monthly net generation is negative, replace the generation profile with a flat profile
        shaped_partial_plants.loc[
            shaped_partial_plants["net_generation_mwh"] < 0, "generation_profile"
        ] = shaped_partial_plants.loc[
            shaped_partial_plants["net_generation_mwh"] < 0, "flat_profile"
        ]

        # check that no profiles contain NA values
        missing_profiles = shaped_partial_plants[
            shaped_partial_plants["generation_profile"].isna()
            | shaped_partial_plants["fuel_profile"].isna()
        ]

        if len(missing_profiles) > 0:
            logger.warning(
                "Certain partial CEMS plants are missing hourly profile data. This will result in inaccurate results"
            )
        # check that all profiles add to 1 for each month
        incorrect_profiles = (
            shaped_partial_plants.groupby(SUBPLANT_KEYS)[
                ["generation_profile", "fuel_profile"]
            ]
            .sum()
            .reset_index()
        )
        incorrect_profiles = incorrect_profiles[
            (~np.isclose(incorrect_profiles["generation_profile"], 1))
            | (~np.isclose(incorrect_profiles["fuel_profile"], 1))
        ]

        if len(incorrect_profiles) > 0:
            logger.warning(
                "Certain partial CEMS profiles do not add to 100%. This will result in inaccurate results"
            )

        # shape the profiles
        for col in DATA_COLUMNS:
            # use the generation profile to shape net generation data, otherwise use the fuel profile
            if col == "net_generation_mwh":
                profile_to_use = "generation_profile"
            else:
                profile_to_use = "fuel_profile"

            shaped_partial_plants[col] = (
                shaped_partial_plants[col] * shaped_partial_plants[profile_to_use]
            )

        # remove the intermediate columns in-place
        shaped_partial_plants.drop(
            columns=["generation_profile", "fuel_profile", "flat_profile"],
            inplace=True,
        )

        # validate that the shaping process did not alter the data
        validation.validate_shaped_totals(
            shaped_partial_plants,
            eia_data_to_shape,
            year,
            group_keys=["plant_id_eia", "subplant_id"],
        )

    else:
        shaped_partial_plants = pd.DataFrame(
            columns=(
                ["plant_id_eia", "subplant_id", "report_date", "datetime_utc"]
                + DATA_COLUMNS
            )
        )
    return shaped_partial_plants


def shape_partial_cems_subplants(cems, eia923_allocated, year):
    """Scales CEMS subplant data for which there is partial units reporting.

    Returns:
        cems_data: cems dataframe with partial_cems subplant-months removed
        partial_cems_shaped: dataframe with hourly data from EIA scaled using partial
        cems data
    """
    SUBPLANT_KEYS = ["report_date", "plant_id_eia", "subplant_id"]

    # identify all of the partial cems plants and group by subplant-month
    eia_data_to_shape = eia923_allocated.loc[
        eia923_allocated.hourly_data_source == "partial_cems_subplant"
    ].copy()
    # if there is no data in the partial cems dataframe, skip.
    if len(eia_data_to_shape) > 0:
        eia_data_to_shape = (
            eia_data_to_shape.groupby(SUBPLANT_KEYS, dropna=False)[DATA_COLUMNS]
            .sum()
            .reset_index()
        )

        # split the cems data into partial cems and cems
        cems = cems.merge(
            eia_data_to_shape[SUBPLANT_KEYS],
            how="outer",
            on=SUBPLANT_KEYS,
            indicator="data_source",
            validate="m:1",
        )
        if len(cems[cems["data_source"] == "right_only"]) > 0:
            raise UserWarning(
                " At least one subplant-month identified as partial_cems does not exist in the cems data."
            )
        partial_cems_data = cems[cems["data_source"] == "both"].copy()
        partial_cems_data.drop(columns=["data_source"], inplace=True)
        cems = cems[cems["data_source"] == "left_only"]
        cems.drop(columns=["data_source"], inplace=True)

        # merge cems gross generation and fuel consumption totals into the EIA totals
        # these will be used to scale the EIA data
        partial_cems_totals = (
            partial_cems_data.groupby(SUBPLANT_KEYS, dropna=False)[
                ["gross_generation_mwh", "fuel_consumed_mmbtu"]
            ]
            .sum()
            .reset_index()
        )
        eia_data_to_shape = eia_data_to_shape.merge(
            partial_cems_totals,
            how="left",
            on=SUBPLANT_KEYS,
            validate="1:1",
            suffixes=(None, "_cems"),
        )

        partial_cems_shaped = partial_cems_data.copy()

        # shape the cems data
        for index, row in eia_data_to_shape.iterrows():
            plant_id = row.plant_id_eia
            subplant_id = row.subplant_id
            report_date = row.report_date

            for eia_column in DATA_COLUMNS:
                # we will shape net generation data based on the cems gross gen profile
                if eia_column == "net_generation_mwh":
                    cems_total_column = "gross_generation_mwh"
                    cems_column = "gross_generation_mwh"
                # all other fuel and emissions data will be shaped using the fuel profile
                else:
                    cems_total_column = "fuel_consumed_mmbtu_cems"
                    cems_column = "fuel_consumed_mmbtu"

                # if both values are zero, do nothing since the profile is already zero
                if (row[eia_column] == 0) & (row[cems_total_column] == 0):
                    partial_cems_shaped = set_value_to_zero(
                        partial_cems_shaped,
                        partial_cems_data,
                        report_date,
                        plant_id,
                        subplant_id,
                        eia_column,
                    )
                # if the eia data is positive, but the cems data is zero, use the fuel data to shape it
                elif (
                    (row[eia_column] > 0)
                    & (row[cems_total_column] == 0)
                    & (row["fuel_consumed_mmbtu_cems"] != 0)
                ):
                    cems_column = "fuel_consumed_mmbtu"
                    scaling_factor = row[eia_column] / row["fuel_consumed_mmbtu_cems"]
                    partial_cems_shaped = scale_data(
                        partial_cems_shaped,
                        partial_cems_data,
                        report_date,
                        plant_id,
                        subplant_id,
                        eia_column,
                        cems_column,
                        scaling_factor,
                    )
                # if the eia data is positive but all cems data is zero, assign a flat profile
                elif (
                    (row[eia_column] > 0)
                    & (row[cems_total_column] == 0)
                    & (row["fuel_consumed_mmbtu_cems"] == 0)
                ):
                    shift_factor = row[eia_column] - row[cems_total_column]
                    partial_cems_shaped = shift_data(
                        partial_cems_shaped,
                        partial_cems_data,
                        report_date,
                        plant_id,
                        subplant_id,
                        eia_column,
                        cems_column,
                        shift_factor,
                    )
                # if the eia value is negative (should only be for net generation), shift data
                elif row[eia_column] < 0:
                    shift_factor = row[eia_column] - row[cems_total_column]
                    partial_cems_shaped = shift_data(
                        partial_cems_shaped,
                        partial_cems_data,
                        report_date,
                        plant_id,
                        subplant_id,
                        eia_column,
                        cems_column,
                        shift_factor,
                    )
                # if the eia net generation is zero and cems gross generation is positive, shift the data
                elif (
                    (eia_column == "net_generation_mwh")
                    & (row[eia_column] == 0)
                    & (row[cems_total_column] > 0)
                ):
                    shift_factor = row[eia_column] - row[cems_total_column]
                    partial_cems_shaped = shift_data(
                        partial_cems_shaped,
                        partial_cems_data,
                        report_date,
                        plant_id,
                        subplant_id,
                        eia_column,
                        cems_column,
                        shift_factor,
                    )
                # if both values are positive, scale the data
                elif (row[eia_column] >= 0) & (row[cems_total_column] > 0):
                    scaling_factor = row[eia_column] / row[cems_total_column]
                    partial_cems_shaped = scale_data(
                        partial_cems_shaped,
                        partial_cems_data,
                        report_date,
                        plant_id,
                        subplant_id,
                        eia_column,
                        cems_column,
                        scaling_factor,
                    )
                else:
                    raise UserWarning(
                        f"Uncategorized combination of {eia_column} data for plant {plant_id} subplant {subplant_id} in {report_date}:\n   EIA data is {row[eia_column]} and CEMS data is {row[cems_total_column]}"
                    )

        # validate that the scaled totals match
        validation.validate_shaped_totals(
            partial_cems_shaped,
            eia_data_to_shape,
            year,
            group_keys=["plant_id_eia", "subplant_id"],
        )

        partial_cems_shaped.drop(
            columns=["steam_load_1000_lb", "gross_generation_mwh"],
            errors="ignore",
            inplace=True,
        )

        partial_cems_shaped = apply_dtypes(partial_cems_shaped)
    else:
        partial_cems_shaped = pd.DataFrame(
            columns=(
                ["plant_id_eia", "subplant_id", "report_date", "datetime_utc"]
                + DATA_COLUMNS
            )
        )

    return cems, partial_cems_shaped


def set_value_to_zero(
    partial_cems_shaped,
    partial_cems_data,
    report_date,
    plant_id,
    subplant_id,
    eia_column,
):
    """Used by `scale_partial_cems_data` to set a shaped value to all zeros."""
    partial_cems_shaped.loc[
        (partial_cems_shaped.report_date == report_date)
        & (partial_cems_shaped.plant_id_eia == plant_id)
        & (partial_cems_shaped.subplant_id == subplant_id),
        eia_column,
    ] = 0

    return partial_cems_shaped


def scale_data(
    partial_cems_shaped,
    partial_cems_data,
    report_date,
    plant_id,
    subplant_id,
    eia_column,
    cems_column,
    scaling_factor,
):
    """Used by `scale_partial_cems_data` to shape data using a scaling factor."""
    partial_cems_shaped.loc[
        (partial_cems_shaped.report_date == report_date)
        & (partial_cems_shaped.plant_id_eia == plant_id)
        & (partial_cems_shaped.subplant_id == subplant_id),
        eia_column,
    ] = (
        partial_cems_data.loc[
            (partial_cems_data.report_date == report_date)
            & (partial_cems_data.plant_id_eia == plant_id)
            & (partial_cems_data.subplant_id == subplant_id),
            cems_column,
        ]
        * scaling_factor
    )

    return partial_cems_shaped


def shift_data(
    partial_cems_shaped,
    partial_cems_data,
    report_date,
    plant_id,
    subplant_id,
    eia_column,
    cems_column,
    shift_factor,
):
    """Used by `scale_partial_cems_data` to shape data using a shift factor."""
    # get a count of the number of hours
    number_of_hours = len(
        partial_cems_data.loc[
            (partial_cems_data.report_date == report_date)
            & (partial_cems_data.plant_id_eia == plant_id)
            & (partial_cems_data.subplant_id == subplant_id),
            cems_column,
        ]
    )
    # divide the scaling factor by the number of hours to get the hourly shift
    hourly_shift = shift_factor / number_of_hours
    # add the shift factor to the hourly profile
    partial_cems_shaped.loc[
        (partial_cems_shaped.report_date == report_date)
        & (partial_cems_shaped.plant_id_eia == plant_id)
        & (partial_cems_shaped.subplant_id == subplant_id),
        eia_column,
    ] = (
        partial_cems_data.loc[
            (partial_cems_data.report_date == report_date)
            & (partial_cems_data.plant_id_eia == plant_id)
            & (partial_cems_data.subplant_id == subplant_id),
            cems_column,
        ]
        + hourly_shift
    )

    return partial_cems_shaped
