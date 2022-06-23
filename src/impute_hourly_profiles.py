import src.load_data as load_data
from src.column_checks import apply_dtypes
import pandas as pd
import numpy as np

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


def aggregate_for_residual(
    data,
    time_key: str = "datetime_utc",
    ba_key: str = "ba_code",
    transmission: bool = False,
):
    """
        aggregate_for_residual()
    Inputs: 
        data: dataframe with time_key, ba_key, "fuel_category_eia930", "net_generation_mwh" and "distribution_flag" columns 

    Utility function for trying different BA aggregations in 930 and 923 data
    """
    if transmission:
        data = data[data["distribution_flag"] == False]

    data = (
        data.groupby([ba_key, "fuel_category_eia930", time_key], dropna=False)[
            "net_generation_mwh"
        ]
        .sum()
        .reset_index()
        .rename(columns={"ba_code_physical": "ba_code",})
    )

    return data


def calculate_residual(cems, eia930_data, plant_attributes, year: int):
    """
        calculate_residual
    Inputs: 
        `cems`: dataframe of CEMS hourly plant-level data containing columns
            `plant_id_eia`
        `eia930`: dataframe of 930 hourly BA-level data containing columns 
            `net_generation_mwh_930`
        `plant_frame` dataframe of static plant data containing columns 
            `plant_id_eia`, `ba_code`, `ba_code_physical`
    Returns: 
        Dataframe of hourly profiles, containing columns

    Scaling: 
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
    """
    # Options for how to group. Could make command line arguments if needed.
    # transmission = True and physical BA code is based on EIA-930 instructions
    TRANSMISSION = False  # use only transmission-level connections?
    BA_CODE = "ba_code"  # ba_code or ba_code_physical?

    # Name column same as 930, hourly_profiles.
    cems = cems.merge(plant_attributes, how="left", on="plant_id_eia")

    cems = aggregate_for_residual(
        cems, "datetime_utc", BA_CODE, transmission=TRANSMISSION
    )
    combined_data = eia930_data.merge(
        cems, how="left", on=["ba_code", "fuel_category_eia930", "datetime_utc"]
    )
    # only keep rows where local datetime is in the current year
    combined_data = combined_data[
        combined_data["datetime_local"].apply(lambda x: x.year) == year
    ]

    # if there is no cems data for a ba-fuel, replace missing values with zero
    combined_data["net_generation_mwh"] = combined_data["net_generation_mwh"].fillna(0)

    # Find scaling factor
    # only keep data where the cems data is greater than zero
    scaling_factors = combined_data.copy()[combined_data["net_generation_mwh"] != 0]
    # calculate the ratio of 930 net generation to cems net generation
    # if correct, ratio should be >=1
    scaling_factors["scaling_factor"] = (
        scaling_factors["net_generation_mwh_930"]
        / scaling_factors["net_generation_mwh"]
    )
    # find the minimum ratio for each ba-fuel
    scaling_factors = (
        scaling_factors.groupby(["ba_code", "fuel_category_eia930"], dropna=False)[
            "scaling_factor"
        ]
        .min()
        .reset_index()
    )

    # only keep scaling factors < 1, which means the data needs to be scaled
    scaling_factors = scaling_factors[scaling_factors["scaling_factor"] < 1]

    # merge the scaling factor into the combined data
    # for any BA-fuels without a scaling factor, fill with 1 (scale to 100% of the origina data)
    combined_data = combined_data.merge(
        scaling_factors, how="left", on=["ba_code", "fuel_category_eia930"]
    ).fillna(1)

    # calculate the scaled cems data
    combined_data["cems_scaled"] = (
        combined_data["net_generation_mwh"] * combined_data["scaling_factor"]
    )

    # calculate the residual
    combined_data["profile"] = (
        combined_data["net_generation_mwh_930"] - combined_data["cems_scaled"]
    )

    # identify the method used to calculate the profile
    # if the scaling factor is 1, then the profile was not scaled
    combined_data = combined_data.assign(
        profile_method=lambda x: np.where(
            (x.scaling_factor == 1), "residual", "scaled_residual"
        )
    )

    return combined_data[
        [
            "ba_code",
            "fuel_category_eia930",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "profile",
            "profile_method",
        ]
    ]


def create_flat_profile(report_date, ba, fuel):
    year = report_date.year

    df_temporary = pd.DataFrame(
        index=pd.date_range(
            start=f"{year-1}-12-31 00:00:00",
            end=f"{year+1}-01-01 23:00:00",
            freq="H",
            tz="UTC",
            name="datetime_utc",
        ),
        columns=["ba_code", "fuel_category"],
    ).reset_index()

    df_temporary["profile"] = 1.0
    df_temporary["datetime_local"] = df_temporary["datetime_utc"]
    df_temporary["datetime_local"] = df_temporary["datetime_utc"].dt.tz_convert(
        load_data.ba_timezone(ba=ba, type="local")
    )
    # only keep data for which the local datetime is in the current year
    df_temporary = df_temporary[df_temporary["datetime_local"].dt.year == year]

    # create a report date column
    df_temporary["report_date"] = df_temporary["datetime_local"].astype(str).str[:7]
    df_temporary["report_date"] = pd.to_datetime(df_temporary["report_date"])

    # only keep the report dates that match
    df_temporary = df_temporary[df_temporary["report_date"] == report_date]

    df_temporary["ba_code"] = ba
    df_temporary["fuel_category"] = fuel

    return df_temporary


def impute_missing_hourly_profiles(
    monthly_eia_data_to_shape, residual_profiles, plant_attributes, year
):
    """Identify and estimate hourly profiles for missing BA-fuels."""
    # change local datetime column to string due to challenges with mixed timezones
    residual_profiles["datetime_local"] = residual_profiles["datetime_local"].astype(
        str
    )

    residual_profiles = residual_profiles.rename(
        columns={"fuel_category_eia930": "fuel_category"}
    )

    # round the data to the nearest tenth
    residual_profiles["profile"] = residual_profiles["profile"].round(1)

    missing_profiles, residual_profiles = identify_missing_profiles(
        monthly_eia_data_to_shape, residual_profiles, plant_attributes
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
                    (dibas.ba_code == ba)
                    & (dibas.ba_region == dibas.diba_region)
                    & (dibas.timezone_local == dibas.timezone_local_diba),
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
            df_temporary["profile_method"] = "assumed_flat"

        # For now assume hydro is dispatched with a flat profile
        # TODO improve this assumption
        elif fuel in ["hydro"]:
            df_temporary = create_flat_profile(report_date, ba, fuel)
            df_temporary["profile_method"] = "assumed_flat"
        # for any other fossil resources, use a flat profile
        # NOTE: we need to improve this method
        elif fuel in ["natural_gas", "coal", "petroleum"]:
            df_temporary = create_flat_profile(report_date, ba, fuel)
            df_temporary["profile_method"] = "assumed_flat"
        elif fuel in ["other"]:
            # assign a flat profile
            df_temporary = create_flat_profile(report_date, ba, fuel)
            df_temporary["profile_method"] = "assumed_flat"
        else:
            raise UserWarning(f"Fuel category {fuel} not recognized.")

        hourly_profiles_to_add.append(df_temporary)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )

    hourly_profiles = pd.concat([residual_profiles, hourly_profiles_to_add], axis=0)

    # round the data to the nearest tenth
    hourly_profiles["profile"] = hourly_profiles["profile"].round(1)

    print("Summary of methods used to estimate missing hourly profiles:")
    print(
        hourly_profiles[["ba_code", "fuel_category", "report_date", "profile_method"]]
        .drop_duplicates()
        .drop(columns=["ba_code"])
        .pivot_table(index="fuel_category", columns="profile_method", aggfunc="count")
        .fillna(0)
        .astype(int)
    )

    return hourly_profiles


def identify_missing_profiles(
    monthly_eia_data_to_shape, residual_profiles, plant_attributes
):
    # drop ba fuel months where the reported profile is all zeros
    MONTHLY_GROUP_COLUMNS = [
        "ba_code",
        "fuel_category",
        "report_date",
    ]
    ba_fuel_months_with_zero_generation = (
        residual_profiles.groupby(MONTHLY_GROUP_COLUMNS, dropna=False)
        .sum()
        .reset_index()
    )
    # identify the BA fuel months that have zero reported generation
    ba_fuel_months_with_zero_generation = ba_fuel_months_with_zero_generation.loc[
        ba_fuel_months_with_zero_generation["profile"] == 0, MONTHLY_GROUP_COLUMNS
    ]

    # remove residual profiles that are all zeros
    residual_profiles = residual_profiles.merge(
        ba_fuel_months_with_zero_generation,
        how="outer",
        on=MONTHLY_GROUP_COLUMNS,
        indicator="source",
    )
    residual_profiles = residual_profiles[residual_profiles["source"] == "left_only"]
    residual_profiles = residual_profiles.drop(columns="source")

    # determine for which BA-fuels we are missing residual profiles
    available_profiles = residual_profiles[MONTHLY_GROUP_COLUMNS].drop_duplicates()
    # add ba and fuel data to the plant-level monthly data
    monthly_eia_data_to_shape = monthly_eia_data_to_shape.merge(
        plant_attributes[["plant_id_eia", "fuel_category", "ba_code"]],
        how="left",
        on="plant_id_eia",
    )
    ba_fuel_to_distribute = monthly_eia_data_to_shape[
        MONTHLY_GROUP_COLUMNS
    ].drop_duplicates()
    missing_profiles = ba_fuel_to_distribute.merge(
        available_profiles, how="outer", on=MONTHLY_GROUP_COLUMNS, indicator="source",
    )
    # identify ba fuel months where there is no data in the available residual profiles
    missing_profiles = missing_profiles[missing_profiles.source == "left_only"]
    missing_profiles = missing_profiles.sort_values(by=MONTHLY_GROUP_COLUMNS)

    return missing_profiles, residual_profiles


def average_diba_wind_solar_profiles(
    residual_profiles, ba, fuel, report_date, ba_dibas
):

    # calculate the average generation profile for the fuel in all neighboring DIBAs
    df_temporary = residual_profiles.copy()[
        (residual_profiles["ba_code"].isin(ba_dibas))
        & (residual_profiles["fuel_category"] == fuel)
        & (residual_profiles["report_date"] == report_date)
    ]
    if len(df_temporary) == 0:
        # if this error is raised, we might have to implement an approach that uses average values for the wider region
        raise UserWarning(
            f"There is no {fuel} data for the balancing authorities interconnected to {ba}"
        )
    else:
        df_temporary = (
            df_temporary.groupby(
                ["fuel_category", "datetime_utc", "datetime_local", "report_date",],
                dropna=False,
            )
            .mean()
            .reset_index()
        )
        df_temporary["ba_code"] = ba
        df_temporary["profile_method"] = "DIBA_average"

    return df_temporary


def average_national_wind_solar_profiles(residual_profiles, ba, fuel, report_date):
    df_temporary = residual_profiles.copy()[
        (residual_profiles["fuel_category"] == fuel)
        & (residual_profiles["report_date"] == report_date)
    ]
    # strip the time zone information so we can group by local time
    df_temporary["datetime_local"] = df_temporary["datetime_local"].str[:-6]
    df_temporary = (
        df_temporary.groupby(
            ["fuel_category", "datetime_local", "report_date",], dropna=False,
        )
        .mean()
        .reset_index()
    )
    df_temporary["ba_code"] = ba
    df_temporary["profile_method"] = "national_average"

    # re-localize the datetime_local
    local_tz = load_data.ba_timezone(ba, "local")
    df_temporary["datetime_local"] = pd.to_datetime(df_temporary["datetime_local"])
    df_temporary["datetime_local"] = (
        df_temporary["datetime_local"]
        .dt.tz_localize(local_tz, nonexistent="NaT", ambiguous="NaT")
        .fillna(method="ffill")
    )
    df_temporary["datetime_utc"] = df_temporary["datetime_local"].dt.tz_convert("UTC")

    return df_temporary


def convert_profile_to_percent(hourly_profiles):
    """converts hourly timeseries profiles from absolute mwh to percentage of monthly total mwh."""
    # convert the profile so that each hour is a percent of the monthly total
    MONTHLY_GROUP_COLUMNS = [
        "ba_code",
        "fuel_category",
        "report_date",
        "profile_method",
    ]

    monthly_total = (
        hourly_profiles.groupby(MONTHLY_GROUP_COLUMNS, dropna=False).sum().reset_index()
    )

    hourly_profiles = hourly_profiles.merge(
        monthly_total,
        how="left",
        on=MONTHLY_GROUP_COLUMNS,
        suffixes=(None, "_monthly_total"),
    )
    hourly_profiles["profile"] = (
        hourly_profiles["profile"] / hourly_profiles["profile_monthly_total"]
    )
    hourly_profiles = hourly_profiles.drop(columns="profile_monthly_total")

    return hourly_profiles


def get_synthetic_plant_id_from_ba_fuel(df):
    """
        Return artificial plant code. Max real plant is 64663
        Our codes look like 9BBBFF where BBB is the three digit BA number and FF is the 
        two-digit fuel number

        df must contain `ba_code` and `fuel_category`
    """

    # load the ba reference table with all of the ba number ids
    ba_numbers = pd.read_csv("../data/manual/ba_reference.csv")[
        ["ba_code", "ba_number"]
    ]
    # reformat the number with leading zeros
    ba_numbers["ba_number"] = ba_numbers["ba_number"].astype(str).str.zfill(3)
    # convert to a dictionary
    ba_numbers = dict(zip(ba_numbers["ba_code"], ba_numbers["ba_number"]))

    # make sure the ba codes are strings
    df["ba_code"] = df["ba_code"].astype(str)
    # create a new column with the synthetic plant ids
    df["plant_id_eia"] = df.apply(
        lambda row: f"9{ba_numbers[row['ba_code']]}{FUEL_NUMBERS[row['fuel_category']]}",
        axis=1,
    )
    # convert to an int32 column
    df["plant_id_eia"] = df["plant_id_eia"].astype("Int32")

    return df


def aggregate_eia_data_to_ba_fuel(monthly_eia_data_to_shape, plant_attributes):
    """
        Given cleaned monthly EIA-923 data and plant attributes, aggregate to BA-fuel
        using artificial plant IDs 9XXXYYY where XXX=BA code (see `ba_reference.csv`)
        and YY=fuel (see `impute_hourly_profiles.get_synthetic_plant_id_from_ba_fuel`)

        Add new artificial plants to plant_attributes frame.
    """

    # Note: currently using ba_code, could alternatively use ba_code_physical
    # Add plant attributes for grouping
    eia_agg = monthly_eia_data_to_shape.merge(
        plant_attributes[["plant_id_eia", "ba_code", "fuel_category"]],
        how="left",
        on="plant_id_eia",
    )

    # Group
    eia_agg = (
        eia_agg.groupby(["ba_code", "report_date", "fuel_category"], dropna=False)
        .sum()
        .reset_index()
        .drop(columns=["plant_id_eia", "subplant_id"])
    )

    # Make nan BA "None" so equality test works
    eia_agg.ba_code = eia_agg.ba_code.replace(np.nan, None)

    eia_agg = get_synthetic_plant_id_from_ba_fuel(eia_agg)

    # Add BA code and fuel type for synthetic plants into plant_attributes.
    # Note: We leave nan/NA values for all other columns, as the plants composing each synthetic plant may have a combination of values
    to_add = (
        eia_agg.groupby("plant_id_eia")
        .first()
        .reset_index()[["plant_id_eia", "ba_code", "fuel_category"]]
    )
    plant_attributes = pd.concat([plant_attributes, to_add])

    return eia_agg, plant_attributes


def shape_monthly_eia_data_as_hourly(monthly_eia_data_to_shape, hourly_profiles):
    """
    Uses monthly-level EIA data and assigns an hourly profile
    Intended for calling after `monthly_eia_data_to_ba` 
    Inputs: 
        shaped_monthly_data: a dataframe that contains monthly total net generation, 
            fuel consumption, and co2 data, along with columns for report_date and ba_code
    """
    # specify columns containing monthly data that should be distributed to hourly
    DATA_COLUMNS = [
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
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
    ]

    # merge the hourly profiles into each plant-month
    shaped_monthly_data = monthly_eia_data_to_shape.merge(
        hourly_profiles, how="left", on=["report_date", "fuel_category", "ba_code"]
    )

    # plant-months where there is negative net generation, assign a flat profile
    # to do this, we assign each hour an equal share of the total number of hours in the month
    shaped_monthly_data.loc[
        shaped_monthly_data["net_generation_mwh"] < 0, "profile"
    ] = 1 / (shaped_monthly_data["report_date"].dt.daysinmonth * 24)
    # update the method column
    shaped_monthly_data.loc[
        shaped_monthly_data["net_generation_mwh"] < 0, "profile_method"
    ] = "flat_negative_generation"

    # shape the data
    for column in DATA_COLUMNS:
        shaped_monthly_data[column] = (
            shaped_monthly_data[column] * shaped_monthly_data["profile"]
        )

    # re order the columns
    column_order = [
        "plant_id_eia",
        "ba_code",
        "fuel_category",
        "datetime_utc",
        "report_date",
        "profile_method",
    ] + DATA_COLUMNS

    # re-order and drop intermediate columns
    shaped_monthly_data = shaped_monthly_data[column_order]

    return shaped_monthly_data


def scale_partial_cems_data(cems, eia923_allocated):
    """Scales CEMS subplant data for which there is partial units reporting"""
    subplant_keys = ["report_date", "plant_id_eia", "subplant_id"]
    data_columns = [
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "net_generation_mwh",
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
        "co2_mass_lb_adjusted",
        "ch4_mass_lb_adjusted",
        "n2o_mass_lb_adjusted",
        "nox_mass_lb_adjusted",
        "so2_mass_lb_adjusted",
    ]

    # identify all of the partial cems plants and group by subplant-month
    partial_cems = eia923_allocated.loc[
        eia923_allocated.hourly_data_source == "partial_cems"
    ]
    partial_cems = (
        partial_cems.groupby(subplant_keys, dropna=False)
        .sum()[data_columns]
        .reset_index()
    )

    # group the cems data by subplant month and merge into partial cems
    cems_monthly = (
        cems.groupby(subplant_keys, dropna=False).sum()[data_columns].reset_index()
    )
    partial_cems = partial_cems.merge(
        cems_monthly, how="left", on=subplant_keys, suffixes=("_eia", "_cems")
    )

    # compare the fuel consumption values from each source. If the CEMS data
    # actually represents only partial data, the CEMS-reported values should be
    # less than the EIA values. Where the CEMS-reported values are greater than
    # the EIA-reported values, change the flag back to cems only, and remove
    # these from the partial cems list
    complete_cems = partial_cems.loc[
        (partial_cems.fuel_consumed_mmbtu_cems >= partial_cems.fuel_consumed_mmbtu_eia),
        subplant_keys,
    ]
    eia923_allocated = eia923_allocated.merge(
        complete_cems, how="outer", on=subplant_keys, indicator="source"
    )
    eia923_allocated.loc[
        eia923_allocated["source"] == "both", "hourly_data_source"
    ] = "cems"
    eia923_allocated = eia923_allocated.drop(columns="source")
    partial_cems = partial_cems[
        (partial_cems.fuel_consumed_mmbtu_cems < partial_cems.fuel_consumed_mmbtu_eia)
    ]

    # create a version of the cems data that is aggregated at the subplant level
    #  and filtered to include only the subplant-months that need to be scaled
    cems_scaled = cems.merge(
        partial_cems[subplant_keys], how="outer", on=subplant_keys, indicator="source"
    )
    cems_scaled = cems_scaled[cems_scaled["source"] == "both"].drop(columns=["source"])
    cems_scaled = (
        cems_scaled.groupby(subplant_keys + ["datetime_utc"], dropna=False)
        .sum()[data_columns]
        .reset_index()
    )

    # scale the cems data
    for index, row in partial_cems.iterrows():
        plant_id = row.plant_id_eia
        subplant_id = row.subplant_id
        report_date = row.report_date

        for column in data_columns:
            try:
                scaling_factor = row[f"{column}_eia"] / row[f"{column}_cems"]
                scaling_method = "multiply_by_cems_value"

                if scaling_factor < 0:
                    scaling_factor = row[f"{column}_eia"] - row[f"{column}_cems"]
                    scaling_method = "shift_negative_profile"

            except ZeroDivisionError:
                if (row[f"{column}_eia"] == 0) & (row[f"{column}_cems"] == 0):
                    scaling_factor = 1
                    scaling_method = "multiply_by_cems_value"
                elif (row[f"{column}_eia"] > 0) & (row[f"{column}_cems"] == 0):
                    scaling_factor = (
                        row[f"{column}_eia"] / row["fuel_consumed_mmbtu_cems"]
                    )
                    scaling_method = "multiply_by_cems_fuel"
                elif (row[f"{column}_eia"] < 0) & (row[f"{column}_cems"] == 0):
                    scaling_factor = row[f"{column}_eia"] - row[f"{column}_cems"]
                    scaling_method = "shift_negative_profile"

            if scaling_method == "multiply_by_cems_value":
                cems_scaled.loc[
                    (cems_scaled.report_date == report_date)
                    & (cems_scaled.plant_id_eia == plant_id)
                    & (cems_scaled.subplant_id == subplant_id),
                    column,
                ] = (
                    cems_scaled.loc[
                        (cems_scaled.report_date == report_date)
                        & (cems_scaled.plant_id_eia == plant_id)
                        & (cems_scaled.subplant_id == subplant_id),
                        column,
                    ]
                    * scaling_factor
                )
            elif scaling_method == "multiply_by_cems_fuel":
                cems_scaled.loc[
                    (cems_scaled.report_date == report_date)
                    & (cems_scaled.plant_id_eia == plant_id)
                    & (cems_scaled.subplant_id == subplant_id),
                    column,
                ] = (
                    cems_scaled.loc[
                        (cems_scaled.report_date == report_date)
                        & (cems_scaled.plant_id_eia == plant_id)
                        & (cems_scaled.subplant_id == subplant_id),
                        "fuel_consumed_mmbtu",
                    ]
                    * scaling_factor
                )
            elif scaling_method == "shift_negative_profile":
                # get a count of the number of hours
                number_of_hours = len(
                    cems_scaled.loc[
                        (cems_scaled.report_date == report_date)
                        & (cems_scaled.plant_id_eia == plant_id)
                        & (cems_scaled.subplant_id == subplant_id),
                        column,
                    ]
                )
                # divide the scaling factor by the number of hours to get the hourly shift
                hourly_shift = scaling_factor / number_of_hours
                # add the shift factor to the hourly profile
                cems_scaled.loc[
                    (cems_scaled.report_date == report_date)
                    & (cems_scaled.plant_id_eia == plant_id)
                    & (cems_scaled.subplant_id == subplant_id),
                    column,
                ] = (
                    cems_scaled.loc[
                        (cems_scaled.report_date == report_date)
                        & (cems_scaled.plant_id_eia == plant_id)
                        & (cems_scaled.subplant_id == subplant_id),
                        "fuel_consumed_mmbtu",
                    ]
                    + hourly_shift
                )

    cems_scaled = apply_dtypes(cems_scaled)

    return cems_scaled, eia923_allocated
