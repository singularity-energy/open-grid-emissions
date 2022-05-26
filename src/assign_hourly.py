"""
Functions for assigning hourly profiles 
"""

import pandas as pd


# TODO: moved from data_pipeline. Refactor data_pipeline so it calls this version.
def create_flat_profile(year):
    df_temp = pd.DataFrame(
        index=pd.date_range(
            start=f"{year-1}-12-31 00:00:00",
            end=f"{year+1}-01-01 23:00:00",
            freq="H",
            tz="UTC",
            name="datetime_utc",
        ),
        columns=["ba_code", "fuel_category"],
    ).reset_index()

    df_temp["net_generation_mwh_930"] = 1.0
    df_temp["datetime_local"] = df_temp["datetime_utc"]
    df_temp["datetime_local"] = df_temp["datetime_utc"].dt.tz_convert(
        data_cleaning.ba_timezone(ba=ba, type="local")
    )
    # create a report date column
    df_temp["report_date"] = df_temp["datetime_local"].astype(str).str[:7]
    df_temp["report_date"] = pd.to_datetime(df_temp["report_date"])

    return df_temp


def add_missing_profiles(
    monthly_eia_data_to_distribute, hourly_profiles, profile_column_name: str
):
    """
    Adds default 930 profiles for missing BA/fuel combos in EIA-923 data.
    TODO: may be better assumption than flat for some profiles
    Note: this is not directly generating default residual profiles.

    Inputs: 
        monthly_eia_data_to_distribute: a dataframe that contains fuel_category,
        report_date and ba_code for a single BA. (other columns allowed but not required)

        hourly_profiles

    """
    # for fuel categories that exist in the EIA-923 data but not in EIA-930, create flat profiles to add to the hourly profiles from 930
    # TODO: Identify for which BA-fuels a flat profile was created
    # TODO: Is there a better assumption than flat?
    ba_list = list(monthly_eia_data_to_distribute["ba_code"].dropna().unique())

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    # for each ba
    for ba in ba_list:
        # get a list of fuels categories that exist in that BA
        ba_fuel_list = list(
            monthly_eia_data_to_distribute.loc[
                monthly_eia_data_to_distribute["ba_code"] == ba, "fuel_category"
            ].unique()
        )
        for fuel in ba_fuel_list:
            # if there is no data for that fuel type in the eia930 data, create a flat profile
            if (
                len(
                    hourly_profiles[
                        (hourly_profiles["ba_code"] == ba)
                        & (hourly_profiles["fuel_category"] == fuel)
                    ]
                )
                == 0
            ):
                print(f"Adding flat profile for {ba} {fuel}")
                df_temp = create_flat_profile(year)
                df_temp["ba_code"] = ba
                df_temp["fuel_category"] = fuel
                hourly_profiles_to_add.append(df_temp)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )
    # concat the flat profiles to the hourly profiles
    return pd.concat([hourly_profiles, hourly_profiles_to_add], axis=0)


def distribute_monthly_eia_data_to_hourly(
    monthly_eia_data_to_distribute, hourly_profiles, profile_column_name
):
    """
    Uses monthly-level EIA data and assigns an hourly profile
    Inputs: 
        monthly_eia_data_to_distribute: a dataframe that contains monthly total net generation, 
        fuel consumption, and co2 data, along with columns for report_date and ba_code
        for a single BA 

        hourly_profiles: a dataframe that contains hourly residual profiles for all fuel types 

        profile_column_name: name of that profile in hourly_profiles
    """
    columns_to_shape = [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_tons",
        "co2_mass_tons_adjusted",
    ]

    # calculate totals by BA, Fuel Group, and Month
    monthly_eia_ba_fuel = (
        monthly_eia_data_to_distribute.groupby(
            ["ba_code", "fuel_category", "report_date"]
        )
        .sum()[columns_to_shape]
        .reset_index()
    )

    # calculate the total monthly net generation profile by BA and fuel group
    monthly_profile_total = (
        hourly_profiles.groupby(["ba_code", "fuel_category", "report_date"])
        .sum()
        .reset_index()
    )

    # merge the total monthly profile into the monthly totals
    monthly_eia_ba_fuel = monthly_eia_ba_fuel.merge(
        monthly_profile_total,
        how="left",
        on=["ba_code", "fuel_category", "report_date"],
    )

    # calculate how much net generation, fuel, and co2 should be assigned to each unit of net generation in the profile
    for col in columns_to_shape:
        monthly_eia_ba_fuel[col] = (
            monthly_eia_ba_fuel[col] / monthly_eia_ba_fuel[profile_column_name]
        )

    # drop the profile column and merge the hourly generation, fuel, and co2 factors back into the profile timeseries data
    monthly_eia_ba_fuel = monthly_eia_ba_fuel.drop(columns=profile_column_name)
    hourly_eia_data = hourly_profiles.merge(
        monthly_eia_ba_fuel, how="left", on=["ba_code", "fuel_category", "report_date"]
    )

    # multiply each factor by the profile to calculate the hourly shape
    for col in columns_to_shape:
        hourly_eia_data[col] = (
            hourly_eia_data[col] * hourly_eia_data[profile_column_name]
        )

    # create a column identifying the source of the data
    hourly_eia_data["data_source"] = "EIA"

    return hourly_eia_data
