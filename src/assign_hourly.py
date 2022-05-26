"""
Functions for assigning hourly profiles 
"""


def distribute_monthly_eia_data_to_hourly(
    monthly_eia_data_to_distribute, hourly_profiles, profile_column_name
):
    """
    Uses monthly-level EIA data and assigns an hourly profile
    Inputs: 
        monthly_eia_data_to_distribute: a dataframe that contains monthly total net generation, 
        fuel consumption, and co2 data, along with columns for report_date and ba_code
        for a single BA 

        hourly_profiles: a dataframe that contains an hourly profile 

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
