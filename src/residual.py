import src.data_cleaning as data_cleaning
import pandas as pd


def create_flat_profile(year, ba):
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


def assign_flat_profiles(monthly_eia_data_to_distribute, hourly_profiles, year):
    """
     for fuel categories that exist in the EIA-923 data but not in EIA-930, 
     create flat profiles to add to the hourly profiles from 930
     TODO: Identify for which BA-fuels a flat profile was created
     TODO: Is there a better assumption than flat?
    """
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
                df_temp = create_flat_profile(year, ba)
                df_temp["ba_code"] = ba
                df_temp["fuel_category"] = fuel
                hourly_profiles_to_add.append(df_temp)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )

    return pd.concat([hourly_profiles, hourly_profiles_to_add], axis=0)

