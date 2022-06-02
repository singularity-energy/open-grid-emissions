import src.data_cleaning as data_cleaning
import src.load_data as load_data
import pandas as pd


def create_flat_profile(year, ba, fuel):
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

    df_temp["ba_code"] = ba
    df_temp["fuel_category"] = fuel

    return df_temp

def load_hourly_profiles(monthly_eia_data_to_distribute, year):
    # load the residual hourly profiles
    residual_profiles = pd.read_csv(
        "../data/outputs/residual_profiles.csv", parse_dates=["report_date"]
    )

    # determine for which BA-fuels we are missing residual profiles
    available_profiles = residual_profiles[['ba_code','fuel_category']].drop_duplicates()
    ba_fuel_to_distribute = monthly_eia_data_to_distribute[['ba_code','fuel_category']].drop_duplicates().dropna()
    missing_profiles = ba_fuel_to_distribute.merge(available_profiles, how='outer', on=['ba_code','fuel_category'], indicator='source')
    missing_profiles = missing_profiles[missing_profiles.source == 'left_only']
    missing_profiles.sort_values(by=['fuel_category','ba_code'])

    # load information about directly interconnected balancing authorities (DIBAs)
    # this will help us fill profiles using data from nearby BAs
    dibas = load_data.load_diba_data(year)

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    for index, row in missing_profiles.iterrows():
        ba = row['ba_code']
        fuel = row['fuel_category']

        # if geothermal or nuclear, assign a flat profile as baseload
        if fuel in ['geoethermal','nuclear']:
            print(f"Adding flat baseload profile for {ba} {fuel}")
            df_temp = create_flat_profile(year, ba, fuel)
        # use the profile from "other" if available
        elif fuel in ['biomass','waste']:
            if len(residual_profiles[
                (residual_profiles["ba_code"] == ba)
                & (residual_profiles["fuel_category"] == 'other')
            ]) >= 8760:
                print(f"Adding profile for {ba} {fuel} based on `other` fuel profile")
                df_temp = residual_profiles[
                    (residual_profiles["ba_code"] == ba)
                    & (residual_profiles["fuel_category"] == 'other')
                ]
            else:
                # assign a flat profile
                print(f"Adding flat baseload profile for {ba} {fuel}")
                df_temp = create_flat_profile(year, ba, fuel)
        elif fuel in ['wind','solar']:
            # get a list of diba located in the same region and located in the same time zone
            ba_dibas = list(dibas.loc[(dibas.ba_code == ba) & (dibas.ba_region == dibas.diba_region) & (dibas.timezone_local == dibas.timezone_local_diba), 'diba_code'].unique())
            if len(ba_dibas) > 0:
                # calculate the average generation profile for the fuel in all neighboring DIBAs
                df_temp = residual_profiles[
                                    (residual_profiles["ba_code"].isin(ba_dibas))
                                    & (residual_profiles["fuel_category"] == fuel)
                                ]
                if len(df_temp) == 0:
                    # if this error is raised, we might have to implement an approach that uses average values for the wider region
                    raise UserWarning(f'There is no {fuel} data for the balancing authorities interconnected to {ba}')
                else:
                    df_temp = df_temp.groupby(['fuel_category','datetime_utc','datetime_local','report_date',]).mean().reset_index()
                    # check that the length is less than 8784
                    if len(df_temp) > 8784:
                        raise UserWarning(f'Length of {fuel} profile is {len(df_temp)}, expected 8760 or 8784. Check that local timezones of DIBAs are the same as {ba}')
                    df_temp['ba_code'] = ba
            else:
                raise UserWarning(f'There are no balancing authorities directly interconnected to {ba} in the same time zone')

        hourly_profiles_to_add.append(df_temp)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )

    hourly_profiles = pd.concat([residual_profiles, hourly_profiles_to_add], axis=0)

    return hourly_profiles



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
                df_temp = create_flat_profile(year, ba, fuel)
                hourly_profiles_to_add.append(df_temp)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )

    return pd.concat([hourly_profiles, hourly_profiles_to_add], axis=0)

