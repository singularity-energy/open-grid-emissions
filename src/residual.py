import src.data_cleaning as data_cleaning
import src.load_data as load_data
import pandas as pd


def create_flat_profile(year, ba, fuel):
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
        data_cleaning.ba_timezone(ba=ba, type="local")
    )
    # create a report date column
    df_temporary["report_date"] = df_temporary["datetime_local"].astype(str).str[:7]
    df_temporary["report_date"] = pd.to_datetime(df_temporary["report_date"])

    df_temporary["ba_code"] = ba
    df_temporary["fuel_category"] = fuel

    return df_temporary


def load_hourly_profiles(monthly_eia_data_to_distribute, year):
    # load the residual hourly profiles
    residual_profiles = pd.read_csv(
        "../data/outputs/residual_profiles.csv", parse_dates=["report_date"]
    ).rename(columns={"residual_scaled": "profile"})
    residual_profiles["profile_method"] = "residual"

    # determine for which BA-fuels we are missing residual profiles
    available_profiles = residual_profiles[
        ["ba_code", "fuel_category"]
    ].drop_duplicates()
    ba_fuel_to_distribute = (
        monthly_eia_data_to_distribute[["ba_code", "fuel_category"]]
        .drop_duplicates()
        .dropna()
    )
    missing_profiles = ba_fuel_to_distribute.merge(
        available_profiles,
        how="outer",
        on=["ba_code", "fuel_category"],
        indicator="source",
    )
    missing_profiles = missing_profiles[missing_profiles.source == "left_only"]
    missing_profiles.sort_values(by=["fuel_category", "ba_code"])

    # load information about directly interconnected balancing authorities (DIBAs)
    # this will help us fill profiles using data from nearby BAs
    dibas = load_data.load_diba_data(year)

    # create an hourly datetime series in local time for each ba/fuel type
    hourly_profiles_to_add = []

    for index, row in missing_profiles.iterrows():
        ba = row["ba_code"]
        fuel = row["fuel_category"]

        # For fuels that should be reported in 930, but are missing, assign a flat profile
        if fuel in ["nuclear", "hydro"]:
            df_temporary = create_flat_profile(year, ba, fuel)
            df_temporary["profile_method"] = "assumed_flat"
        # for fuels that would be grouped into the "other" category in EIA-930,
        # use the other profile if available. Otherwise, set a flat profile
        elif fuel in ["geothermal", "biomass", "waste"]:
            if (
                len(
                    residual_profiles[
                        (residual_profiles["ba_code"] == ba)
                        & (residual_profiles["fuel_category"] == "other")
                    ]
                )
                >= 8760
            ):
                df_temporary = residual_profiles.copy()[
                    (residual_profiles["ba_code"] == ba)
                    & (residual_profiles["fuel_category"] == "other")
                ]
                df_temporary["fuel_category"] = fuel
                df_temporary["profile_method"] = "other_profile"
            else:
                # assign a flat profile
                df_temporary = create_flat_profile(year, ba, fuel)
                df_temporary["profile_method"] = "assumed_flat"
        # for wind and solar, average the wind and solar generation profiles from
        # nearby interconnected BAs
        elif fuel in ["wind", "solar"]:
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
                # calculate the average generation profile for the fuel in all neighboring DIBAs
                df_temporary = residual_profiles.copy()[
                    (residual_profiles["ba_code"].isin(ba_dibas))
                    & (residual_profiles["fuel_category"] == fuel)
                ]
                if len(df_temporary) == 0:
                    # if this error is raised, we might have to implement an approach that uses average values for the wider region
                    raise UserWarning(
                        f"There is no {fuel} data for the balancing authorities interconnected to {ba}"
                    )
                else:
                    df_temporary = (
                        df_temporary.groupby(
                            [
                                "fuel_category",
                                "datetime_utc",
                                "datetime_local",
                                "report_date",
                            ]
                        )
                        .mean()
                        .reset_index()
                    )
                    # check that the length is less than 8784
                    if len(df_temporary) > 8784:
                        raise UserWarning(
                            f"Length of {fuel} profile is {len(df_temporary)}, expected 8760 or 8784. Check that local timezones of DIBAs are the same as {ba}"
                        )
                    df_temporary["ba_code"] = ba
                    df_temporary["profile_method"] = "DIBA_average"
        # for any other resources, use a flat profile
        # NOTE: all other resources should be available using other methods first
        else:
            df_temporary = create_flat_profile(year, ba, fuel)
            df_temporary["profile_method"] = "assumed_flat"

        hourly_profiles_to_add.append(df_temporary)

    hourly_profiles_to_add = pd.concat(
        hourly_profiles_to_add, axis=0, ignore_index=True
    )

    hourly_profiles = pd.concat([residual_profiles, hourly_profiles_to_add], axis=0)

    print("Summary of methods used to estimate missing hourly profiles:")
    print(
        hourly_profiles[["ba_code", "fuel_category", "profile_method"]]
        .drop_duplicates()
        .pivot_table(index="fuel_category", columns="profile_method", aggfunc="count")
        .fillna(0)
        .astype(int)
    )

    return hourly_profiles

