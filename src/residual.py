import src.data_cleaning as data_cleaning
import src.load_data as load_data
import pandas as pd



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
        data.groupby([ba_key, "fuel_category_eia930", time_key])["net_generation_mwh"]
        .sum()
        .reset_index()
        .rename(
            columns={
                "fuel_category_eia930": "fuel_category",
                "ba_code_physical": "ba_code",
            }
        )
    )

    return data


def calculate_residual(cems, eia930, plant_attributes, year: int):
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
    TRANSMISSION = True  # use only transmission-level connections?
    BA_CODE = "ba_code_physical"  # ba_code or ba_code_physical?

    # Name column same as 930, hourly_profiles. TODO: Choose just one datetime column name
    cems = cems.rename(columns={"operating_datetime_utc": "datetime_utc"})
    cems = cems.merge(
        plant_attributes, how="left", on="plant_id_eia", suffixes=("_orig", "")
    )

    cems = aggregate_for_residual(
        cems, "datetime_utc", BA_CODE, transmission=TRANSMISSION
    )
    combined_data = eia930.merge(
        cems, how="left", on=["ba_code", "fuel_category", "datetime_utc"]
    )
    # only keep rows where local datetime is in the current year
    combined_data = combined_data[
        combined_data["datetime_local"].apply(lambda x: x.year) == year
    ]

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
        scaling_factors.groupby(["ba_code", "fuel_category"])["scaling_factor"]
        .min()
        .reset_index()
    )

    # merge the scaling factor into the combined data
    # for any BA-fuels without a scaling factor, fill with 1 (scale to 100% of the origina data)
    combined_data = combined_data.merge(
        scaling_factors, how="left", on=["ba_code", "fuel_category"]
    ).fillna(1)

    # calculate the scaled cems data
    combined_data["cems_scaled"] = (
        combined_data["net_generation_mwh"] * combined_data["scaling_factor"]
    )

    # calculate a scaled residual
    combined_data["residual_scaled"] = (
        combined_data["net_generation_mwh_930"] - combined_data["cems_scaled"]
    )

    return combined_data[
        [
            "ba_code",
            "fuel_category",
            "datetime_utc",
            "datetime_local",
            "report_date",
            "residual_scaled",
        ]
    ]


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

