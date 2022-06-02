import src.data_cleaning as data_cleaning
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

