import pandas as pd
import numpy as np

import load_data
import emissions
import validation
from filepaths import outputs_folder, manual_folder, results_folder
from column_checks import get_dtypes

oge_version_number = "0.2.1"

zero_carbon_fuels = [
    "solar",
    "wind",
    "hydro",
    "nuclear",
    "variable_renewables",
    "batteries",
    "power_storage",
    "storage",
]


def create_egrid_factors_for_eia_regions(egrid_years):
    """This function runs the data processing pipeline to create EFs for EIA regions from eGRID"""
    # load energy source group mapping and merge with egrid
    eia_fuel_categories = pd.read_csv(
        manual_folder("energy_source_groups.csv"), dtype=get_dtypes()
    )[["energy_source_code", "fuel_category_eia930"]].rename(
        columns={"fuel_category_eia930": "fuel_category"}
    )

    # load the egrid data and calculate all relevant emissions
    egrid_plant_all = []
    for year in egrid_years:
        # load the plant data
        egrid_plant = load_egrid_plant_data_for_year(
            year, ipcc_version="AR4", gwp_horizon=100
        )

        # Assign each plant to to a fuel category
        egrid_plant = egrid_plant.merge(
            eia_fuel_categories,
            how="left",
            left_on="plant_primary_fuel",
            right_on="energy_source_code",
            validate="m:1",
        )

        egrid_plant = remove_anomalous_plant_values(egrid_plant)

        # add national totals
        egrid_plant = add_national_total_values(egrid_plant)

        egrid_plant_all.append(egrid_plant)

    # concat all of the years together
    egrid_plant_all = pd.concat(egrid_plant_all, axis=0)

    # aggregate the data to regions and calculate emission factors
    egrid_factors = calculate_fleet_emission_factors(egrid_plant_all)

    egrid_factors = format_df_for_output(
        egrid_factors,
        fuel_mix_source="EIA",
        id_columns=[
            "year",
            "fuel_mix_source",
            "region",
            "fuel_category",
            "emission_factor_source",
        ],
    )

    # validate factors
    validate_factors(egrid_factors)

    return egrid_factors


def create_egrid_factors_for_iso_regions(iso_list, egrid_years, fuel_mix_source):
    """This function runs the data processing pipeline to create EFs for ISO regions from eGRID"""
    iso_egrid_factors = []
    for ba in iso_list:
        print(ba)
        # load and concat multiple years of data
        ba_year_data = []
        for year in egrid_years:

            # load the egrid data and only keep data for plants in the ba
            egrid_plant = load_egrid_plant_data_for_year(
                year, ipcc_version="AR4", gwp_horizon=100
            )
            egrid_plant = egrid_plant[egrid_plant["ba_code"] == ba]

            # load iso-specific fuel categories and map ids
            fuel_category_table = load_plant_fuel_category_for_iso(
                ba, fuel_mix_source, year
            )

            # merge fuel categories into egrid data
            egrid_plant = egrid_plant.merge(
                fuel_category_table,
                how="left",
                on="plant_id_eia",
                validate="m:1",
            )

            egrid_plant = remove_anomalous_plant_values(egrid_plant)

            ba_year_data.append(egrid_plant)

        # concat all of the years together
        ba_year_data = pd.concat(ba_year_data, axis=0)

        # append the data for each region-year to the larger dataframe
        iso_egrid_factors.append(ba_year_data)

    # concat all of the data together
    iso_egrid_factors = pd.concat(iso_egrid_factors, axis=0)

    # aggregate the data to regions and calculate emission factors
    iso_egrid_factors = calculate_fleet_emission_factors(iso_egrid_factors)

    iso_egrid_factors = format_df_for_output(
        iso_egrid_factors,
        fuel_mix_source=fuel_mix_source,
        id_columns=[
            "year",
            "fuel_mix_source",
            "region",
            "fuel_category",
            "emission_factor_source",
        ],
    )

    # validate factors
    validate_factors(iso_egrid_factors)

    return iso_egrid_factors


def create_oge_factors_for_eia_regions(years_to_load):
    """This function runs the data processing pipeline to create EFs for EIA regions from OGE"""
    # get a list of all balancing areas for which there is data in EIA-930
    ba_reference = load_data.load_ba_reference()
    bas_in_eia930 = ba_reference[
        (~ba_reference["timezone_reporting_eia930"].isna())
        & (ba_reference["us_ba"] == "Yes")
        & (~(ba_reference["ba_category"] == "miscellaneous"))
        & (~(ba_reference["retirement_date"].dt.year < min(years_to_load)))
    ]
    ba_list = list(bas_in_eia930.ba_code.unique())

    plant_data = []
    for year in years_to_load:
        # load plant attributes
        plant_attributes = pd.read_csv(
            outputs_folder(f"{year}/plant_static_attributes_{year}.csv"),
            dtype=get_dtypes(),
        )[["plant_id_eia", "fuel_category_eia930", "ba_code"]].rename(
            columns={"fuel_category_eia930": "fuel_category"}
        )

        # load the annual average plant data
        plant_data_year = pd.read_csv(
            results_folder(f"{year}/plant_data/annual/us_units/plant_data.csv")
        )
        # add a year column to the data
        plant_data_year["year"] = year
        plant_data_year["emission_factor_source"] = f"OGEv{oge_version_number}"

        # merge the plant atributes into the plant data
        plant_data_year = plant_data_year.merge(
            plant_attributes, how="left", on="plant_id_eia", validate="1:1"
        )

        # only keep data for the relevant BAs
        plant_data_year = plant_data_year[plant_data_year["ba_code"].isin(ba_list)]

        plant_data_year = remove_anomalous_plant_values(plant_data_year)

        # add national totals
        plant_data_year = add_national_total_values(plant_data_year)

        plant_data.append(plant_data_year)

    plant_data = pd.concat(plant_data, axis=0)

    # aggregate the data to regions and calculate emission factors
    oge_factors_eia = calculate_fleet_emission_factors(plant_data)

    oge_factors_eia = format_df_for_output(
        oge_factors_eia,
        fuel_mix_source="EIA",
        id_columns=[
            "year",
            "fuel_mix_source",
            "region",
            "fuel_category",
            "emission_factor_source",
        ],
    )

    # validate factors
    validate_factors(oge_factors_eia)

    return oge_factors_eia


def create_oge_factors_for_iso_regions(iso_list, years_to_load, fuel_mix_source):
    """This function runs the data processing pipeline to create EFs for ISO regions from OGE"""
    oge_factors_iso = []
    for ba in iso_list:
        print(ba)
        # load and concat multiple years of data
        ba_year_data = []
        for year in years_to_load:

            # load the annual average plant data
            plant_data_year = pd.read_csv(
                results_folder(f"{year}/plant_data/annual/us_units/plant_data.csv")
            )
            # add a year column to the data
            plant_data_year["year"] = year
            plant_data_year["emission_factor_source"] = f"OGEv{oge_version_number}"

            # map each plant to a ba and filter
            plant_ba = pd.read_csv(
                outputs_folder(f"{year}/plant_static_attributes_{year}.csv"),
                dtype=get_dtypes(),
            )[["plant_id_eia", "ba_code"]]
            # merge the plant atributes into the plant data
            plant_data_year = plant_data_year.merge(
                plant_ba, how="left", on="plant_id_eia", validate="1:1"
            )
            plant_data_year = plant_data_year[plant_data_year["ba_code"] == ba]

            # load iso-specific fuel categories and map ids
            fuel_category_table = load_plant_fuel_category_for_iso(
                ba, fuel_mix_source, year
            )

            # merge fuel categories into plant data
            plant_data_year = plant_data_year.merge(
                fuel_category_table,
                how="left",
                on="plant_id_eia",
                validate="m:1",
            )

            plant_data_year = remove_anomalous_plant_values(plant_data_year)

            ba_year_data.append(plant_data_year)

        # concat all of the years together
        ba_year_data = pd.concat(ba_year_data, axis=0)

        # append the data for each region-year to the larger dataframe
        oge_factors_iso.append(ba_year_data)

    # concat all of the data together
    oge_factors_iso = pd.concat(oge_factors_iso, axis=0)

    # aggregate the data to regions and calculate emission factors
    oge_factors_iso = calculate_fleet_emission_factors(oge_factors_iso)

    oge_factors_iso = format_df_for_output(
        oge_factors_iso,
        fuel_mix_source=fuel_mix_source,
        id_columns=[
            "year",
            "fuel_mix_source",
            "region",
            "fuel_category",
            "emission_factor_source",
        ],
    )

    # validate factors
    validate_factors(oge_factors_iso)

    return oge_factors_iso


def remove_anomalous_plant_values(df):
    """Checks for and removes plants with a co2 factor > 20,000 lb/MWh"""
    df["co2_factor"] = df["co2_mass_lb_for_electricity"] / df["net_generation_mwh"]

    anomalies = df[df["co2_factor"] >= 20000]
    if len(anomalies) > 0:
        print("Anomalous plant data detected. Removing from data.")
        print(
            anomalies[
                [
                    "plant_id_eia",
                    "ba_code",
                    "fuel_category",
                    "year",
                    "co2_factor",
                    "co2_mass_lb_for_electricity",
                    "net_generation_mwh",
                ]
            ].sort_values(
                by=[
                    "year",
                    "ba_code",
                    "fuel_category",
                    "plant_id_eia",
                ]
            )
        )
        df = df[df["co2_factor"] < 20000]

    df = df.drop(columns="co2_factor")

    return df


def load_egrid_plant_data_for_year(year, ipcc_version, gwp_horizon):
    egrid_plant = validation.load_egrid_plant_file(year)

    # add a column for eia plant id
    # create a map of eia plant ids to egrid plant ids
    egrid_crosswalk = pd.read_csv(
        manual_folder("eGRID2020_crosswalk_of_EIA_ID_to_EPA_ID.csv")
    )
    egrid_to_eia_id = dict(
        zip(
            list(egrid_crosswalk["plant_id_egrid"]),
            list(egrid_crosswalk["plant_id_eia"]),
        )
    )
    egrid_plant["plant_id_eia"] = egrid_plant["plant_id_egrid"]
    egrid_plant["plant_id_eia"].update(egrid_plant["plant_id_eia"].map(egrid_to_eia_id))

    # if there is a missing value for electric allocation factor, fill with 100%
    egrid_plant["chp_electric_allocation_factor"] = egrid_plant[
        "chp_electric_allocation_factor"
    ].fillna(1)

    # calculate _for_electricity values
    for pol in ["co2", "ch4", "n2o", "nox", "so2"]:
        egrid_plant[f"{pol}_mass_lb_for_electricity"] = (
            egrid_plant[f"{pol}_mass_lb"]
            * egrid_plant["chp_electric_allocation_factor"]
        )

    # egrid 2018-2020 uses the AR4 GWP
    # code adapted from `emissions.calculate_co2e_mass()`
    df_gwp = load_data.load_ipcc_gwp()
    gwp_to_use = df_gwp[df_gwp.ipcc_version == ipcc_version]
    ch4_gwp = gwp_to_use.loc[
        (gwp_to_use.gwp_horizon == gwp_horizon) & (gwp_to_use.gas == "ch4"),
        "gwp",
    ].item()
    n2o_gwp = gwp_to_use.loc[
        (gwp_to_use.gwp_horizon == gwp_horizon) & (gwp_to_use.gas == "n2o"),
        "gwp",
    ].item()

    egrid_plant["co2e_mass_lb_for_electricity"] = (
        egrid_plant["co2_mass_lb_for_electricity"]
        + (ch4_gwp * egrid_plant["ch4_mass_lb_for_electricity"].fillna(0))
        + (n2o_gwp * egrid_plant["n2o_mass_lb_for_electricity"].fillna(0))
    )

    egrid_plant["year"] = year
    if year == 2018:
        egrid_plant["emission_factor_source"] = f"eGRID{year}v2"
    else:
        egrid_plant["emission_factor_source"] = f"eGRID{year}"

    return egrid_plant


def add_national_total_values(plant_data):
    # calculate national average values for each year
    national_data = (
        plant_data.dropna()
        .groupby(["year", "emission_factor_source", "fuel_category"], dropna=False)
        .sum()
        .reset_index()
    )
    national_data["ba_code"] = "national"

    # add the national factors
    plant_data = pd.concat([plant_data, national_data], axis=0)

    return plant_data


def calculate_fleet_emission_factors(egrid_plant):

    data_columns = [
        "net_generation_mwh",
        "co2_mass_lb_for_electricity",
        "ch4_mass_lb_for_electricity",
        "n2o_mass_lb_for_electricity",
        "co2e_mass_lb_for_electricity",
        "nox_mass_lb_for_electricity",
        "so2_mass_lb_for_electricity",
        "co2_mass_lb_for_electricity_adjusted",
        "co2e_mass_lb_for_electricity_adjusted",
        "nox_mass_lb_for_electricity_adjusted",
        "so2_mass_lb_for_electricity_adjusted",
    ]

    # groupby and do same calculations as below
    # explicitly drop na ba and fuel values
    egrid_plant = (
        egrid_plant.groupby(
            ["year", "emission_factor_source", "ba_code", "fuel_category"], dropna=True
        )[data_columns]
        .sum()
        .reset_index()
    )

    # if there are any negative generation values, replace with zero
    egrid_plant.loc[egrid_plant["net_generation_mwh"] < 0, "net_generation_mwh"] = 0

    # calculate emission factors
    for pol in ["co2", "co2e", "nox", "so2"]:
        for pol_type in ["for_electricity", "for_electricity_adjusted"]:
            egrid_plant[f"generated_{pol}_rate_lb_per_mwh_{pol_type}"] = (
                egrid_plant[f"{pol}_mass_lb_{pol_type}"]
                / egrid_plant["net_generation_mwh"]
            )

    factor_columns = [
        "generated_co2_rate_lb_per_mwh_for_electricity",
        "generated_co2e_rate_lb_per_mwh_for_electricity",
        "generated_nox_rate_lb_per_mwh_for_electricity",
        "generated_so2_rate_lb_per_mwh_for_electricity",
        "generated_co2_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_co2e_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_nox_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_so2_rate_lb_per_mwh_for_electricity_adjusted",
    ]

    # set all emissions for zero carbon fuels to zero
    egrid_plant.loc[
        egrid_plant["fuel_category"].isin(zero_carbon_fuels), factor_columns
    ] = 0

    # replace inf values with na
    egrid_plant = egrid_plant.replace(np.inf, np.NaN)

    return egrid_plant


def calculate_emissions_for_ba_year_fuel(
    region,
    ba_year_data,
    fuel_category_table,
    year,
    time_groupby_columns=["datetime_utc"],
):
    # assign a fuel category to each plant
    ba_year_data = ba_year_data.merge(
        fuel_category_table[["plant_id_eia", "fuel_category"]],
        how="left",
        on="plant_id_eia",
        validate="m:1",
    )

    # check that there are no missing fuel categories
    if len(ba_year_data[ba_year_data["fuel_category"].isna()]) > 0:
        print(f"Warning: there are missing fuel categories in {region}")

    # replace negative generation values with zero
    ba_year_data.loc[ba_year_data["net_generation_mwh"] < 0, "net_generation_mwh"] = 0

    # aggregate by fuel category
    ba_year_data = (
        ba_year_data.groupby(["fuel_category"] + time_groupby_columns, dropna=False)
        .sum()
        .reset_index()
    )

    # create columns for adjusted emissions for all pollutants
    # biomass adjustment does not affect these pollutants, only co2
    for pol in ["ch4", "n2o", "nox", "so2"]:
        ba_year_data[f"{pol}_mass_lb_for_electricity_adjusted"] = ba_year_data[
            f"{pol}_mass_lb_for_electricity"
        ]

    # calculate co2eq emissions
    ba_year_data = emissions.calculate_co2e_mass(
        ba_year_data,
        year,
        gwp_horizon=100,
        ar5_climate_carbon_feedback=True,
    )

    return ba_year_data


def calculate_month_hour_emission_factors(
    region, ba, ba_data, local_tz_to_use, datetime_to_use, ba_reference
):
    # groupby to make sure we don't have duplicate timestamps
    ba_data = (
        ba_data.groupby(["fuel_category", "datetime_utc"], dropna=False)
        .sum()
        .reset_index()
    )

    # convert the datetime column to a datetime dtype
    ba_data["datetime_utc"] = pd.to_datetime(ba_data["datetime_utc"])

    # load the local timezone to which each ba reports data to EIA-930
    ba_local_tz = ba_reference.loc[
        ba_reference["ba_code"] == ba, local_tz_to_use
    ].values[0]

    # create a local datetime column
    ba_data["datetime_local"] = ba_data["datetime_utc"].dt.tz_convert(ba_local_tz)

    # create columns for month and hour
    ba_data["year"] = ba_data[datetime_to_use].dt.year
    ba_data["month"] = ba_data[datetime_to_use].dt.month
    ba_data["hour"] = ba_data[datetime_to_use].dt.hour

    # groupby month-hour
    ba_data = (
        ba_data.groupby(["fuel_category", "year", "month", "hour"]).sum().reset_index()
    )

    # calculate emission factors
    for pol in ["co2", "co2e", "nox", "so2"]:
        for pol_type in ["for_electricity", "for_electricity_adjusted"]:
            ba_data[f"generated_{pol}_rate_lb_per_mwh_{pol_type}"] = (
                ba_data[f"{pol}_mass_lb_{pol_type}"] / ba_data["net_generation_mwh"]
            )

    # add a column for ba code
    ba_data["ba_code"] = region

    # only keep relevant columns
    key_columns = [
        "ba_code",
        "fuel_category",
        "year",
        "month",
        "hour",
    ]

    factor_columns = [
        "generated_co2_rate_lb_per_mwh_for_electricity",
        "generated_co2e_rate_lb_per_mwh_for_electricity",
        "generated_nox_rate_lb_per_mwh_for_electricity",
        "generated_so2_rate_lb_per_mwh_for_electricity",
        "generated_co2_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_co2e_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_nox_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_so2_rate_lb_per_mwh_for_electricity_adjusted",
    ]

    ba_data = ba_data[key_columns + factor_columns]

    # set all emissions for zero carbon fuels to zero
    ba_data.loc[ba_data["fuel_category"].isin(zero_carbon_fuels), factor_columns] = 0

    return ba_data


def calculate_emission_factors(region, ba_data):

    # calculate emission factors
    for pol in ["co2", "co2e", "nox", "so2"]:
        for pol_type in ["for_electricity", "for_electricity_adjusted"]:
            ba_data[f"generated_{pol}_rate_lb_per_mwh_{pol_type}"] = (
                ba_data[f"{pol}_mass_lb_{pol_type}"] / ba_data["net_generation_mwh"]
            )

    # add a column for ba code
    ba_data["ba_code"] = region

    # only keep relevant columns
    key_columns = [
        "ba_code",
        "fuel_category",
        "year",
    ]

    factor_columns = [
        "generated_co2_rate_lb_per_mwh_for_electricity",
        "generated_co2e_rate_lb_per_mwh_for_electricity",
        "generated_nox_rate_lb_per_mwh_for_electricity",
        "generated_so2_rate_lb_per_mwh_for_electricity",
        "generated_co2_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_co2e_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_nox_rate_lb_per_mwh_for_electricity_adjusted",
        "generated_so2_rate_lb_per_mwh_for_electricity_adjusted",
    ]

    ba_data = ba_data[key_columns + factor_columns]

    # set all emissions for zero carbon fuels to zero
    ba_data.loc[ba_data["fuel_category"].isin(zero_carbon_fuels), factor_columns] = 0

    return ba_data


def load_plant_fuel_category_for_iso(ba, fuel_mix_source, year):
    """Loads the ISO-specific fuel category for each plant in a single BA."""

    # load plant primary fuel and region
    if year < 2019:  # since there is not OGE data prior to 2019, use 2019
        year = 2019
    plant_attributes = pd.read_csv(
        outputs_folder(f"{year}/plant_static_attributes_{year}.csv"),
        dtype=get_dtypes(),
    )[["plant_id_eia", "plant_primary_fuel"]]

    # merge special flags into plant attributes
    special_categories = load_special_category_flags(year)
    plant_attributes = plant_attributes.merge(
        special_categories, how="left", on="plant_id_eia", validate="1:1"
    )

    # merge fuel category into plant attributes
    energy_source_groups = pd.read_csv(
        manual_folder("energy_source_groups.csv"), dtype=get_dtypes()
    )[["energy_source_code", f"fuel_category_{ba}_{fuel_mix_source}"]].rename(
        columns={f"fuel_category_{ba}_{fuel_mix_source}": "fuel_category"}
    )

    plant_attributes = plant_attributes.merge(
        energy_source_groups,
        how="left",
        left_on="plant_primary_fuel",
        right_on="energy_source_code",
        validate="m:1",
    ).drop(columns="energy_source_code")

    # update special fuel categories
    if ba == "PJM":
        # identify multi-fuel plants in PJM
        plant_attributes.loc[
            (plant_attributes["multiple_fuels"] == 1), "fuel_category"
        ] = "multi_fuel"
    if ba == "NYIS":
        # identify multi-fuel plants
        plant_attributes.loc[
            (plant_attributes["multiple_fuels"] == 1), "fuel_category"
        ] = "dual_fuel"
    if (ba == "ERCO") & (fuel_mix_source == "ISO-H"):
        # identify combined cycle gas plants
        plant_attributes.loc[
            (plant_attributes["combined_cycle"] == 1)
            & (plant_attributes["fuel_category"] == "natural_gas"),
            "fuel_category",
        ] = "natural_gas_cc"

    return plant_attributes[["plant_id_eia", "fuel_category"]]


def load_special_category_flags(year):
    """Loads a dataframe identifying whether any plants are multi-fuel or combined-cycle plants"""
    # identify multi fuel and combined cycle gas plants
    gen_types = load_data.load_pudl_table("generators_eia860", year=year)[
        ["plant_id_eia", "generator_id", "multiple_fuels", "prime_mover_code"]
    ]

    # create a column identifying where a generator is part of a combined cycle plant
    cc_codes = ["CA", "CC", "CS", "CT"]
    gen_types["combined_cycle"] = 0
    gen_types.loc[gen_types["prime_mover_code"].isin(cc_codes), "combined_cycle"] = 1

    # if there is a missing multiple fuel flag, assume that the generator does not burn multiple fuels
    gen_types["multiple_fuels"] = gen_types["multiple_fuels"].fillna(0).astype(int)

    # sum the boolean flags
    gen_types = (
        gen_types.groupby(["plant_id_eia"], dropna=False)[
            ["multiple_fuels", "combined_cycle"]
        ]
        .sum()
        .reset_index()
    )

    # replace the summed values with 1 if greater than zero
    gen_types.loc[gen_types["multiple_fuels"] > 0, "multiple_fuels"] = 1
    gen_types.loc[gen_types["combined_cycle"] > 0, "combined_cycle"] = 1

    return gen_types


def fill_missing_egrid_region_fuels(fuel_categories, egrid_factors, year):
    # make sure that all fuel-regions are represented
    # create a dataframe with all expected region-fuels
    expected_fuel_categories = list(fuel_categories.fuel_category.unique())

    if year == 2018:
        ef_source = f"eGRID{year}v2"
    else:
        ef_source = f"eGRID{year}"

    # identify a list of all BAs that exist in that year
    all_egrid_bas = list(egrid_factors["ba_code"].unique())

    complete_categories = []
    for region in all_egrid_bas:
        for fuel in expected_fuel_categories:
            row_dict = pd.DataFrame(
                {
                    "emission_factor_source": [ef_source],
                    "ba_code": [region],
                    "fuel_category": [fuel],
                }
            )
            complete_categories.append(row_dict)

    complete_categories = pd.concat(complete_categories, axis=0)

    egrid_factors = egrid_factors.merge(
        complete_categories,
        how="outer",
        on=["emission_factor_source", "ba_code", "fuel_category"],
    )

    return egrid_factors


def add_national_average_oge_values(oge_factors):
    # calculate national average values for each year
    national_factors = (
        oge_factors.dropna()
        .groupby(
            ["emission_factor_source", "fuel_category", "year", "month", "hour"],
            dropna=False,
        )
        .mean()
        .reset_index()
    )
    national_factors["ba_code"] = "national"

    # add the national factors
    oge_factors = pd.concat([oge_factors, national_factors], axis=0)

    return oge_factors


def format_df_for_output(df, fuel_mix_source, id_columns):
    df["fuel_mix_source"] = fuel_mix_source
    # change "ba_code" to region
    df = df.rename(columns={"ba_code": "region"})

    # only keep relevant columns
    df = df[id_columns + list(df.filter(like="generated_").columns)]

    # change to long format
    df = df.melt(
        id_vars=id_columns,
        var_name="column_name",
        value_name="emission_factor_lb_per_mwh",
    )

    df["pollutant"] = df["column_name"].str.split("_", expand=True)[1]
    df = df.assign(
        emission_factor_adjustment=lambda x: np.where(
            x.column_name.str.contains("_adjusted"),
            "for_electricity_adjusted",
            "for_electricity",
        )
    )

    # update the ERCO-H iso to ISO-H
    df.loc[df["region"] == "ERCO-H", "fuel_mix_source"] = "ISO-H"
    df.loc[df["region"] == "ERCO-H", "region"] = "ERCO"

    # re order columns
    df = df[
        id_columns
        + ["pollutant", "emission_factor_adjustment", "emission_factor_lb_per_mwh"]
    ]

    return df


def validate_factors(formatted_df):
    # check for negative values
    negative_check = formatted_df[(formatted_df["emission_factor_lb_per_mwh"] < 0)]
    if len(negative_check) > 0:
        print("WARNING: Negative EFs detected")
        print(negative_check)

    # check for missing national backstop values
    backstop_check = formatted_df[
        (formatted_df["region"] == "national")
        & (formatted_df["emission_factor_lb_per_mwh"].isna())
    ]
    if len(backstop_check) > 0:
        print("WARNING: Missing national backstop values detected")
        print(backstop_check)


def screen_for_anomalous_factors(formatted_df, filter_factor):
    """
    This function screens for anomalous values by identifying values that are more than a factor of n larger or smaller than the median.
    """

    anomalous_screen = formatted_df.copy()
    # add a median value
    anomalous_screen["median_factor"] = anomalous_screen.groupby(
        [
            "fuel_mix_source",
            "region",
            "fuel_category",
            "pollutant",
            "emission_factor_adjustment",
        ]
    )["emission_factor_lb_per_mwh"].transform(np.median)
    # add a min value
    anomalous_screen["min_factor"] = anomalous_screen.groupby(
        [
            "fuel_mix_source",
            "region",
            "fuel_category",
            "pollutant",
            "emission_factor_adjustment",
        ]
    )["emission_factor_lb_per_mwh"].transform(np.min)
    # add a max value
    anomalous_screen["max_factor"] = anomalous_screen.groupby(
        [
            "fuel_mix_source",
            "region",
            "fuel_category",
            "pollutant",
            "emission_factor_adjustment",
        ]
    )["emission_factor_lb_per_mwh"].transform(np.max)

    anomalous_screen = anomalous_screen[
        (
            (
                anomalous_screen["emission_factor_lb_per_mwh"]
                < (anomalous_screen["median_factor"] / filter_factor)
            )
            | (
                anomalous_screen["emission_factor_lb_per_mwh"]
                > (anomalous_screen["median_factor"] * filter_factor)
            )
        )
        & (
            (
                anomalous_screen["emission_factor_lb_per_mwh"]
                + anomalous_screen["median_factor"]
            )  # ignore tiny differences
            > 1
        )
    ].sort_values(
        by=[
            "fuel_mix_source",
            "region",
            "fuel_category",
            "pollutant",
            "emission_factor_adjustment",
            "year",
        ]
    )

    if len(anomalous_screen) > 0:
        print("Potentially anomalous values detected")

    return anomalous_screen
