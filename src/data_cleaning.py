import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from pandas import DataFrame
import sqlalchemy as sa
import warnings

import src.load_data as load_data

import pudl.analysis.allocate_net_gen as allocate_gen_fuel


def crosswalk_epa_eia_plant_ids(cems, year):
    """
    Adds a column to the CEMS data that matches the EPA plant ID to the EIA plant ID
    Inputs:
        cems: pandas dataframe with hourly emissions data and columns for "plant_id_epa" and "unitid"
    Returns:
        cems: pandas dataframe with an additional column for "plant_id_eia"
    """

    psdc = load_data.load_epa_eia_crosswalk(year)

    # create a table that matches EPA plant and unit IDs to an EIA plant ID
    plant_id_crosswalk = psdc[
        ["plant_id_epa", "unitid", "plant_id_eia", "generator_id"]
    ].drop_duplicates()

    # only keep plant ids where the two are different
    plant_id_crosswalk = plant_id_crosswalk[
        plant_id_crosswalk["plant_id_epa"] != plant_id_crosswalk["plant_id_eia"]
    ].dropna()

    # match plant_id_eia on plant_id_epa and unitid
    cems = cems.merge(plant_id_crosswalk, how="left", on=["plant_id_epa", "unitid"])

    # if the merge resulted in any missing plant_id associations, fill with the plant_id_epa, assuming that they are the same
    cems["plant_id_eia"] = cems["plant_id_eia"].fillna(cems["plant_id_epa"])

    # change the id column from float dtype to int
    cems["plant_id_eia"] = cems["plant_id_eia"].astype(int)

    return cems


def remove_plants(
    df,
    non_grid_connected=False,
    remove_states=[],
    steam_only_plants=False,
    distribution_connected_plants=False,
):
    """
    Coordinating function to remove specific plants based on specified options
    Each function should identify how many plants are being removed
    Args:
        df: dataframe containing plant_id_eia column
        non_grid_connected: if True, remove all plants that are not grid connected
        remove_states: list of two-letter state codes for which plants should be removed if located within
        steam_only_plants: if True, remove plants that only generate heat and no electricity (not yet implemented)
        distribution_connected_plants: if True, remove plants that are connected to the distribution grid (not yet implemented)
    """
    if non_grid_connected:
        df = remove_non_grid_connected_plants(df)
    if len(remove_states) > 0:
        plant_states = (
            load_data.initialize_pudl_out()
            .plants_eia860()
            .loc[:, ["plant_id_eia", "state"]]
        )
        plants_in_states_to_remove = list(
            plant_states[
                plant_states["state"].isin(remove_states)
            ].plant_id_eia.unique()
        )
        print(
            f"Removing {len(plants_in_states_to_remove)} plants located in the following states: {remove_states}"
        )
        df = df[~df["plant_id_eia"].isin(plants_in_states_to_remove)]
    if steam_only_plants:
        pass
    if distribution_connected_plants:
        pass

    return df


def remove_non_grid_connected_plants(df):
    """
    Removes any records from a dataframe associated with plants that are not connected to the electricity grid
    Inputs: 
        df: any pandas dataframe containing the column 'plant_id_eia'
    Returns:
        df: pandas dataframe with non-grid connected plants removed
    """

    # get the list of plant_id_eia from the static table
    ngc_plants = list(
        pd.read_csv(
            "../data/egrid/egrid_static_tables/table_4-2_plants_not_connected_to_grid.csv"
        )["Plant ID"]
    )

    num_plants = len(
        df[df["plant_id_eia"].isin(ngc_plants)]["plant_id_eia"].unique()
    ) + len(
        df[(df["plant_id_eia"] >= 880000) & (df["plant_id_eia"] < 890000)][
            "plant_id_eia"
        ].unique()
    )
    print(f"Removing {num_plants} plants that are not grid-connected")

    df = df[~df["plant_id_eia"].isin(ngc_plants)]

    # according to the egrid documentation, any plants that have an id of 88XXXX are not grid connected
    # only keep plants that dont have an id of 88XXXX
    df = df[(df["plant_id_eia"] < 880000) | (df["plant_id_eia"] >= 890000)]

    return df


def manually_remove_steam_units(df):
    """
    Removes any records from CEMS that we've identified as being steam only plants that need to be removed
    """

    # get the list of plant_id_eia from the static table
    units_to_remove = list(
        pd.read_csv("../data/egrid/egrid_static_tables/steam_units_to_remove.csv")[
            "cems_id"
        ]
    )

    print(
        f"Removing {len(units_to_remove)} units that only produce steam and do not report to EIA"
    )

    df = df[~df["cems_id"].isin(units_to_remove)]

    return df


def remove_heating_only_plants(cems):
    """
    Removes plants from the cems data that only report steam generation and no electrical generation
    Inputs:
        cems: pandas dataframe containing hourly CEMS data
    Returns:
        cems: pandas dataframe with steam-only plants removed

    """

    # create a list of plants that report only steam generation but no electrical generation
    cems_annual = cems.groupby(["plant_id_eia"]).sum()
    steam_only_cems_plant_ids = list(
        cems_annual[
            (cems_annual["gross_load_mw"] == 0)
            & (cems_annual["steam_load_1000_lbs"] > 0)
        ].index
    )

    # remove these plants from the cems data
    num_plants = len(
        cems[cems["plant_id_eia"].isin(steam_only_cems_plant_ids)][
            "plant_id_eia"
        ].unique()
    )
    print(f"Removing {num_plants} plants that only produce heat and no power")
    cems = cems[~cems["plant_id_eia"].isin(steam_only_cems_plant_ids)]

    return cems


def determine_cems_reporting_status(cems):
    """
    Determines whether a plant that reports to CEMS reports for the entire year, or only partial year
    Inputs:
        cems: pandas dataframe with hourly cems data
    Returns:
        cems: pandas dataframe with additional column added for cems_reporting_category
    """
    # sum CEMS data by month for each unit
    cems_monthly = (
        cems.groupby(["cems_id", "report_date"])
        .sum()[
            [
                "operating_time_hours",
                "gross_load_mw",
                "steam_load_1000_lbs",
                "co2_mass_tons",
                "fuel_consumed_mmbtu",
            ]
        ]
        .reset_index()
    )

    # identify all of the plants that report to CEMS in all 12 months
    full_year_reporters = (
        cems_monthly.groupby(["cems_id"])
        .count()
        .query("report_date == 12")
        .reset_index()
    )
    full_year_reporters["cems_reporting_category"] = "full_year"

    # add this data to the cems data
    cems = cems.merge(
        full_year_reporters[["cems_id", "cems_reporting_category"]],
        how="left",
        on=["cems_id"],
    )

    cems["cems_reporting_category"] = cems["cems_reporting_category"].fillna(
        "partial_year"
    )

    return cems


def get_epa_unit_fuel_types(year):
    """
    Loads the energy source code assigned to each CAMD unit in the EPA-EIA crosswalk.

    If the EIA fuel type is missing, uses the CAMD fuel type to fill.
    """
    # get a unique list of plant unit fuels
    fuel_types = load_data.load_epa_eia_crosswalk(year)[
        ["plant_id_epa", "unitid", "energy_source_code_eia"]
    ].drop_duplicates()

    fuel_types = fuel_types.dropna(subset="energy_source_code_eia")

    # remove any entries where there are multiple fuel types listed
    fuel_types = fuel_types[
        ~fuel_types[["plant_id_epa", "unitid"]].duplicated(keep=False)
    ]

    # rename the column
    fuel_types = fuel_types.rename(
        columns={"energy_source_code_eia": "energy_source_code"}
    )

    return fuel_types


def update_energy_source_codes(df):
    """
    Manually update fuel source codes
    """
    # refinery with energy source = OTH
    df.loc[
        (df["plant_id_eia"] == 50626) & (df["generator_id"] == "GEN1"),
        "energy_source_code",
    ] = "OG"
    df.loc[
        (df["plant_id_eia"] == 56139) & (df["generator_id"] == "NPCG"),
        "energy_source_code",
    ] = "OG"

    return df


def calculate_geothermal_emission_factors(year):
    """
    Updates the list of geothermal plants provided by EPA using EIA data
    Calculates a weighted average EF for each plant-month based on the fraction 
    of fuel consumed from each type of prime mover (steam, binary, flash)
    """
    pudl_db = "sqlite:///../data/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)

    # load the eia generation fuel data
    generation_fuel_eia923 = pd.read_sql(
        f"SELECT * FROM generation_fuel_eia923 WHERE report_date >= '{year}-01-01' AND report_date <= '{year}-12-01'",
        pudl_engine,
    )

    # create a dataframe of total heat input by prime mover for each geothermal plant
    geo_in_eia = (
        generation_fuel_eia923[generation_fuel_eia923["energy_source_code"] == "GEO"]
        .groupby(["plant_id_eia", "prime_mover_code", "report_date"])
        .sum()["fuel_consumed_mmbtu"]
        .reset_index()
    )
    # remove prime movers for which there was no heat input
    geo_in_eia = geo_in_eia[geo_in_eia["fuel_consumed_mmbtu"] > 0]

    # merge in the EPA's assigned Geotype
    geothermal_geotype = pd.read_csv(
        "../data/egrid/egrid_static_tables/table_geothermal_geotype.csv"
    )
    geo_in_eia = geo_in_eia.merge(
        geothermal_geotype[["plant_id_eia", "geotype_code"]],
        how="left",
        on="plant_id_eia",
    )

    # identify plants with multiple prime mover types
    multi_type_plants = (
        geo_in_eia.groupby(["plant_id_eia", "prime_mover_code"])
        .count()
        .reset_index()
        .groupby("plant_id_eia")
        .count()["prime_mover_code"]
    )
    multi_type_plants = multi_type_plants[multi_type_plants > 1]
    multi_type_plants = list(multi_type_plants.index)

    # update the geotype codes for plants with multiple types
    # for plants identified as flash steam that also have a binary component, update to binary
    geo_in_eia.loc[
        (geo_in_eia["plant_id_eia"].isin(multi_type_plants))
        & (geo_in_eia["geotype_code"] == "F")
        & (geo_in_eia["prime_mover_code"] == "BT"),
        "geotype_code",
    ] = "B"
    # for plants identified as binary that also have a steam component, update to flash (it seems that all other multi-types are F/B combinatioms)
    geo_in_eia.loc[
        (geo_in_eia["plant_id_eia"].isin(multi_type_plants))
        & (geo_in_eia["geotype_code"] == "B")
        & (geo_in_eia["prime_mover_code"] == "ST"),
        "geotype_code",
    ] = "F"

    # if EPA assigned a plant as flash or steam, but EIA identified it as binary, re-assign as binary
    geo_in_eia.loc[
        (geo_in_eia["prime_mover_code"] == "BT")
        & (geo_in_eia["geotype_code"].isin(["F", "S"])),
        "geotype_code",
    ] = "B"

    # if EPA assigned a plant as binary, but EIA identified it as a steam turbine, re-assign as flash
    # we use flash instead of steam, b/c flash is more common than steam according to EIA
    # Source: https://www.eia.gov/energyexplained/geothermal/geothermal-power-plants.php
    geo_in_eia.loc[
        (geo_in_eia["prime_mover_code"] == "ST")
        & (geo_in_eia["geotype_code"].isin(["B"])),
        "geotype_code",
    ] = "F"

    # where plants are missing a geotype code, assign based on the EIA-identified prime mover
    geo_in_eia.loc[
        (geo_in_eia["geotype_code"].isna()) & (geo_in_eia["prime_mover_code"] == "BT"),
        "geotype_code",
    ] = "B"
    geo_in_eia.loc[
        (geo_in_eia["geotype_code"].isna()) & (geo_in_eia["prime_mover_code"] == "ST"),
        "geotype_code",
    ] = "F"

    # calculate the fraction of heat input from each prime mover in each month
    fuel_frac = (
        geo_in_eia.set_index(["plant_id_eia", "report_date", "geotype_code"])[
            ["fuel_consumed_mmbtu"]
        ]
        / geo_in_eia.groupby(["plant_id_eia", "report_date"]).sum()
    ).reset_index()
    fuel_frac = fuel_frac.rename(columns={"fuel_consumed_mmbtu": "fuel_frac"})
    geo_in_eia = geo_in_eia.merge(
        fuel_frac, how="left", on=["plant_id_eia", "report_date", "geotype_code"]
    )

    # calculate a weighted average emission factor for each plant

    # load geothermal efs
    geothermal_efs = pd.read_csv(
        "../data/egrid/egrid_static_tables/table_C6_geothermal_emission_factors.csv"
    )[["geotype_code", "co2_lb_per_mmbtu"]]
    # convert lb to ton
    geothermal_efs["co2_tons_per_mmbtu"] = geothermal_efs["co2_lb_per_mmbtu"] / 2000
    geothermal_efs = geothermal_efs[["geotype_code", "co2_tons_per_mmbtu"]]
    # merge in the emission factor
    geo_in_eia = geo_in_eia.merge(geothermal_efs, how="left", on="geotype_code")
    # multiply the emission factor by the fraction
    geo_in_eia["co2_tons_per_mmbtu"] = (
        geo_in_eia["fuel_frac"] * geo_in_eia["co2_tons_per_mmbtu"]
    )

    # groupby plant and month to get the weighted emission factor
    geo_in_eia = (
        geo_in_eia.groupby(["plant_id_eia", "report_date"])
        .sum()["co2_tons_per_mmbtu"]
        .reset_index()
    )

    # if there are any plants missing from our list, add them back in

    # identify the plants that are in the epa geotype table but not the EIA-derived one
    epa_geo_plants = list(geothermal_geotype.plant_id_eia.unique())
    plants_from_eia = list(geo_in_eia.plant_id_eia.unique())
    missing_plants = list(set(epa_geo_plants) - set(plants_from_eia))

    # create a dataframe with the geotype of all misisng plants
    missing_plants = geothermal_geotype.loc[
        geothermal_geotype["plant_id_eia"].isin(missing_plants),
        ["plant_id_eia", "geotype_code"],
    ]

    # merge in the efs
    missing_plants = missing_plants.merge(geothermal_efs, how="left", on="geotype_code")

    # drop the geotype code
    missing_plants = missing_plants.drop(columns=["geotype_code"])

    # create a record for each month of the year
    missing_plants = create_monthly_gens_records(missing_plants, year)

    # concat the missing plants to the other dataframe
    geo_efs = pd.concat([geo_in_eia, missing_plants], axis=0)

    geo_efs["report_date"] = pd.to_datetime(geo_efs["report_date"])

    return geo_efs


def calculate_co2_from_fuel_consumption(df, year):
    """
    Inputs:
        df: pandas dataframe containing the following columns: ['plant_id_eia', 'report_date,'fuel_consumed_mmbtu','energy_source_code']
    """

    # get emission factors
    emission_factors = load_data.load_emission_factors()[
        ["energy_source_code", "co2_tons_per_mmbtu"]
    ]

    # add emission factor to  df
    df = df.merge(emission_factors, how="left", on="energy_source_code")

    geothermal_efs = calculate_geothermal_emission_factors(year).rename(
        columns={"co2_tons_per_mmbtu": "co2_tons_per_mmbtu_geo"}
    )

    # add geothermal emission factor to df
    df = df.merge(geothermal_efs, how="left", on=["plant_id_eia", "report_date"])

    # update missing efs using the geothermal efs if available
    df["co2_tons_per_mmbtu"] = df["co2_tons_per_mmbtu"].fillna(
        df["co2_tons_per_mmbtu_geo"]
    )

    # create a new column with the  co2 mass in tons
    df["co2_mass_tons"] = df["fuel_consumed_mmbtu"] * df["co2_tons_per_mmbtu"]

    # if there is a column for fuel_consumed for electricity, add an adjusted co2 column
    if "fuel_consumed_for_electricity_mmbtu" in df.columns:
        df["co2_mass_tons_adjusted"] = (
            df["fuel_consumed_for_electricity_mmbtu"] * df["co2_tons_per_mmbtu"]
        )

    # drop intermediate columns
    df = df.drop(columns=["co2_tons_per_mmbtu", "co2_tons_per_mmbtu_geo"])

    return df


def calculate_co2e_from_fuel_consumption(df, year):
    """
    Calculate CO2e emissions from fuel consumption data.

    Inputs:
        df: Should contain the following columns:
            [`plant_id_eia`, `report_date`, `fuel_consumed_units`, `energy_source_code`, `prime_mover_code`]
    
    If the `fuel_consumed_for_electricity_units` column is available, we also
    compute the adjusted emissions.
    """
    emission_factors = load_data.load_emission_factors()[
        ["energy_source_code", "co2_tons_per_mmbtu", "ch4_lbs_per_mmbtu", "n2o_lbs_per_mmbtu"]]

    df = df.merge(emission_factors, how="left", on="energy_source_code")

    geothermal_efs = calculate_geothermal_emission_factors(year).rename(
        columns={"co2_tons_per_mmbtu": "co2_tons_per_mmbtu_geo"})

    # add geothermal emission factor to df
    df = df.merge(geothermal_efs, how="left", on=["plant_id_eia", "report_date"])

    # update missing efs using the geothermal efs if available
    df["co2_tons_per_mmbtu"] = df["co2_tons_per_mmbtu"].fillna(df["co2_tons_per_mmbtu_geo"])

    # eGRID was updated to use AR4 in 2018
    if year < 2018:
        gwp_100_ch4 = 23.0
        gwp_100_n2o = 296.0
    else:
        gwp_100_ch4 = 25.0
        gwp_100_n2o = 298.0

    # Compute CO2-eq mass using the same GWP factors as eGRID.
    df["co2e_mass_tons"] = df["fuel_consumed_mmbtu"] * \
        (df["co2_tons_per_mmbtu"] + \
         gwp_100_ch4 * df["ch4_lbs_per_mmbtu"] / 2000 + \
         gwp_100_n2o * df["n2o_lbs_per_mmbtu"] / 2000)

    # if there is a column for fuel_consumed for electricity, add an adjusted co2 column
    if "fuel_consumed_for_electricity_mmbtu" in df.columns:
        df["co2e_mass_tons_adjusted"] = df["fuel_consumed_for_electricity_mmbtu"] * \
            (df["co2_tons_per_mmbtu"] + \
            gwp_100_ch4 * df["ch4_lbs_per_mmbtu"] / 2000 + \
            gwp_100_n2o * df["n2o_lbs_per_mmbtu"] / 2000)

    # drop intermediate columns
    df = df.drop(columns=["co2_tons_per_mmbtu", "ch4_lbs_per_mmbtu", "n2o_lbs_per_mmbtu", "co2_tons_per_mmbtu_geo"])

    return df


def calculate_nox_from_fuel_consumption(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate NOx emissions from fuel consumption data.

    Inputs:
        df: Should contain the following columns:
            [`plant_id_eia`, `report_date`, `fuel_consumed_units`, `energy_source_code`, `prime_mover_code`]
    
    If the `fuel_consumed_for_electricity_units` column is available, we also
    compute the adjusted emissions.
    """
    emission_factors = load_data.load_emission_factors_nox()[
        ["Prime Mover", "Primary Fuel Type", "Emission Factor"]]

    # Add the emission factors to 'df' as a column. Note the the columns names are different.
    df = df.merge(emission_factors, how="left",
                  left_on=["energy_source_code", "prime_mover_code"], 
                  right_on=["Primary Fuel Type", "Prime Mover"])
    df = df.rename(columns={"Emission Factor": "nox_lbs_per_unit", "B": "c"})

    # Create a new column with the NOx mass in lbs. Note that the NOx emission
    # rate numerator units are reported in 'lbs' already.
    df["nox_mass_lbs"] = df["fuel_consumed_units"] * df["nox_lbs_per_unit"]

    # If there is a column for electricity-related fuel consumption, add an adjusted NOx column.
    if "fuel_consumed_for_electricity_units" in df.columns:
        df["nox_mass_lbs_adjusted"] = (df["fuel_consumed_for_electricity_units"] * df["nox_lbs_per_unit"])

    # Drop intermediate columns.
    df = df.drop(columns=["nox_lbs_per_unit"])

    return df


def calculate_so2_from_fuel_consumption(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate SO2 emissions from fuel consumption data and fuel sulfur content.

    Inputs:
        df: Should contain the following columns:
            [`plant_id_eia`, `report_date`, `fuel_consumed_units`, `energy_source_code`, `prime_mover_code`]
    
    If the `fuel_consumed_for_electricity_units` column is available, we also
    compute the adjusted emissions.
    """
    emission_factors = load_data.load_emission_factors_so2()[
        ["Prime Mover", "Primary Fuel Type", "emission_factor_coeff", "multiply_by_sulfur_content"]]

    boiler_fuel_923 = load_data.load_pudl_table('boiler_fuel_eia923')[
        ["plant_id_eia", "report_date", "energy_source_code", "sulfur_content_pct"]]

    # Compute average sulfur contents for each energy_source_code to use as a default.
    mean_sulfur_content_pct = boiler_fuel_923.groupby(["energy_source_code"]).mean()

    default_sulfur_each_row = mean_sulfur_content_pct.loc[boiler_fuel_923.energy_source_code]["sulfur_content_pct"]

    # Need to have aligned indices for fillna to work.
    default_sulfur_each_row.index = boiler_fuel_923.index
    boiler_fuel_923["sulfur_content_pct"] = boiler_fuel_923["sulfur_content_pct"].fillna(default_sulfur_each_row)

    # Add the emission factors to 'df' as a column. Note the the columns names are different.
    df = df.merge(emission_factors, how="left",
                  left_on=["energy_source_code", "prime_mover_code"], 
                  right_on=["Primary Fuel Type", "Prime Mover"])

    # Add the sulfur content to 'df'.
    df = df.merge(boiler_fuel_923, how="left",
                  on=["plant_id_eia", "report_date", "energy_source_code"])

    df["so2_lbs_per_unit"] = df["multiply_by_sulfur_content"] * df["sulfur_content_pct"] * df["emission_factor_coeff"]

    # Create a new column with the SO2 mass in lbs. Note that the SO2 emission
    # rate numerator units are reported in 'lbs' already.
    df["so2_mass_lbs"] = df["fuel_consumed_units"] * df["so2_lbs_per_unit"]

    # If there is a column for electricity-related fuel consumption, add an adjusted SO2 column.
    if "fuel_consumed_for_electricity_units" in df.columns:
        df["so2_mass_lbs_adjusted"] = (df["fuel_consumed_for_electricity_units"] * df["so2_lbs_per_unit"])

    # Drop intermediate columns.
    df = df.drop(columns=["so2_lbs_per_unit", "sulfur_content_pct", "multiply_by_sulfur_content", "emission_factor_coeff"])

    return df


def assign_fuel_type_to_cems(cems, year):
    "Assigns a fuel type to each observation in CEMS"

    fuel_types = get_epa_unit_fuel_types(year)

    # merge in the reported fuel type
    cems = cems.merge(fuel_types, how="left", on=["plant_id_epa", "unitid"])

    # TODO: fill fuel codes for plants that only have a single fossil type identified in EIA
    cems = fill_missing_fuel_for_single_fuel_plant_months(cems, year)

    return cems


def fill_cems_missing_co2(cems, year):
    """
    Fills missing hourly CO2 data in CEMS based on a two-tiered approach.

    CO2 data is considered missing if reported CO2 is zero and fuel consumption is positive.
    If a unit has a unit-specific fuel type identified by the EPA-EIA crosswalk, calculate co2 using a fuel-specific emission factor.
    If not, fill missing data by calculating a plant-month weighted average emission factor of all fuels burned in that plant-month.
    """
    # replace all "missing" CO2 values with zero
    cems["co2_mass_tons"] = cems["co2_mass_tons"].fillna(0)

    # replace 0 reported CO2 values with missing values, if there was reported heat input
    cems.loc[
        (cems["co2_mass_tons"] == 0) & (cems["fuel_consumed_mmbtu"] > 0),
        "co2_mass_tons",
    ] = np.NaN

    # create a new df with all observations with missing co2 data
    missing_co2 = cems[cems["co2_mass_tons"].isnull()]

    #### First round of filling using fuel types in PSDC

    # for rows that have a successful fuel code match, move to a temporary dataframe to hold the data
    co2_to_fill = missing_co2.copy()[~missing_co2["energy_source_code"].isna()]
    fill_index = co2_to_fill.index

    # remove these from the missing co2 dataframe. We'll need to apply a different method for these remaining plants
    missing_co2 = missing_co2[missing_co2["energy_source_code"].isna()]
    missing_index = missing_co2.index

    # calculate emissions based on fuel type
    co2_to_fill = calculate_co2_from_fuel_consumption(co2_to_fill, year).set_index(
        fill_index
    )

    # fill this data into the original cems data
    cems.update(co2_to_fill[["co2_mass_tons"]])

    #### Second round of data filling using weighted average EF based on EIA-923 heat input data

    # get a list of plant ids in the missing data
    missing_plants = list(missing_co2["plant_id_eia"].unique())

    # load 923 data
    generation_fuel_eia923 = load_data.load_pudl_table(
        "generation_fuel_eia923", year=year
    )

    # get monthly fuel data for each of the missing plants
    missing_gf = generation_fuel_eia923[
        generation_fuel_eia923["plant_id_eia"].isin(missing_plants)
    ]

    # calculate total fuel consumed of each fuel type in each month
    missing_gf = missing_gf.groupby(
        ["plant_id_eia", "report_date", "energy_source_code"]
    ).sum()[["fuel_consumed_for_electricity_mmbtu"]]

    # calculate the percent of heat input from each fuel in each month
    missing_gf = (
        missing_gf
        / missing_gf.reset_index().groupby(["plant_id_eia", "report_date"]).sum()
    )

    missing_gf = missing_gf.fillna(1)

    emission_factors = load_data.load_emission_factors()[
        ["energy_source_code", "co2_tons_per_mmbtu"]
    ]

    # merge in the emission factor
    missing_gf = missing_gf.reset_index().merge(
        emission_factors, how="left", on="energy_source_code"
    )

    # calculate weighted emission factor
    missing_gf["weighted_ef"] = (
        missing_gf["fuel_consumed_for_electricity_mmbtu"]
        * missing_gf["co2_tons_per_mmbtu"]
    )
    missing_gf = (
        missing_gf.groupby(["plant_id_eia", "report_date"])
        .sum()["weighted_ef"]
        .reset_index()
    )

    # convert report date back to datetime
    missing_gf["report_date"] = pd.to_datetime(missing_gf["report_date"])

    # merge the weighted ef into the missing data
    missing_co2 = missing_co2.merge(
        missing_gf, how="left", on=["plant_id_eia", "report_date"]
    ).set_index(missing_index)

    # calculate missing co2 data
    missing_co2["co2_mass_tons"] = (
        missing_co2["fuel_consumed_mmbtu"] * missing_co2["weighted_ef"]
    )

    # update in CEMS table
    cems.update(missing_co2[["co2_mass_tons"]])

    return cems


def fill_missing_fuel_for_single_fuel_plant_months(df, year):
    """
    Identifies all plant-months where a single fuel was burned based on the EIA-923 generation fuel table
    Uses this to fill in the energy source code if a match was not made based on the PSDC
    """

    # identify plant-months for which there is a single fossil fuel type reported
    gf = load_data.load_pudl_table("generation_fuel_eia923", year=year)[
        ["plant_id_eia", "report_date", "energy_source_code", "fuel_consumed_mmbtu"]
    ]

    # remove any rows for clean fuels
    clean_fuels = ["SUN", "MWH", "WND", "WAT", "WH", "PUR", "NUC"]
    gf = gf[~gf["energy_source_code"].isin(clean_fuels)]

    # group the data by plant, month, and fuel
    gf = (
        gf.groupby(["plant_id_eia", "report_date", "energy_source_code"])
        .sum()
        .reset_index()
    )

    # only keep rows with >0 fuel consumed reported
    gf = gf[gf["fuel_consumed_mmbtu"] > 0]

    # identify which plant months have multiple fuels reported
    multi_fuels = (
        gf.groupby(["plant_id_eia", "report_date"])["energy_source_code"]
        .count()
        .reset_index()
        .rename(columns={"energy_source_code": "num_fuels"})
    )

    # merge this information back into the other dataframe
    gf = gf.merge(multi_fuels, how="left", on=["plant_id_eia", "report_date"])

    # only keep rows that have a single fuel type in a month
    gf = gf[gf["num_fuels"] == 1]

    # clean up the columns
    gf = gf.rename(columns={"energy_source_code": "energy_source_code_single"}).drop(
        columns=["num_fuels", "fuel_consumed_mmbtu"]
    )
    gf["report_date"] = pd.to_datetime(gf["report_date"])

    # merge this data into the df
    df = df.merge(gf, how="left", on=["plant_id_eia", "report_date"])

    # fill missing fuel types with this data
    df["energy_source_code"] = df["energy_source_code"].fillna(
        df["energy_source_code_single"]
    )

    # remove the intermediate column
    df = df.drop(columns=["energy_source_code_single"])

    return df


def crosswalk_epa_unit_to_eia_generator_id(df, year, unique_gen_match=False):
    """
    Crosswalks the EPA unitid to the EIA generator_id. NOTE: there may be multiple generators associated with each unit
    Inputs:
        df: pandas dataframe with the columns ['plant_id_eia','unitid']
        unique_gen_match: T/F, whether to only keep where one or more units map to a single generator
    Returns:
        df with new column for 'generator_id' (May have duplicate records for each unitid)
    """

    # load the power sector data crosswalk
    psdc = load_data.load_epa_eia_crosswalk(year)

    # create a table that matches EPA plant and unit IDs to an EIA plant ID
    unit_generator_crosswalk = psdc[
        ["plant_id_eia", "plant_id_epa", "unitid", "generator_id"]
    ].drop_duplicates()

    # fill any missing eia plant ids with epa plant ids
    unit_generator_crosswalk["plant_id_eia"] = unit_generator_crosswalk[
        "plant_id_eia"
    ].fillna(unit_generator_crosswalk["plant_id_epa"])

    # drop the plant_id_epa column
    unit_generator_crosswalk = unit_generator_crosswalk.drop(columns="plant_id_epa")

    if unique_gen_match == True:
        unit_generator_crosswalk = unit_generator_crosswalk.drop_duplicates(
            subset=["plant_id_eia", "unitid"], keep=False
        )

    df = df.merge(unit_generator_crosswalk, how="left", on=["plant_id_eia", "unitid"])

    return df


def remove_cems_with_zero_monthly_data(cems):
    """
    Identifies months where zero generation, heat input, or emissions are reported 
    from each unit and removes associated hours from CEMS so that these can be filled using the eia923 data
    Inputs:
        cems: pandas dataframe of hourly cems data containing columns "cems_id" and "report_date"
    Returns:
        cems df with hourly observations for months when no emissions reported removed
    """
    # calculate teh totals reported in each month
    cems_with_zero_monthly_emissions = cems.groupby(["cems_id", "report_date"]).sum()[
        ["co2_mass_tons", "gross_generation_mwh", "fuel_consumed_mmbtu"]
    ]
    # identify unit-months where zero emissions reported
    cems_with_zero_monthly_emissions = cems_with_zero_monthly_emissions[
        cems_with_zero_monthly_emissions.sum(axis=1) == 0
    ]
    # add a flag to these observations
    cems_with_zero_monthly_emissions["missing_data_flag"] = "remove"

    # merge the missing data flag into the cems data
    cems = cems.merge(
        cems_with_zero_monthly_emissions.reset_index()[
            ["cems_id", "report_date", "missing_data_flag"]
        ],
        how="left",
        on=["cems_id", "report_date"],
    )
    # remove any observations with the missing data flag
    print(
        f"removing {len(cems[cems['missing_data_flag'] == 'remove'])} observations from cems for unit-months where no data reported"
    )
    cems = cems[cems["missing_data_flag"] != "remove"]
    # drop the missing data flag column
    cems = cems.drop(columns="missing_data_flag")

    return cems


def identify_emissions_data_source(cems, gen_fuel_allocated, year):
    """
    For each generator-month record in gen_fuel_allocated, identify whether hourly cems data exists
    The monthly records that don't have cems data are what we will need to assign an hourly profile to
    """

    # Step 1: Match based on the EPA EIA Crosswalk
    ##############################################

    # aggregate cems data to plant-unit-month
    cems_monthly = (
        cems.groupby(["plant_id_eia", "unitid", "report_date"])
        .sum()[["gross_generation_mwh", "co2_mass_tons", "fuel_consumed_mmbtu"]]
        .reset_index()
    )

    # crosswalk this data with each generator id
    cems_monthly = crosswalk_epa_unit_to_eia_generator_id(cems_monthly, year)

    # rename the columns
    # cems_monthly = cems_monthly.rename(columns={'gross_generation_mwh':'cems_gross_generation_mwh','fuel_consumed_mmbtu':'cems_fuel_consumed_mmbtu','co2_mass_tons':'cems_co2_mass_tons'})

    # create a dataframe containing all generator-months with data reported to cems
    generator_months_in_cems = cems_monthly[
        ["plant_id_eia", "generator_id", "report_date"]
    ].drop_duplicates()
    generator_months_in_cems["data_source"] = "cems"

    # identify which generation and fuel data is not reported in cems
    gen_fuel_allocated = gen_fuel_allocated.merge(
        generator_months_in_cems,
        how="left",
        on=["plant_id_eia", "generator_id", "report_date"],
    )

    # identify all generators that report to cems in at least one month
    generator_months_in_cems["plant_gen_id"] = (
        generator_months_in_cems["plant_id_eia"].astype(str)
        + "_"
        + generator_months_in_cems["generator_id"].astype(str)
    )
    gens_in_cems = list(generator_months_in_cems["plant_gen_id"].unique())

    # for all months where a generator that reports to cems does not report data, fill the data type as eia
    # this prevents accidental identification of cems as the data source in the next step
    gen_fuel_allocated["plant_gen_id"] = (
        gen_fuel_allocated["plant_id_eia"].astype(str)
        + "_"
        + gen_fuel_allocated["generator_id"].astype(str)
    )
    gen_fuel_allocated.loc[
        (gen_fuel_allocated["plant_gen_id"].isin(gens_in_cems))
        & (gen_fuel_allocated["data_source"].isna()),
        "data_source",
    ] = "eia_only"

    # drop the plant gen id column
    gen_fuel_allocated = gen_fuel_allocated.drop(columns=["plant_gen_id"])

    # Step 2: Match based on amount of heat input reported in each source
    ####################################################################

    # aggregate data by plant, month, and fuel type
    eia_plant_fuel = (
        gen_fuel_allocated.groupby(
            ["plant_id_eia", "energy_source_code", "report_date"], dropna=False
        )
        .sum()
        .reset_index()
    )
    cems_plant_fuel = (
        cems.groupby(
            ["plant_id_eia", "energy_source_code", "report_date"], dropna=False
        )
        .sum()[["gross_generation_mwh", "fuel_consumed_mmbtu", "co2_mass_tons"]]
        .reset_index()
    )

    # merge the data together so that we can compare heat input reported by the two sources
    columns_to_match = [
        "plant_id_eia",
        "energy_source_code",
        "report_date",
        "fuel_consumed_mmbtu",
    ]
    fuel_comparison = eia_plant_fuel[columns_to_match].merge(
        cems_plant_fuel[columns_to_match],
        how="left",
        on=["plant_id_eia", "energy_source_code", "report_date"],
        suffixes=("_eia", "_cems"),
    )

    # if the heat input reported in cems is >= 90% of the heat input reported in EIA, mark CEMS as the data source, even if we dont have a direct unit to generator match
    fuel_comparison = fuel_comparison[
        fuel_comparison["fuel_consumed_mmbtu_cems"]
        >= (fuel_comparison["fuel_consumed_mmbtu_eia"] * 0.9)
    ]

    # create a data source column
    fuel_comparison["data_source_fuel"] = "cems"

    # drop heat content columns
    fuel_comparison = fuel_comparison.drop(
        columns=["fuel_consumed_mmbtu_eia", "fuel_consumed_mmbtu_cems"]
    )

    # merge this information back into gen_fuel_allocated based on matching plant month fuel
    gen_fuel_allocated = gen_fuel_allocated.merge(
        fuel_comparison,
        how="left",
        on=["plant_id_eia", "energy_source_code", "report_date"],
    )

    # fill missing data source with the fuel-based match
    gen_fuel_allocated["data_source"] = gen_fuel_allocated["data_source"].fillna(
        gen_fuel_allocated["data_source_fuel"]
    )

    # remove intermediate columns
    gen_fuel_allocated = gen_fuel_allocated.drop(columns=["data_source_fuel"])

    # the remainder of the data sources is likely EIA-only
    gen_fuel_allocated["data_source"] = gen_fuel_allocated["data_source"].fillna(
        "eia_only"
    )

    return gen_fuel_allocated


def assign_ba_code_to_plant(df, year):
    """
    Assigns a balancing authority code and state to each plant based on the plant id
    Inputs:
        df: a pandas dataframe containing a 'plant_id_eia' column
        year: four digit year number for the data
    Returns:
        df with a new column for 'ba_code' and 'state'
    """

    pudl_out = load_data.initialize_pudl_out(year=year)

    plant_ba = pudl_out.plants_eia860().loc[
        :,
        [
            "plant_id_eia",
            "balancing_authority_code_eia",
            "balancing_authority_name_eia",
            "utility_name_eia",
            "transmission_distribution_owner_name",
            "state",
        ],
    ]

    # convert the dtype of the balancing authority code column from string to object
    # this will allow for missing values to be filled
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].astype(object)

    # specify a ba code for certain utilities
    utility_as_ba_code = {
        "Anchorage Municipal Light and Power": "AMPL",
        "Arizona Public Service Co": "AZPS",
        "Associated Electric Coop, Inc": "AECI",
        "Avista Corp": "AVA",
        "Avangrid Renewables Inc": "AVRN",
        "Bonneville Power Administration": "BPAT",
        "Bonneville Power Admin": "BPAT",
        "Chugach Electric Assn Inc": "CEA",
        "Duke Energy Carolinas, LLC": "DUK",
        "Duke Energy Florida, Inc": "FPC",
        "Duke Energy Florida, LLC": "FPC",
        "Duke Energy Progress - (NC)": "CPLE",
        "El Paso Electric Co": "EPE",
        "Florida Power & Light Co": "FPL",
        "Florida Power &amp; Light Co": "FPL",
        "Gainesville Regional Utilities": "GVL",
        "Hawaiian Electric Co Inc": "HECO",
        "Hawaii Electric Light Co Inc": "HECO",
        "City of Homestead - (FL)": "HST",
        "Imperial Irrigation District": "IID",
        "JEA": "JEA",
        "Kentucky Utilities Co": "LGEE",
        "Los Angeles Department of Water & Power": "LDWP",
        "Louisville Gas & Electric Co": "LGEE",
        "Nevada Power Co": "NEVP",
        "New Smyrna Beach City of": "NSB",
        "NorthWestern Corporation": "NWMT",
        "NorthWestern Energy": "NWMT",
        "NorthWestern Energy - (SD)": "NWMT",
        "NorthWestern Energy LLC - (MT)": "NWMT",
        "Ohio Valley Electric Corp": "OVEC",
        "Portland General Electric Co": "PGE",
        "Portland General Electric Company": "PGE",
        "PowerSouth Energy Cooperative": "AEC",
        "Public Service Co of Colorado": "PSCO",
        "Public Service Co of NM": "PNM",
        "PUD No 1 of Chelan County": "CHPD",
        "PUD No 1 of Douglas County": "DOPD",
        "PUD No 2 of Grant County": "GCPD",
        "Puget Sound Energy Inc": "PSEI",
        "Sacramento Municipal Util Dist": "BANC",
        "Salt River Project": "SRP",
        "Seminole Electric Cooperative Inc": "SEC",
        "South Carolina Electric&Gas Company": "SCEG",
        "South Carolina Electric & Gas Co": "SCEG",
        "South Carolina Electric &amp; Gas Co": "SCEG",
        "South Carolina Electric&amp;Gas Company": "SCEG",
        "South Carolina Public Service Authority": "SC",
        "South Carolina Public Service Auth": "SC",
        "Southwestern Power Administration": "SPA",
        "Tacoma City of": "TPWR",
        "Tampa Electric Co": "TEC",
        "Tennessee Valley Authority": "TVA",
        "Tucson Electric Power Co": "TEPC",
        "Turlock Irrigation District": "TIDC",
    }

    # fill missing BA codes first based on the BA name, then utility name, then on the transmisison owner name
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["balancing_authority_name_eia"].map(utility_as_ba_code))
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["utility_name_eia"].map(utility_as_ba_code))
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["transmission_distribution_owner_name"].map(utility_as_ba_code))

    # use this to explore plants without an assigned ba
    # sorted(plant_ba[plant_ba['balancing_authority_code_eia'].isna()]['utility_name_eia'].unique().astype(str))

    # rename the ba column
    plant_ba = plant_ba.rename(columns={"balancing_authority_code_eia": "ba_code"})

    # TODO: Remove this once the PUDL issue is fixed
    # As of 4/16/22, there are currently a few incorrect BA assignments in the pudl tables (see https://github.com/catalyst-cooperative/pudl/issues/1584)
    # thus, we will manually correct some of the BA codes based on data in the most recent EIA forms
    manual_ba_corrections = {
        57698: "BANC",
        7966: "SWPP",
        6292: "None",
        7367: "None",
        55966: "None",
        6283: "None",
        57206: "None",
        10093: "None",
    }  # TODO: Tesoro Hawaii has no BA assigned, but is connected to the HECO transmission grid - investigate further

    plant_ba["ba_code"].update(plant_ba["plant_id_eia"].map(manual_ba_corrections))
    plant_ba["ba_code"] = plant_ba["ba_code"].replace("None", np.NaN)

    # add a physical ba code based on the owner of the transmission system
    plant_ba["ba_code_physical"] = plant_ba["ba_code"]
    plant_ba["ba_code_physical"].update(
        plant_ba["transmission_distribution_owner_name"].map(utility_as_ba_code)
    )

    # update based on mapping table when ambiguous
    physical_ba = pd.read_csv("../data/manual/physical_ba.csv")
    plant_ba = plant_ba.merge(
        physical_ba,
        how="left",
        on=["ba_code", "transmission_distribution_owner_name"],
        suffixes=("", "_map"),
    )
    plant_ba["ba_code_physical"].update(plant_ba["ba_code_physical_map"])

    # merge the ba code into the dataframe
    df = df.merge(
        plant_ba.loc[:, ["plant_id_eia", "ba_code", "ba_code_physical", "state"]],
        how="left",
        on="plant_id_eia",
    )

    return df


def create_monthly_gens_records(df, year):
    """
    Creates a duplicate record for each month of the year in the gens file
    """
    # If we want to allocate net generation at the monthly level, we need to ensure that the gens file has monthly records
    # to do this, we can duplicate the records in gens 11 times for each month, so that there is a record for each month of the year
    # duplicate the entries for each month

    if "report_date" not in df.columns:
        # create a report date column with the first month of the year
        df["report_date"] = f"{year}-01-01"
        df["report_date"] = pd.to_datetime(df["report_date"])

    df_month = df.copy()

    month = 2
    while month <= 12:
        # add one month to the copied data each iteration
        df_month["report_date"] = df_month["report_date"] + pd.DateOffset(months=1)
        # concat this data to the gens file
        df = pd.concat([df, df_month], axis=0)
        month += 1

    return df


def clean_cems(year):
    """
    Coordinating function for all of the cems data cleaning
    """
    # load the CEMS data
    cems = load_data.load_cems_data(year)

    # remove non-grid connected plants
    cems = remove_plants(
        cems,
        non_grid_connected=True,
        remove_states=["PR"],
        steam_only_plants=False,
        distribution_connected_plants=False,
    )

    # manually remove steam-only units
    cems = manually_remove_steam_units(cems)

    # remove plants that only report steam generation and no electrical generation
    # NOTE: keeping steam only plants for now
    # cems = remove_heating_only_plants(cems)

    # add a report date
    cems = add_report_date(cems)

    # identify cems reporting status (full year or partial year)
    # NOTE: this information is not really useful yet, so we are not going to run this to save time
    # cems = determine_cems_reporting_status(cems)

    # TODO: identify and remove any hourly values that appear to be outliers

    # add a fuel type to each observation
    cems = assign_fuel_type_to_cems(cems, year)

    # fill in missing hourly emissions data using the fuel type and heat input
    cems = fill_cems_missing_co2(cems, year)

    # remove any observations from cems where zero operation is reported for an entire month
    # although this data could be considered to be accurately reported, let's remove it so that we can double check against the eia data
    # TODO: check if any of these observations are from geothermal generators
    cems = remove_cems_with_zero_monthly_data(cems)

    # calculated CHP-adjusted emissions
    cems = adjust_cems_for_CHP(cems)

    # identify any remaining missing values
    # TODO: Try to identify fuel types
    still_missing_co2_data = list(
        cems[cems["co2_mass_tons"].isnull()]["cems_id"].unique()
    )
    print(
        f"Unable to calculate emissions for the following plants_units: {still_missing_co2_data}"
    )

    # For now, lets drop these from the data
    cems = cems[~cems["cems_id"].isin(still_missing_co2_data)]

    return cems


def adjust_cems_for_CHP(cems, drop_interim_columns=True):
    """
    Calculates the portion of fuel consumption and CO2 emissions for electricity for each hour in CEMS.
    """
    # factors to convert to MMBTU
    mwh_to_mmbtu = 3.412142
    klb_to_mmbtu = 1.194  # NOTE: this might differ for each plant

    # calculate total heat output
    cems["heat_output_mmbtu"] = (cems["gross_generation_mwh"] * mwh_to_mmbtu) + (
        cems["steam_load_1000_lbs"] * klb_to_mmbtu
    )

    # calculate the fraction of heat input for electricity
    cems["frac_electricity"] = (cems["gross_generation_mwh"] * mwh_to_mmbtu) / cems[
        "heat_output_mmbtu"
    ]
    # where both of these terms are zero, change the fraction to 1
    cems.loc[
        (cems.gross_generation_mwh == 0) & (cems.heat_output_mmbtu == 0),
        "frac_electricity",
    ] = cems.loc[
        (cems.gross_generation_mwh == 0) & (cems.heat_output_mmbtu == 0),
        "frac_electricity",
    ].fillna(
        1
    )

    # calculate fuel consumed for electricity and co2 adjusted
    cems["fuel_consumed_for_electricity_mmbtu"] = (
        cems["fuel_consumed_mmbtu"] * cems["frac_electricity"]
    )
    cems["co2_mass_tons_adjusted"] = cems["co2_mass_tons"] * cems["frac_electricity"]

    if drop_interim_columns:
        cems = cems.drop(columns=["heat_output_mmbtu", "frac_electricity"])

    return cems


def model_gross_to_net(df):
    """
    Create a linear regression model of monthly gross to net generation

    Args:
        arg
    Returns:
        output
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # get a linear model for the data points
        model = smf.ols("net_generation_mwh ~ gross_generation_mwh", data=df).fit()

        # find and remove any outliers
        try:
            outliers = model.outlier_test()
            corrected = df[~df.index.isin(outliers[outliers["bonf(p)"] < 0.5].index)]

            # get a linear model of the corrected data
            model = smf.ols(
                "net_generation_mwh ~ gross_generation_mwh", data=corrected
            ).fit()
        except ValueError:
            pass
        slope = model.params[1]
        rsquared = model.rsquared
        rsquared_adj = model.rsquared_adj
        number_observations = model.nobs

    return slope, rsquared, rsquared_adj, number_observations


def convert_gross_to_net_generation(cems, gen_fuel_allocated):
    """
    Converts hourly gross generation in CEMS to hourly net generation by calculating a gross to net generation ratio
    Inputs:

    Returns: 
        cems df with an added column for net_generation_mwh and a column indicated the method used to calculate net generation
    """

    # add a placeholder column that assumes a 1:1 gross to net generation ratio
    # if for some reason we are not able to calculate a gross to net generation ratio, this will be used as the default assumption
    cems["net_generation_mwh"] = cems["gross_generation_mwh"]

    # load the allocated eia data for each month where there is corresponding cems data
    eia_plant_month_net_gen = gen_fuel_allocated[
        (gen_fuel_allocated["data_source"] == "cems")
        & ~(gen_fuel_allocated["net_generation_mwh"].isna())
    ]
    # aggregate at the plant month level
    eia_plant_month_net_gen = (
        eia_plant_month_net_gen.groupby(["plant_id_eia", "report_date"])
        .sum()["net_generation_mwh"]
        .reset_index()
    )

    # calculate the total gross generation for each plant month in cems
    cems_plant_month_gross_gen = (
        cems.groupby(["plant_id_eia", "report_date"])
        .sum()["gross_generation_mwh"]
        .reset_index()
    )

    # merge the net generation data into the gross generation data
    monthly_gtn_ratio = cems_plant_month_gross_gen.merge(
        eia_plant_month_net_gen, how="left", on=["plant_id_eia", "report_date"]
    )

    # calculate the gtn
    monthly_gtn_ratio["gross_to_net_ratio"] = (
        monthly_gtn_ratio["net_generation_mwh"]
        / monthly_gtn_ratio["gross_generation_mwh"]
    )

    # only keep values where the monthly ratio is greater than zero
    monthly_gtn_ratio.loc[
        (monthly_gtn_ratio["gross_to_net_ratio"] < 0), "gross_to_net_ratio"
    ] = np.NaN

    # Set up the regression analysis for missing values

    # only keep values where there are not missing values
    gtn_regression = monthly_gtn_ratio.copy()[
        ~(monthly_gtn_ratio["gross_to_net_ratio"].isna())
    ]
    # calculate the ratio for each plant and create a dataframe
    gtn_regression = gtn_regression.groupby("plant_id_eia").apply(model_gross_to_net)
    gtn_regression = pd.DataFrame(
        gtn_regression.tolist(),
        index=gtn_regression.index,
        columns=["gtn_linear", "rsquared", "rsquared_adj", "observations"],
    ).reset_index()
    # only keep the results with adjusted rsquared values greater than 0.70
    gtn_regression = gtn_regression[gtn_regression["rsquared_adj"] >= 0.7]

    # merge in regression results
    monthly_gtn_ratio = monthly_gtn_ratio.merge(
        gtn_regression[["plant_id_eia", "gtn_linear"]], how="left", on="plant_id_eia"
    )

    # add a status column for how the net generation was calculated
    monthly_gtn_ratio["net_gen_method"] = "monthly_ratio"
    monthly_gtn_ratio.loc[
        (monthly_gtn_ratio["gross_to_net_ratio"].isna())
        & ~(monthly_gtn_ratio["gtn_linear"].isna()),
        "net_gen_method",
    ] = "annual_regression"
    monthly_gtn_ratio.loc[
        (monthly_gtn_ratio["gross_to_net_ratio"].isna())
        & (monthly_gtn_ratio["gtn_linear"].isna()),
        "net_gen_method",
    ] = "net_equals_gross"

    # fill missing values using the ratio from the regression results
    monthly_gtn_ratio["gross_to_net_ratio"] = monthly_gtn_ratio[
        "gross_to_net_ratio"
    ].fillna(monthly_gtn_ratio["gtn_linear"])

    # merge the gtn ratio into the cems data
    cems = cems.merge(
        monthly_gtn_ratio[
            ["plant_id_eia", "report_date", "gross_to_net_ratio", "net_gen_method"]
        ],
        how="left",
        on=["plant_id_eia", "report_date"],
    )
    # calculate hourly net generation
    cems["net_generation_mwh_calculated"] = (
        cems["gross_generation_mwh"] * cems["gross_to_net_ratio"]
    )

    # update the net generation column using the calculated values
    cems["net_generation_mwh"].update(cems["net_generation_mwh_calculated"])

    # update the method column to indicate which used the default assumption
    cems["net_gen_method"] = cems["net_gen_method"].fillna("net_equals_gross")

    # drop the calculated column
    cems = cems.drop(columns=["net_generation_mwh_calculated"])

    return cems


def impute_missing_hourly_net_generation(cems, gen_fuel_allocated):
    """
    Where CEMS reports hourly heat content but no hourly generation, and EIA reports that there is positive net generation,
    impute the missing net generation values based on the hourly heat input.
    We calculate "heat rates" for both steam load and electricity load from EIA data
    We then adjust the hourly heat data to account for steam load

    NOTE: 4/21/22: this has been modified to only inpute net generation for plants that also do not report steam input
    In some cases, these are actually just steam only generators
    Unless we match on units using the crosswalk, we should not use this to impute data
    """
    # calculate total values for each plant-month in CEMS and EIA
    cems_monthly = (
        cems.groupby(["plant_id_eia", "report_date"])
        .sum()[
            [
                "gross_load_mw",
                "gross_generation_mwh",
                "steam_load_1000_lbs",
                "fuel_consumed_mmbtu",
                "co2_mass_tons",
            ]
        ]
        .reset_index()
    )
    eia_monthly = (
        gen_fuel_allocated.groupby(["plant_id_eia", "report_date"]).sum().reset_index()
    )

    # identify all plant months that report no generation but have heat input and have no steam production
    missing_generation = cems_monthly[
        (cems_monthly["gross_generation_mwh"] == 0)
        & (cems_monthly["fuel_consumed_mmbtu"] > 0)
        & (cems_monthly["steam_load_1000_lbs"] == 0)
    ]
    # merge in the EIA net generation data
    missing_generation = missing_generation.merge(
        eia_monthly[
            [
                "plant_id_eia",
                "report_date",
                "net_generation_mwh",
                "fuel_consumed_mmbtu",
                "fuel_consumed_for_electricity_mmbtu",
            ]
        ],
        how="left",
        on=["plant_id_eia", "report_date"],
        suffixes=("_cems", "_eia"),
    )

    # only keep plant-months where zero gross generation was reported, but EIA reports positive net generation
    # ignore negative net generation for now
    missing_generation = missing_generation[
        missing_generation["net_generation_mwh"] > 0
    ]

    # scale the EIA-reported heat input by the heat input reported in CEMS
    missing_generation["heat_scaling_factor"] = (
        missing_generation["fuel_consumed_mmbtu_cems"]
        / missing_generation["fuel_consumed_mmbtu_eia"]
    )

    # calculate a heat rate for steam
    # missing_generation['steam_heat_rate_mmbtu_per_klb'] = (missing_generation['fuel_consumed_mmbtu'] - missing_generation['fuel_consumed_for_electricity_mmbtu']) / missing_generation['steam_load_1000_lbs'] * missing_generation['heat_scaling_factor']
    # missing_generation['steam_heat_rate_mmbtu_per_klb'] = missing_generation['steam_heat_rate_mmbtu_per_klb'].fillna(0)

    # calculate a heat rate for electricity
    missing_generation["heat_to_netgen_mwh_per_mmbtu"] = missing_generation[
        "net_generation_mwh"
    ] / (
        missing_generation["fuel_consumed_for_electricity_mmbtu"]
        * missing_generation["heat_scaling_factor"]
    )

    # get a list of all of the plants that have missing net gen that needs to be imputed
    plants_missing_net_gen = list(missing_generation.plant_id_eia.unique())

    cems_ng_imputation = cems[cems["plant_id_eia"].isin(plants_missing_net_gen)]

    # save the original index to assist with updating the original data later
    cems_index = cems_ng_imputation.index

    # merge in the factors
    # cems_ng_imputation = cems_ng_imputation.merge(missing_generation[['plant_id_eia','report_date','steam_heat_rate_mmbtu_per_klb','heat_to_netgen_mwh_per_mmbtu']], how='left', on=['plant_id_eia','report_date'])
    cems_ng_imputation = cems_ng_imputation.merge(
        missing_generation[
            ["plant_id_eia", "report_date", "heat_to_netgen_mwh_per_mmbtu"]
        ],
        how="left",
        on=["plant_id_eia", "report_date"],
    )

    # calculate hourly heat content for electricity
    # cems_ng_imputation['heat_content_for_electricity_mmbtu'] = cems_ng_imputation['fuel_consumed_mmbtu'] - (cems_ng_imputation['steam_load_1000_lbs'].fillna(0) * cems_ng_imputation['steam_heat_rate_mmbtu_per_klb'])

    # convert heat input to net generation
    # cems_ng_imputation['net_generation_mwh'] = cems_ng_imputation['heat_content_for_electricity_mmbtu'] * cems_ng_imputation['heat_to_netgen_mwh_per_mmbtu']
    cems_ng_imputation["net_generation_mwh"] = (
        cems_ng_imputation["fuel_consumed_mmbtu"]
        * cems_ng_imputation["heat_to_netgen_mwh_per_mmbtu"]
    )

    # drop intermediate columns
    # cems_ng_imputation = cems_ng_imputation.drop(columns=['steam_heat_rate_mmbtu_per_klb','heat_to_netgen_mwh_per_mmbtu'])
    cems_ng_imputation = cems_ng_imputation.drop(
        columns=["heat_to_netgen_mwh_per_mmbtu"]
    )

    # if there are missing values where heat content is also zero, fill with zero
    cems_ng_imputation.loc[
        cems_ng_imputation["fuel_consumed_mmbtu"] == 0, "net_generation_mwh"
    ] = cems_ng_imputation.loc[
        cems_ng_imputation["fuel_consumed_mmbtu"] == 0, "net_generation_mwh"
    ].fillna(
        0
    )

    # change the net gen method
    cems_ng_imputation["net_gen_method"] = "imputed_from_fuel_consumption"
    cems_ng_imputation.loc[
        cems_ng_imputation["net_generation_mwh"].isna(), "net_gen_method"
    ] = np.NaN

    # add the original index and update
    cems_ng_imputation.index = cems_index
    cems["net_generation_mwh"].update(cems_ng_imputation["net_generation_mwh"])
    cems["net_gen_method"].update(cems_ng_imputation["net_gen_method"])

    return cems


def clean_eia_930(df: DataFrame):
    """
    Args:
       df (pd.DataFrame): dataframe containing rows of EIA-930 in the format provided by balance
       sheets.
    Returns:
       cleaned df with same format as input
    """
    ## Remove bad data (negative and zero fossil fuel generation)
    fossil_cols = [
        "Net Generation (MW) from Coal",
        "Net Generation (MW) from Natural Gas",
        "Net Generation (MW) from All Petroleum Products",
    ]
    for col in fossil_cols:
        df[df[col] < 0] = np.nan

    # TODO other forms of cleaning as needed

    return df


def add_report_date(df):
    """
    Add a report date column to the cems data based on the plant's local timezone

    Args:
        df (pd.Dataframe): dataframe containing 'plant_id_eia' and 'operating_datetime_utc' columns
    Returns:
        Original dataframe with 'report_date' column added
    """
    plants_entity_eia = load_data.load_pudl_table("plants_entity_eia")

    # get timezone
    df = df.merge(
        plants_entity_eia[["plant_id_eia", "timezone"]], how="left", on="plant_id_eia"
    )

    # create a datetimeindex from the operating_datetime_utc column
    datetime_utc = pd.DatetimeIndex(df["operating_datetime_utc"])

    # create blank column to hold local datetimes
    df["report_date"] = np.NaN

    # get list of unique timezones
    timezones = list(df["timezone"].unique())

    # convert UTC to the local timezone
    for tz in timezones:
        tz_mask = df["timezone"] == tz  # find all rows where the tz matches
        df.loc[tz_mask, "report_date"] = (
            datetime_utc[tz_mask]
            .tz_convert(tz)  # convert to local time
            .to_series(index=df[tz_mask].index)  # convert to a series
            .dt.to_period("M")
            .dt.to_timestamp()  # convert to a YYYY-MM-01 stamp
        )

    df["report_date"] = pd.to_datetime(df["report_date"])

    # drop the operating_datetime_local column
    df = df.drop(columns=["timezone"])

    return df


def ba_timezone(ba, type):
    """
    Retrieves the UTC Offset (for standard time) for a single balancing area.
    Args:
        ba: string containing the ba_code
        type: either 'reporting_eia930' or 'local'. Reporting will return the TZ used by the BA when reporting to EIA-930, local will return the actual local tz
    """

    tz = pd.read_csv(
        "../data/manual/ba_reference.csv", usecols=["ba_code", f"timezone_{type}"]
    )
    tz = tz.loc[tz["ba_code"] == ba, f"timezone_{type}"].item()

    return tz


def distribute_monthly_eia_data_to_hourly(
    monthly_eia_data_to_distribute, hourly_profiles, profile_column_name
):
    """
    Uses monthly-level EIA data and assigns an hourly profile
    Inputs: 
        monthly_eia_data_to_distribute: a dataframe that contains monthly total net generation, fuel consumption, and co2 data, along with columns for report_date and ba_code
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


def create_primary_fuel_table(gen_fuel_allocated):
    """
    Identifies the primary fuel for each generator and plant
    Gen primary fuel is identified based on the "energy source code 1" identified in EIA-860
    Plant primary fuel is based on the most-consumed fuel at a plant based on allocated heat input
    """
    # get a table of primary energy source codes
    gen_primary_fuel = gen_fuel_allocated[
        gen_fuel_allocated["energy_source_code_num"] == "energy_source_code_1"
    ].drop_duplicates(subset=["plant_id_eia", "generator_id"])[
        ["plant_id_eia", "generator_id", "energy_source_code"]
    ]
    # rename the energy source code column to gen primary fuel
    # gen_primary_fuel = gen_primary_fuel.rename(columns={'energy_source_code':'generator_primary_fuel'})

    # calculate the total annual heat input by fuel type for each plant
    plant_primary_fuel = (
        gen_fuel_allocated.groupby(["plant_id_eia", "energy_source_code"])
        .sum()[["fuel_consumed_mmbtu"]]
        .reset_index()
    )

    # identify the energy source code with the greatest fuel consumption for each plant
    plant_primary_fuel = plant_primary_fuel[
        plant_primary_fuel.groupby("plant_id_eia")["fuel_consumed_mmbtu"].transform(max)
        == plant_primary_fuel["fuel_consumed_mmbtu"]
    ][["plant_id_eia", "energy_source_code"]]

    # rename the column to plant primary fuel
    plant_primary_fuel = plant_primary_fuel.rename(
        columns={"energy_source_code": "plant_primary_fuel"}
    )

    # merge the plant primary fuel into the gen primary fuel
    primary_fuel_table = gen_primary_fuel.merge(
        plant_primary_fuel, how="left", on="plant_id_eia"
    )

    return primary_fuel_table


def identify_distribution_connected_plants(df, year, voltage_threshold_kv=60):
    """
    Identifies which plant_id_eia are "distribution grid connected" based on a voltage threshold.

    The distribution grid is generally considered to operate at 60-69kV and below.
    Thus, any plants that have a grid voltage under this threshold will be flagged as distribution connected.
    Args:
        df: pandas dataframe with a column for plant_id_eia
        voltage_threshold_kv: the voltage (kV) under which a plant will be considered to be a distribution asset
    Returns:
        df: with additional binary column `distribution_flag`
    """
    # load the EIA-860 data
    pudl_out = load_data.initialize_pudl_out(year=year)

    plant_voltage = pudl_out.plants_eia860().loc[:, ["plant_id_eia", "grid_voltage_kv"]]

    plant_voltage = plant_voltage.assign(
        distribution_flag=lambda x: np.where(
            x.grid_voltage_kv <= voltage_threshold_kv, True, False
        )
    )

    df = df.merge(
        plant_voltage[["plant_id_eia", "distribution_flag"]],
        how="left",
        on="plant_id_eia",
    )

    return df


def assign_fuel_category_to_ESC(
    df, fuel_category_name, esc_column="energy_source_code"
):
    """
    Assigns a fuel category to each energy source code in a dataframe.
    Args:
        df: pandas dataframe with column name that matches fuel_category_name and contains energy source codes
        fuel_category_name: name of the column in energy_source_groups.csv that contains the desired category mapping
        esc_column: name of the column in df that contains the energy source codes to assign a category to
    Returns:
        df with additional column for fuel category
    """
    # load the fuel category table
    energy_source_groups = pd.read_csv("../data/manual/energy_source_groups.csv")[
        ["energy_source_code", fuel_category_name]
    ].rename(
        columns={"energy_source_code": esc_column, fuel_category_name: "fuel_category"}
    )
    # assign a fuel category to the monthly eia data
    df = df.merge(
        energy_source_groups[[esc_column, "fuel_category"]], how="left", on=esc_column
    )

    return df


def clean_eia923(year,
                 include_nox=False,
                 include_so2=False,
                 include_co2e=False):
    """
    This is the coordinating function for cleaning and allocating generation and fuel data in EIA-923.
    """
    # Distribute net generation and heat input data reported by the three different EIA-923 tables
    pudl_out = load_data.initialize_pudl_out(year)

    # allocate net generation and heat input to each generator-fuel grouping
    gen_fuel_allocated = allocate_gen_fuel.allocate_gen_fuel_by_generator_energy_source(
        pudl_out, drop_interim_cols=True
    )

    # create a table that identifies the primary fuel of each generator and plant
    primary_fuel_table = create_primary_fuel_table(gen_fuel_allocated)

    # calculate co2 emissions for each generator-fuel based on allocated fuel consumption
    gen_fuel_allocated = calculate_co2_from_fuel_consumption(gen_fuel_allocated, year)

    cols_to_sum_by_generator = [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_tons",
        "co2_mass_tons_adjusted"]

    if include_nox:
        gen_fuel_allocated = calculate_nox_from_fuel_consumption(gen_fuel_allocated)
        cols_to_sum_by_generator.append("nox_mass_lbs")
        if "fuel_consumed_for_electricity_units" in gen_fuel_allocated:
            cols_to_sum_by_generator.append("nox_mass_lbs_adjusted")

    if include_so2:
        gen_fuel_allocated = calculate_so2_from_fuel_consumption(gen_fuel_allocated)
        cols_to_sum_by_generator.append("so2_mass_lbs")
        if "fuel_consumed_for_electricity_units" in gen_fuel_allocated:
            cols_to_sum_by_generator.append("so2_mass_lbs_adjusted")
    
    if include_co2e:
        gen_fuel_allocated = calculate_co2e_from_fuel_consumption(gen_fuel_allocated, year)
        cols_to_sum_by_generator.append("co2e_mass_tons")
        if "fuel_consumed_for_electricity_units" in gen_fuel_allocated:
            cols_to_sum_by_generator.append("co2e_mass_tons_adjusted")

    # aggregate the allocated data to the generator level
    gen_fuel_allocated = allocate_gen_fuel.agg_by_generator(
        gen_fuel_allocated,
        sum_cols=cols_to_sum_by_generator
    )

    # remove any plants that we don't want in the data
    gen_fuel_allocated = remove_plants(
        gen_fuel_allocated,
        non_grid_connected=True,
        remove_states=["PR"],
        steam_only_plants=False,
        distribution_connected_plants=False,
    )

    return gen_fuel_allocated, primary_fuel_table
