import pandas as pd
import numpy as np
import warnings

import src.load_data as load_data

import pudl.analysis.allocate_net_gen as allocate_gen_fuel

from src.column_checks import get_dtypes, apply_dtypes


def clean_eia923(year, small):
    """
    This is the coordinating function for cleaning and allocating generation and fuel data in EIA-923.
    """
    # Distribute net generation and heat input data reported by the three different EIA-923 tables

    pudl_out = load_data.initialize_pudl_out(year=year)

    # allocate net generation and heat input to each generator-fuel grouping
    gen_fuel_allocated = allocate_gen_fuel.allocate_gen_fuel_by_generator_energy_source(
        pudl_out, drop_interim_cols=True
    )

    # manually update energy source code when OTH
    gen_fuel_allocated = update_energy_source_codes(gen_fuel_allocated)

    # round all values to the nearest tenth of a unit
    gen_fuel_allocated.loc[
        :,
        [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ],
    ] = gen_fuel_allocated.loc[
        :,
        [
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ],
    ].round(
        1
    )

    # create a table that identifies the primary fuel of each generator and plant
    primary_fuel_table = create_primary_fuel_table(gen_fuel_allocated, pudl_out)

    if small:
        gen_fuel_allocated = smallerize_test_data(df=gen_fuel_allocated, random_seed=42)

    # calculate co2 emissions for each generator-fuel based on allocated fuel consumption
    gen_fuel_allocated = calculate_ghg_emissions_from_fuel_consumption(
        df=gen_fuel_allocated,
        year=year,
        include_co2=True,
        include_ch4=True,
        include_n2o=True,
    )

    # Calculate NOx and SO2 emissions
    gen_fuel_allocated = calculate_nox_from_fuel_consumption(
        gen_fuel_allocated, pudl_out, year
    )
    gen_fuel_allocated = calculate_so2_from_fuel_consumption(
        gen_fuel_allocated, pudl_out, year
    )

    # adjust emissions for CHP
    gen_fuel_allocated = adjust_emissions_for_CHP(gen_fuel_allocated)

    # adjust emissions for biomass
    gen_fuel_allocated = adjust_emissions_for_biomass(gen_fuel_allocated)

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

    # aggregate the allocated data to the generator level
    gen_fuel_allocated = allocate_gen_fuel.agg_by_generator(
        gen_fuel_allocated, sum_cols=DATA_COLUMNS,
    )

    # remove any plants that we don't want in the data
    gen_fuel_allocated = remove_plants(
        gen_fuel_allocated,
        non_grid_connected=True,
        remove_states=["PR"],
        steam_only_plants=False,
        distribution_connected_plants=False,
    )

    # round all values to the nearest tenth of a unit
    gen_fuel_allocated.loc[:, DATA_COLUMNS] = gen_fuel_allocated.loc[
        :, DATA_COLUMNS
    ].round(1)

    # add subplant id
    subplant_crosswalk = pd.read_csv(
        f"../data/outputs/{year}/subplant_crosswalk.csv", dtype=get_dtypes()
    )[["plant_id_eia", "generator_id", "subplant_id"]].drop_duplicates()
    gen_fuel_allocated = gen_fuel_allocated.merge(
        subplant_crosswalk, how="left", on=["plant_id_eia", "generator_id"]
    )

    # add the cleaned prime mover code to the data
    gen_pm = pudl_out.gens_eia860()[
        ["plant_id_eia", "generator_id", "prime_mover_code"]
    ]
    gen_fuel_allocated = gen_fuel_allocated.merge(
        gen_pm, how="left", on=["plant_id_eia", "generator_id"]
    )

    gen_fuel_allocated = apply_dtypes(gen_fuel_allocated)
    primary_fuel_table = apply_dtypes(primary_fuel_table)

    return gen_fuel_allocated, primary_fuel_table


def update_energy_source_codes(df):
    """
    Manually update fuel source codes
    """
    # refinery with energy source = OTH, change to OG
    df.loc[
        (df["plant_id_eia"] == 50626) & (df["energy_source_code"] == "OTH"),
        "energy_source_code",
    ] = "OG"
    df.loc[
        (df["plant_id_eia"] == 56139) & (df["energy_source_code"] == "OTH"),
        "energy_source_code",
    ] = "OG"
    df.loc[
        (df["plant_id_eia"] == 59073) & (df["energy_source_code"] == "OTH"),
        "energy_source_code",
    ] = "OG"

    return df


def create_primary_fuel_table(gen_fuel_allocated, pudl_out):
    """
    Identifies the primary fuel for each generator and plant
    Gen primary fuel is identified based on the "energy source code 1" identified in EIA-860
    Plant primary fuel is based on the most-consumed fuel at a plant based on allocated heat input
    """

    primary_fuel_from_capacity = calculate_capacity_based_primary_fuel(pudl_out)

    # get a table of primary energy source codes
    gen_primary_fuel = gen_fuel_allocated[
        gen_fuel_allocated["energy_source_code_num"] == "energy_source_code_1"
    ].drop_duplicates(subset=["plant_id_eia", "generator_id"])[
        ["plant_id_eia", "generator_id", "energy_source_code"]
    ]

    # create a blank dataframe with all of the plant ids to hold primary fuel data
    plant_primary_fuel = gen_fuel_allocated[["plant_id_eia"]].drop_duplicates()

    # calculate the total annual fuel consumption, generation, and capacity by fuel type
    #  for each plant
    plant_totals_by_fuel = (
        gen_fuel_allocated.groupby(["plant_id_eia", "energy_source_code"], dropna=False)
        .sum()[["fuel_consumed_mmbtu", "net_generation_mwh"]]
        .reset_index()
    )

    # we will calculate primary fuel based on the fuel with the most consumption,
    # generation, and capacity
    for source in ["fuel_consumed_mmbtu", "net_generation_mwh"]:

        # only keep values greater than zero so that these can be filled by other
        # methods if non-zero
        primary_fuel_calc = plant_totals_by_fuel[plant_totals_by_fuel[source] > 0]

        # identify the fuel type with the maximum value for each plant
        primary_fuel_calc = primary_fuel_calc[
            primary_fuel_calc.groupby("plant_id_eia", dropna=False)[source].transform(
                max
            )
            == primary_fuel_calc[source]
        ][["plant_id_eia", "energy_source_code"]]

        # remove duplicate values (if two fuels are both the maximum)
        primary_fuel_calc = primary_fuel_calc.drop_duplicates(
            subset="plant_id_eia", keep=False
        )
        primary_fuel_calc = primary_fuel_calc.rename(
            columns={"energy_source_code": f"primary_fuel_from_{source}"}
        )

        # merge the primary fuel into the main table
        plant_primary_fuel = plant_primary_fuel.merge(
            primary_fuel_calc, how="left", on="plant_id_eia", validate="1:1"
        )

    # merge the primary fuel into the main table
    plant_primary_fuel = plant_primary_fuel.merge(
        primary_fuel_from_capacity, how="left", on="plant_id_eia", validate="1:1"
    )

    # use the fuel-based primary fuel first, then fill using capacit-based primary fuel,
    # then generation based.
    plant_primary_fuel["plant_primary_fuel"] = plant_primary_fuel[
        "primary_fuel_from_fuel_consumed_mmbtu"
    ]
    plant_primary_fuel["plant_primary_fuel"] = plant_primary_fuel[
        "plant_primary_fuel"
    ].fillna(plant_primary_fuel["primary_fuel_from_capacity_mw"])
    plant_primary_fuel["plant_primary_fuel"] = plant_primary_fuel[
        "plant_primary_fuel"
    ].fillna(plant_primary_fuel["primary_fuel_from_net_generation_mwh"])

    if len(plant_primary_fuel[plant_primary_fuel["plant_primary_fuel"].isna()]) > 0:
        raise UserWarning(
            "Plant primary fuel table contains missing primary fuels. Update method of `create_primary_fuel_table()` to fix"
        )

    # merge the plant primary fuel into the gen primary fuel
    primary_fuel_table = gen_primary_fuel.merge(
        plant_primary_fuel[["plant_id_eia", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
        validate="many_to_one",
    )

    return primary_fuel_table


def calculate_capacity_based_primary_fuel(pudl_out):
    # create a table of primary fuel by nameplate capacity
    gen_capacity = pudl_out.gens_eia860().loc[
        :, ["plant_id_eia", "generator_id", "capacity_mw", "energy_source_code_1"]
    ]

    gen_capacity = (
        gen_capacity.groupby(["plant_id_eia", "energy_source_code_1"], dropna=False)
        .sum()["capacity_mw"]
        .reset_index()
    )

    # drop the battery portion of any hybrid plants so that we don't accidentally
    # identify the primary fuel as storage
    gen_capacity = gen_capacity[
        ~(
            (gen_capacity.duplicated(subset="plant_id_eia", keep=False))
            & (gen_capacity.energy_source_code_1 == "MWH")
        )
    ]

    # find the fuel with the greatest capacity
    gen_capacity = gen_capacity[
        gen_capacity.groupby("plant_id_eia", dropna=False)["capacity_mw"].transform(max)
        == gen_capacity["capacity_mw"]
    ][["plant_id_eia", "energy_source_code_1"]].rename(
        columns={"energy_source_code_1": "primary_fuel_from_capacity_mw"}
    )

    # drop any duplicate entries (if two fuel types have the same nameplate capacity)
    gen_capacity = gen_capacity[
        ~(gen_capacity.duplicated(subset="plant_id_eia", keep=False))
    ]

    return gen_capacity


def calculate_ghg_emissions_from_fuel_consumption(
    df, year, include_co2=True, include_ch4=True, include_n2o=True
):
    """
    Inputs:
        df: pandas dataframe containing the following columns: ['plant_id_eia', 'report_date,'fuel_consumed_mmbtu','energy_source_code']
    """

    emissions_to_calc = []
    if include_co2 is True:
        emissions_to_calc.append("co2")
    if include_ch4 is True:
        emissions_to_calc.append("ch4")
    if include_n2o is True:
        emissions_to_calc.append("n2o")

    efs_to_use = [emission + "_lb_per_mmbtu" for emission in emissions_to_calc]

    # get emission factors
    emission_factors = load_data.load_ghg_emission_factors()[
        ["energy_source_code"] + efs_to_use
    ]

    # add emission factor to  df
    df = df.merge(emission_factors, how="left", on="energy_source_code")

    # if there are any geothermal units, load the geothermal EFs
    if df["energy_source_code"].str.contains("GEO").any():
        df = add_geothermal_emission_factors(
            df, year, include_co2=True, include_nox=False, include_so2=False
        )

    # create a new column with the emissions mass
    for e in emissions_to_calc:
        df[f"{e}_mass_lb"] = df["fuel_consumed_mmbtu"] * df[f"{e}_lb_per_mmbtu"]

    # drop intermediate columns
    df = df.drop(columns=efs_to_use)

    return df


def add_geothermal_emission_factors(
    df, year, include_co2=True, include_nox=True, include_so2=True
):
    """"""

    emissions_to_calc = []
    if include_co2 is True:
        emissions_to_calc.append("co2")
    if include_nox is True:
        emissions_to_calc.append("nox")
    if include_so2 is True:
        emissions_to_calc.append("so2")

    efs_to_use = [emission + "_lb_per_mmbtu" for emission in emissions_to_calc]

    geothermal_efs = calculate_geothermal_emission_factors(year).loc[
        :, ["plant_id_eia", "generator_id", "plant_frac"] + efs_to_use
    ]

    for e in emissions_to_calc:
        geothermal_efs = geothermal_efs.rename(
            columns={f"{e}_lb_per_mmbtu": f"{e}_lb_per_mmbtu_geo"}
        )

    # if there is a merge key for generator id, merge in the geothermal EFs on generator id
    if "generator_id" in list(df.columns):
        # add geothermal emission factor to df
        df = df.merge(
            geothermal_efs.drop(columns=["plant_frac"]),
            how="left",
            on=["plant_id_eia", "generator_id"],
        )
    # otherwise, aggregate EF to plant level and merge
    else:
        # multiply the emission factor by the fraction
        for e in emissions_to_calc:
            geothermal_efs[f"{e}_lb_per_mmbtu_geo"] = (
                geothermal_efs["plant_frac"] * geothermal_efs[f"{e}_lb_per_mmbtu_geo"]
            )
        # groupby plant to get the weighted emission factor
        geothermal_efs = (
            geothermal_efs.groupby("plant_id_eia", dropna=False).sum().reset_index()
        ).drop(columns=["plant_frac"])
        # add geothermal emission factor to df
        df = df.merge(geothermal_efs, how="left", on=["plant_id_eia"])

    # update missing efs using the geothermal efs if available
    for e in emissions_to_calc:
        if f"{e}_lb_per_mmbtu" not in df.columns:
            df[f"{e}_lb_per_mmbtu"] = np.NaN
        df[f"{e}_lb_per_mmbtu"] = df[f"{e}_lb_per_mmbtu"].fillna(
            df[f"{e}_lb_per_mmbtu_geo"]
        )

    # drop intermediate columns
    for e in emissions_to_calc:
        df = df.drop(columns=[f"{e}_lb_per_mmbtu_geo"])

    return df


def calculate_geothermal_emission_factors(year):
    """
    Updates the list of geothermal plants provided by EPA using EIA data
    Calculates a weighted average EF for each plant-month based on the fraction
    of fuel consumed from each type of prime mover (steam, binary, flash)
    """
    # load geothermal efs
    geothermal_efs = pd.read_csv(
        "../data/manual/egrid_static_tables/table_C6_geothermal_emission_factors.csv",
        dtype=get_dtypes(),
    ).loc[
        :, ["geotype_code", "co2_lb_per_mmbtu", "nox_lb_per_mmbtu", "so2_lb_per_mmbtu"]
    ]

    geothermal_geotypes = identify_geothermal_generator_geotype(year)

    # merge in the emission factor
    geo_efs = geothermal_geotypes.merge(geothermal_efs, how="left", on="geotype_code")

    return geo_efs


def identify_geothermal_generator_geotype(year):
    """Identifies whether each geothermal generator is binary, flash, or dry steam"""
    pudl_out = load_data.initialize_pudl_out(year)

    geothermal_geotype = pudl_out.gens_eia860()
    geothermal_geotype = geothermal_geotype.loc[
        geothermal_geotype["energy_source_code_1"] == "GEO",
        [
            "plant_id_eia",
            "generator_id",
            "utility_id_eia",
            "prime_mover_code",
            "capacity_mw",
            "nameplate_power_factor",
        ],
    ]

    # default steam turbines to flash steam b/c flash is more common than steam according to EIA
    # Source: https://www.eia.gov/energyexplained/geothermal/geothermal-power-plants.php
    geo_map = {"BT": "B", "ST": "F"}
    geothermal_geotype["geotype_code"] = geothermal_geotype["prime_mover_code"].map(
        geo_map
    )

    # According to the NREL Geothermal report (https://www.nrel.gov/docs/fy21osti/78291.pdf)
    # The Geysers Complex in California contains the only "Dry Steam" (geotype S) plants in the country
    # The Geysers has utility_id_eia = 7160
    geothermal_geotype.loc[
        (geothermal_geotype["utility_id_eia"] == 7160)
        & (geothermal_geotype["geotype_code"] == "F"),
        "geotype_code",
    ] = "S"

    # calculate what fraction of the plant's nameplate capcity each generator is responsible for
    # assume that missing power factor is 100%
    geothermal_geotype["nameplate_power_factor"] = geothermal_geotype[
        "nameplate_power_factor"
    ].fillna(1)
    # calculate adjusted nameplate capacity based on power factor
    geothermal_geotype["pf_adjusted_capacity"] = (
        geothermal_geotype["capacity_mw"] * geothermal_geotype["nameplate_power_factor"]
    )
    # sum adjusted capacity by plant and merge back into generator-level data
    total_plant_capacity = (
        geothermal_geotype.groupby("plant_id_eia", dropna=False)["pf_adjusted_capacity"]
        .sum()
        .reset_index()
        .rename(columns={"pf_adjusted_capacity": "plant_capacity_total"})
    )
    geothermal_geotype = geothermal_geotype.merge(
        total_plant_capacity, how="left", on="plant_id_eia"
    )
    # calculate the fraction
    geothermal_geotype["plant_frac"] = (
        geothermal_geotype["pf_adjusted_capacity"]
        / geothermal_geotype["plant_capacity_total"]
    )
    # drop intermediate columns
    geothermal_geotype = geothermal_geotype.drop(
        columns=[
            "prime_mover_code",
            "utility_id_eia",
            "capacity_mw",
            "nameplate_power_factor",
            "pf_adjusted_capacity",
            "plant_capacity_total",
        ]
    )

    return geothermal_geotype


def adjust_emissions_for_CHP(df):
    """Allocates total emissions for electricity generation."""

    # calculate the electric allocation factor
    df = calculate_electric_allocation_factor(df)

    if "co2_mass_lb" in df.columns:
        df["co2_mass_lb_for_electricity"] = (
            df["co2_mass_lb"] * df["electric_allocation_factor"]
        )
    if "ch4_mass_lb" in df.columns:
        df["ch4_mass_lb_for_electricity"] = (
            df["ch4_mass_lb"] * df["electric_allocation_factor"]
        )
    if "n2o_mass_lb" in df.columns:
        df["n2o_mass_lb_for_electricity"] = (
            df["n2o_mass_lb"] * df["electric_allocation_factor"]
        )
    if "nox_mass_lb" in df.columns:
        df["nox_mass_lb_for_electricity"] = (
            df["nox_mass_lb"] * df["electric_allocation_factor"]
        )
    if "so2_mass_lb" in df.columns:
        df["so2_mass_lb_for_electricity"] = (
            df["so2_mass_lb"] * df["electric_allocation_factor"]
        )

    df = df.drop(columns=["electric_allocation_factor"])

    return df


def calculate_electric_allocation_factor(df):

    mwh_to_mmbtu = 3.412142

    # calculate the useful thermal output
    # 0.8 is an assumed efficiency factor used by eGRID
    df["useful_thermal_output"] = 0.8 * (
        df["fuel_consumed_mmbtu"] - df["fuel_consumed_for_electricity_mmbtu"]
    )

    # convert generation to mmbtu
    try:
        df["generation_mmbtu"] = df["net_generation_mwh"] * mwh_to_mmbtu
    # for CEMS use gross generation
    # TODO: investigate if this works correctly
    except KeyError:
        df["generation_mmbtu"] = df["gross_generation_mwh"] * mwh_to_mmbtu

    # calculate the electric allocation factor
    # 0.75 is an assumed efficiency factor used by eGRID
    df["electric_allocation_factor"] = df["generation_mmbtu"] / (
        df["generation_mmbtu"] + (0.75 * df["useful_thermal_output"])
    )

    # if the allocation factor < 0, set to zero
    df.loc[df["electric_allocation_factor"] < 0, "electric_allocation_factor"] = 0
    # if the allocation factor > 1, set to one
    df.loc[df["electric_allocation_factor"] > 1, "electric_allocation_factor"] = 1
    # fill any missing factors with 1
    df["electric_allocation_factor"] = df["electric_allocation_factor"].fillna(1)

    # remove intermediate columns
    df = df.drop(columns=["useful_thermal_output", "generation_mmbtu"])

    return df


def adjust_emissions_for_biomass(df):
    """Creates a new adjusted co2 emissions column that sets any biomass emissions to zero."""

    # create a column for adjusted biomass emissions, setting these emissions to zero
    biomass_fuels = [
        "AB",
        "BG",
        "BLQ",
        "DG",
        "LFG",
        "MSB",
        "OBG",
        "OBL",
        "OBS",
        "SLW",
        "WDL",
        "WDS",
    ]

    # adjust emissions for co2 for all biomass generators
    if "co2_mass_lb_for_electricity" in df.columns:
        df["co2_mass_lb_adjusted"] = df["co2_mass_lb_for_electricity"]
        df.loc[df["energy_source_code"].isin(biomass_fuels), "co2_mass_lb_adjusted"] = 0
    # for landfill gas (LFG), all other emissions are set to zero
    # this assumes that the gas would have been flared anyway if not used for electricity generation
    if "ch4_mass_lb_for_electricity" in df.columns:
        df["ch4_mass_lb_adjusted"] = df["ch4_mass_lb_for_electricity"]
        df.loc[df["energy_source_code"] == "LFG", "ch4_mass_lb_adjusted"] = 0
    if "n2o_mass_lb_for_electricity" in df.columns:
        df["n2o_mass_lb_adjusted"] = df["n2o_mass_lb_for_electricity"]
        df.loc[df["energy_source_code"] == "LFG", "n2o_mass_lb_adjusted"] = 0
    # nox gets assigned an adjusted value
    # this value is based on using NOx emissions from flaring as a baseline, and subtracting this from the actual emissions
    # to prevent negative emissions, we set the value = 0 if negative
    if "nox_mass_lb_for_electricity" in df.columns:
        df["nox_mass_lb_adjusted"] = df["nox_mass_lb_for_electricity"]
        df.loc[df["energy_source_code"] == "LFG", "nox_mass_lb_adjusted"] = df.loc[
            df["energy_source_code"] == "LFG", "nox_mass_lb_adjusted"
        ] - (
            df.loc[
                df["energy_source_code"] == "LFG", "fuel_consumed_for_electricity_mmbtu"
            ]
            * 0.078
        )
        df.loc[df["nox_mass_lb_adjusted"] < 0, "nox_mass_lb_adjusted"] = 0
    if "so2_mass_lb_for_electricity" in df.columns:
        df["so2_mass_lb_adjusted"] = df["so2_mass_lb_for_electricity"]
        df.loc[df["energy_source_code"] == "LFG", "so2_mass_lb_adjusted"] = 0

    return df


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
            f"   Removing {len(plants_in_states_to_remove)} plants located in the following states: {remove_states}"
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
            "../data/manual/egrid_static_tables/table_4-2_plants_not_connected_to_grid.csv",
            dtype=get_dtypes(),
        )["Plant ID"]
    )

    num_plants = len(
        df[df["plant_id_eia"].isin(ngc_plants)]["plant_id_eia"].unique()
    ) + len(
        df[(df["plant_id_eia"] >= 880000) & (df["plant_id_eia"] < 890000)][
            "plant_id_eia"
        ].unique()
    )
    print(f"   Removing {num_plants} plants that are not grid-connected")

    df = df[~df["plant_id_eia"].isin(ngc_plants)]

    # according to the egrid documentation, any plants that have an id of 88XXXX are not grid connected
    # only keep plants that dont have an id of 88XXXX
    df = df[(df["plant_id_eia"] < 880000) | (df["plant_id_eia"] >= 890000)]

    return df


def clean_cems(year, small):
    """
    Coordinating function for all of the cems data cleaning
    """
    # load the CEMS data
    cems = load_data.load_cems_data(year)

    if small:
        cems = smallerize_test_data(df=cems, random_seed=42)

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

    # add a report date
    cems = add_report_date(cems)

    # TODO: identify and remove any hourly values that appear to be outliers

    # add a fuel type to each observation
    cems = assign_fuel_type_to_cems(cems, year)

    # fill in missing hourly emissions data using the fuel type and heat input
    cems = fill_cems_missing_co2(cems, year)

    # TODO: Add functions for filling missing NOx and SOx

    # calculate ch4 and n2o emissions
    cems = calculate_ghg_emissions_from_fuel_consumption(
        df=cems, year=year, include_co2=False, include_ch4=True, include_n2o=True
    )

    # remove any observations from cems where zero operation is reported for an entire month
    # although this data could be considered to be accurately reported, let's remove it so that we can double check against the eia data
    # TODO: check if any of these observations are from geothermal generators
    cems = remove_cems_with_zero_monthly_data(cems)

    # calculated CHP-adjusted emissions
    cems = calculate_electric_fuel_consumption_for_cems(cems)
    cems = adjust_emissions_for_CHP(cems)

    # calculate biomass-adjusted emissions
    cems = adjust_emissions_for_biomass(cems)

    # add subplant id
    subplant_crosswalk = pd.read_csv(
        f"../data/outputs/{year}/subplant_crosswalk.csv", dtype=get_dtypes()
    )[["plant_id_eia", "unitid", "subplant_id"]].drop_duplicates()
    cems = cems.merge(subplant_crosswalk, how="left", on=["plant_id_eia", "unitid"])

    cems = apply_dtypes(cems)

    return cems


def smallerize_test_data(df, random_seed=None):
    print("   Randomly selecting 5% of plants for faster test run.")
    # Select 5% of plants
    selected_plants = df.plant_id_eia.unique()
    if random_seed is not None:
        np.random.seed(random_seed)
    selected_plants = np.random.choice(
        selected_plants, size=int(len(selected_plants) * 0.05), replace=False
    )
    # Filter for selected plants
    df = df[df.plant_id_eia.isin(selected_plants)]

    return df


def manually_remove_steam_units(df):
    """
    Removes any records from CEMS that we've identified as being steam only plants that need to be removed
    """

    # get the list of plant_id_eia from the static table
    units_to_remove = pd.read_csv(
        "../data/manual/steam_units_to_remove.csv", dtype=get_dtypes(),
    )[["plant_id_eia", "unitid"]]

    print(
        f"   Removing {len(units_to_remove)} units that only produce steam and do not report to EIA"
    )

    df = df.merge(
        units_to_remove, how="outer", on=["plant_id_eia", "unitid"], indicator="source"
    )
    df = df[df["source"] == "left_only"].drop(columns=["source"])

    return df


def add_report_date(df):
    """
    Add a report date column to the cems data based on the plant's local timezone

    Args:
        df (pd.Dataframe): dataframe containing 'plant_id_eia' and 'datetime_utc' columns
    Returns:
        Original dataframe with 'report_date' column added
    """
    plants_entity_eia = load_data.load_pudl_table("plants_entity_eia")

    # get timezone
    df = df.merge(
        plants_entity_eia[["plant_id_eia", "timezone"]], how="left", on="plant_id_eia"
    )

    # create a datetimeindex from the datetime_utc column
    datetime_utc = pd.DatetimeIndex(df["datetime_utc"])

    # create blank column to hold local datetimes
    df["report_date"] = np.NaN

    # get list of unique timezones
    timezones = list(df["timezone"].unique())

    # convert UTC to the local timezone
    for tz in timezones:
        tz_mask = df["timezone"] == tz  # find all rows where the tz matches
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Converting to PeriodArray/Index representation will drop timezone information.",
            )
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


def assign_fuel_type_to_cems(cems, year):
    "Assigns a fuel type to each observation in CEMS"

    fuel_types = get_epa_unit_fuel_types(year)

    # merge in the reported fuel type
    cems = cems.merge(fuel_types, how="left", on=["plant_id_epa", "unitid"])

    # fill missing fuel codes for plants that only have a single fuel type
    single_fuel_plants = (
        fuel_types.drop_duplicates(subset=["plant_id_epa", "energy_source_code"])
        .drop_duplicates(subset=["plant_id_epa"], keep=False)
        .drop(columns="unitid")
    )
    cems = cems.merge(
        single_fuel_plants, how="left", on=["plant_id_epa"], suffixes=(None, "_plant")
    )
    cems["energy_source_code"] = cems["energy_source_code"].fillna(
        cems["energy_source_code_plant"]
    )
    cems = cems.drop(columns="energy_source_code_plant")

    # TODO: fill fuel codes for plants that only have a single fossil type identified in EIA
    cems = fill_missing_fuel_for_single_fuel_plant_months(cems, year)

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
        gf.groupby(["plant_id_eia", "report_date", "energy_source_code"], dropna=False)
        .sum()
        .reset_index()
    )

    # only keep rows with >0 fuel consumed reported
    gf = gf[gf["fuel_consumed_mmbtu"] > 0]

    # identify which plant months have multiple fuels reported
    multi_fuels = (
        gf.groupby(["plant_id_eia", "report_date"], dropna=False)["energy_source_code"]
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


def calculate_co2_eq_mass(
    df, ipcc_version="AR5", gwp_horizon=100, ar5_climate_carbon_feedback=False
):
    """
    Calculate CO2-equivalent emissions from CO2, CH4, and N2O. This is done
    by choosing one of the IPCC's emission factors for CH4 and N2O.

    Inputs:
        df: Should contain at least: ['co2_mass_lb', 'ch4_mass_lb', 'n2o_mass_lb']
    
    If the `fuel_consumed_for_electricity_units` column is available, we also
    compute the adjusted emissions.
    """
    df_gwp = load_data.load_ipcc_gwp()

    if ipcc_version not in ("SAR", "TAR", "AR4", "AR5"):
        raise ValueError("Unsupported option for `ipcc_version`.")
    if gwp_horizon not in (20, 100):
        raise ValueError(
            "Only 20-year and 100-year global warming potentials are supported."
        )
    if ar5_climate_carbon_feedback and ipcc_version not in ("AR5"):
        raise ValueError("Climate carbon feedback (CCF) is only available for AR5.")

    if ar5_climate_carbon_feedback:
        ipcc_version += "f"

    ch4_gwp_factor = df_gwp.loc[ipcc_version][f"ch4_{gwp_horizon}_year"].astype(float)
    n2o_gwp_factor = df_gwp.loc[ipcc_version][f"n2o_{gwp_horizon}_year"].astype(float)

    if (
        "co2_mass_lb" not in df.columns
        or "ch4_mass_lb" not in df.columns
        or "n2o_mass_lb" not in df.columns
    ):
        raise ValueError(
            "Make sure the input dataframe has emissions data for CO2, CH4, and N2O."
        )

    df["co2_eq_mass_lb"] = (
        df["co2_mass_lb"]
        + ch4_gwp_factor * df["ch4_mass_lb"]
        + n2o_gwp_factor * df["n2o_mass_lb"]
    )

    if "co2_mass_lb_adjusted" in df:
        df["co2_eq_mass_lb_adjusted"] = (
            df["co2_mass_lb_adjusted"]
            + ch4_gwp_factor * df["ch4_mass_lb_adjusted"]
            + n2o_gwp_factor * df["n2o_mass_lb_adjusted"]
        )
    if "co2_mass_lb_for_electricity" in df:
        df["co2_eq_mass_lb_for_electricity"] = (
            df["co2_mass_lb_for_electricity"]
            + ch4_gwp_factor * df["ch4_mass_lb_for_electricity"]
            + n2o_gwp_factor * df["n2o_mass_lb_for_electricity"]
        )

    return df


def calculate_nox_from_fuel_consumption(
    df: pd.DataFrame, pudl_out, year
) -> pd.DataFrame:
    """
    Calculate NOx emissions from fuel consumption data.

    Inputs:
        df: Should contain the following columns:
            [`plant_id_eia`, `report_date`, `fuel_consumed_units`, `energy_source_code`, `prime_mover_code`]
    
    If the `fuel_consumed_for_electricity_units` column is available, we also
    compute the adjusted emissions.
    """
    emission_factors = load_data.load_nox_emission_factors()
    # remove emissions factors where the unit is mmbtu
    emission_factors = emission_factors[
        emission_factors["emission_factor_denominator"] != "mmbtu"
    ]
    # for now, we do not have information about the boiler firing type
    # thus, we will average the factors by fuel and prime mover
    emission_factors = (
        emission_factors.groupby(
            ["energy_source_code", "prime_mover_code"], dropna=False
        )
        .mean()
        .reset_index()
    )
    # merge in the emission factor
    df = df.merge(
        emission_factors,
        how="left",
        on=["energy_source_code", "prime_mover_code"],
        validate="m:1",
    )
    # fill missing factors with zero
    df["emission_factor"] = df["emission_factor"].fillna(0)

    # load information about the monthly heat input of fuels
    plant_heat_content = pudl_out.gf_eia923().loc[
        :,
        [
            "plant_id_eia",
            "energy_source_code",
            "prime_mover_code",
            "report_date",
            "fuel_mmbtu_per_unit",
        ],
    ]
    # replace zero heat content with missing values
    plant_heat_content["fuel_mmbtu_per_unit"] = plant_heat_content[
        "fuel_mmbtu_per_unit"
    ].replace(0, np.NaN)
    # calculate the average monthly heat content for a fuel
    fuel_heat_content = (
        plant_heat_content.drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "report_date"], dropna=False)
        .mean()
        .reset_index()
    )

    # change the report date columns back to datetimes
    plant_heat_content["report_date"] = pd.to_datetime(
        plant_heat_content["report_date"]
    )
    fuel_heat_content["report_date"] = pd.to_datetime(fuel_heat_content["report_date"])

    # merge the heat content, starting with plant-specific values, then filling using fuel-specific values
    df = df.merge(
        plant_heat_content,
        how="left",
        on=["plant_id_eia", "energy_source_code", "prime_mover_code", "report_date"],
        validate="m:1",
    )
    df = df.merge(
        fuel_heat_content,
        how="left",
        on=["energy_source_code", "report_date"],
        validate="m:1",
        suffixes=(None, "_generic"),
    )
    df["fuel_mmbtu_per_unit"] = df["fuel_mmbtu_per_unit"].fillna(
        df["fuel_mmbtu_per_unit_generic"]
    )

    # calculate the nox emissions mass
    df["nox_mass_lb"] = (df["fuel_consumed_mmbtu"] / df["fuel_mmbtu_per_unit"]) * df[
        "emission_factor"
    ]

    if df["energy_source_code"].str.contains("GEO").any():
        df = add_geothermal_emission_factors(
            df, year, include_co2=False, include_nox=True, include_so2=False
        )
        df.loc[df["energy_source_code"] == "GEO", "nox_mass_lb"] = (
            df.loc[df["energy_source_code"] == "GEO", "fuel_consumed_mmbtu"]
            * df.loc[df["energy_source_code"] == "GEO", "nox_lb_per_mmbtu"]
        )
        df = df.drop(columns=["nox_lb_per_mmbtu"])

    # Drop intermediate columns.
    df = df.drop(
        columns=[
            "fuel_mmbtu_per_unit",
            "fuel_mmbtu_per_unit_generic",
            "emission_factor",
        ]
    )

    return df


def calculate_so2_from_fuel_consumption(
    df: pd.DataFrame, pudl_out, year
) -> pd.DataFrame:
    """
    Calculate SO2 emissions from fuel consumption data and fuel sulfur content.

    Inputs:
        df: Should contain the following columns:
            [`plant_id_eia`, `report_date`, `fuel_consumed_units`, `energy_source_code`, `prime_mover_code`]
    
    If the `fuel_consumed_for_electricity_units` column is available, we also
    compute the adjusted emissions.
    """
    # load the emission factors
    emission_factors = load_data.load_so2_emission_factors()
    # for now, we do not have information about the boiler firing type
    # thus, we will average the factors by fuel and prime mover
    emission_factors = (
        emission_factors.groupby(
            [
                "energy_source_code",
                "prime_mover_code",
                "emission_factor_denominator",
                "multiply_by_sulfur_content",
            ],
            dropna=False,
        )
        .mean()
        .reset_index()
    )
    # drop all factors for OTH fuel type, since unit is unknown
    emission_factors = emission_factors[emission_factors["energy_source_code"] != "OTH"]
    # move the mmbtu emision factors to a separate df
    emission_factors_mmbtu = emission_factors[
        emission_factors["emission_factor_denominator"] == "mmbtu"
    ]
    emission_factors = emission_factors[
        emission_factors["emission_factor_denominator"] != "mmbtu"
    ]
    # merge in the emission factor
    df = df.merge(
        emission_factors,
        how="left",
        on=["energy_source_code", "prime_mover_code"],
        validate="m:1",
    )
    # merge in the mmbtu emission factors
    df = df.merge(
        emission_factors_mmbtu,
        how="left",
        on=["energy_source_code", "prime_mover_code"],
        validate="m:1",
        suffixes=(None, "_mmbtu"),
    )
    # fill missing factors with mmbtu factors, if available
    df["emission_factor"] = df["emission_factor"].fillna(df["emission_factor_mmbtu"])
    df["multiply_by_sulfur_content"] = df["multiply_by_sulfur_content"].fillna(
        df["multiply_by_sulfur_content_mmbtu"]
    )
    df["emission_factor_denominator"] = df["emission_factor_denominator"].fillna(
        df["emission_factor_denominator_mmbtu"]
    )
    # fill missing factors with zero
    df["emission_factor"] = df["emission_factor"].fillna(0)

    # load the sulfur content
    plant_sulfur_content = pudl_out.bf_eia923().loc[
        :,
        [
            "plant_id_eia",
            "boiler_id",
            "energy_source_code",
            "report_date",
            "sulfur_content_pct",
        ],
    ]
    # merge in the prime mover data
    plant_sulfur_content = plant_sulfur_content.merge(
        pd.read_sql("boilers_entity_eia", pudl_out.pudl_engine),
        how="left",
        on=["plant_id_eia", "boiler_id"],
    ).drop(columns=["boiler_id"])

    # replace zero heat content with missing values
    plant_sulfur_content["sulfur_content_pct"] = plant_sulfur_content[
        "sulfur_content_pct"
    ].replace(0, np.NaN)
    # average the values by plant/PM/ESC
    plant_sulfur_content = (
        plant_sulfur_content.groupby(
            ["plant_id_eia", "energy_source_code", "prime_mover_code", "report_date"],
            dropna=False,
        )
        .mean()
        .reset_index()
    )
    # calculate the average monthly sulfur content for a fuel
    fuel_sulfur_content = (
        plant_sulfur_content.drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "report_date"], dropna=False)
        .mean()
        .reset_index()
    )
    # change the report date columns back to datetimes
    plant_sulfur_content["report_date"] = pd.to_datetime(
        plant_sulfur_content["report_date"]
    )
    fuel_sulfur_content["report_date"] = pd.to_datetime(
        fuel_sulfur_content["report_date"]
    )
    # merge the heat content, starting with plant-specific values, then filling using fuel-specific values
    df = df.merge(
        plant_sulfur_content,
        how="left",
        on=["plant_id_eia", "energy_source_code", "prime_mover_code", "report_date"],
        validate="m:1",
    )
    df = df.merge(
        fuel_sulfur_content,
        how="left",
        on=["energy_source_code", "report_date"],
        validate="m:1",
        suffixes=(None, "_generic"),
    )
    df["sulfur_content_pct"] = df["sulfur_content_pct"].fillna(
        df["sulfur_content_pct_generic"]
    )

    # load information about the monthly heat input of fuels
    plant_heat_content = pudl_out.gf_eia923().loc[
        :,
        [
            "plant_id_eia",
            "energy_source_code",
            "prime_mover_code",
            "report_date",
            "fuel_mmbtu_per_unit",
        ],
    ]
    # replace zero heat content with missing values
    plant_heat_content["fuel_mmbtu_per_unit"] = plant_heat_content[
        "fuel_mmbtu_per_unit"
    ].replace(0, np.NaN)
    # calculate the average monthly heat content for a fuel
    fuel_heat_content = (
        plant_heat_content.drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "report_date"], dropna=False)
        .mean()
        .reset_index()
    )
    # change the report date columns back to datetimes
    plant_heat_content["report_date"] = pd.to_datetime(
        plant_heat_content["report_date"]
    )
    fuel_heat_content["report_date"] = pd.to_datetime(fuel_heat_content["report_date"])
    # merge the heat content, starting with plant-specific values, then filling using fuel-specific values
    df = df.merge(
        plant_heat_content,
        how="left",
        on=["plant_id_eia", "energy_source_code", "prime_mover_code", "report_date"],
        validate="m:1",
    )
    df = df.merge(
        fuel_heat_content,
        how="left",
        on=["energy_source_code", "report_date"],
        validate="m:1",
        suffixes=(None, "_generic"),
    )
    df["fuel_mmbtu_per_unit"] = df["fuel_mmbtu_per_unit"].fillna(
        df["fuel_mmbtu_per_unit_generic"]
    )

    # update the emission factor for those generators where it needs to be multiplied by sulfur content
    df.loc[df["multiply_by_sulfur_content"] == 1, "emission_factor"] = (
        df.loc[df["multiply_by_sulfur_content"] == 1, "emission_factor"]
        * df.loc[df["multiply_by_sulfur_content"] == 1, "sulfur_content_pct"]
    )

    # multiply physical fuel consumption by the emission factor
    df["so2_mass_lb"] = (df["fuel_consumed_mmbtu"] / df["fuel_mmbtu_per_unit"]) * df[
        "emission_factor"
    ]

    # where the emission factor denominator is mmbtu, multiply by fuel consumed instead of physical units
    df.loc[df["emission_factor_denominator"] == "mmbtu", "so2_mass_lb"] = (
        df.loc[df["emission_factor_denominator"] == "mmbtu", "fuel_consumed_mmbtu"]
        * df.loc[df["emission_factor_denominator"] == "mmbtu", "emission_factor"]
    )

    if df["energy_source_code"].str.contains("GEO").any():
        df = add_geothermal_emission_factors(
            df, year, include_co2=False, include_nox=False, include_so2=True
        )
        df.loc[df["energy_source_code"] == "GEO", "so2_mass_lb"] = (
            df.loc[df["energy_source_code"] == "GEO", "fuel_consumed_mmbtu"]
            * df.loc[df["energy_source_code"] == "GEO", "so2_lb_per_mmbtu"]
        )
        df = df.drop(columns=["so2_lb_per_mmbtu"])

    # Drop intermediate columns.
    df = df.drop(
        columns=[
            "emission_factor",
            "emission_factor_mmbtu",
            "multiply_by_sulfur_content",
            "multiply_by_sulfur_content_mmbtu",
            "emission_factor_denominator",
            "emission_factor_denominator_mmbtu",
            "fuel_mmbtu_per_unit",
            "fuel_mmbtu_per_unit_generic",
            "emission_factor",
            "sulfur_content_pct",
            "sulfur_content_pct_generic",
        ]
    )

    return df


def fill_cems_missing_co2(cems, year):
    """
    Fills missing hourly CO2 data in CEMS based on a two-tiered approach.

    CO2 data is considered missing if reported CO2 is zero and fuel consumption is positive.
    If a unit has a unit-specific fuel type identified by the EPA-EIA crosswalk, calculate co2 using a fuel-specific emission factor.
    If not, fill missing data by calculating a plant-month weighted average emission factor of all fuels burned in that plant-month.
    """

    # add a new categorical option to the mass measurement code
    cems["co2_mass_measurement_code"] = cems[
        "co2_mass_measurement_code"
    ].cat.add_categories("Imputed")

    # replace all "missing" CO2 values with zero
    cems["co2_mass_lb"] = cems["co2_mass_lb"].fillna(0)

    # replace 0 reported CO2 values with missing values, if there was reported heat input
    cems.loc[
        (cems["co2_mass_lb"] == 0) & (cems["fuel_consumed_mmbtu"] > 0), "co2_mass_lb",
    ] = np.NaN

    # create a new df with all observations with missing co2 data
    missing_co2 = cems[cems["co2_mass_lb"].isnull()]

    # First round of filling using fuel types in PSDC

    # for rows that have a successful fuel code match, move to a temporary dataframe to hold the data
    co2_to_fill = missing_co2.copy()[~missing_co2["energy_source_code"].isna()]
    fill_index = co2_to_fill.index

    # remove these from the missing co2 dataframe. We'll need to apply a different method for these remaining plants
    missing_co2 = missing_co2[missing_co2["energy_source_code"].isna()]
    missing_index = missing_co2.index

    # calculate emissions based on fuel type
    co2_to_fill = calculate_ghg_emissions_from_fuel_consumption(
        df=co2_to_fill,
        year=year,
        include_co2=True,
        include_ch4=False,
        include_n2o=False,
    ).set_index(fill_index)

    # fill this data into the original cems data
    cems.update(co2_to_fill[["co2_mass_lb"]])

    # update the co2 mass measurement code
    cems.loc[co2_to_fill.index, "co2_mass_measurement_code"] = "Imputed"

    # Second round of data filling using weighted average EF based on EIA-923 heat input data

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
        ["plant_id_eia", "report_date", "energy_source_code"], dropna=False
    ).sum()[["fuel_consumed_for_electricity_mmbtu"]]

    # calculate the percent of heat input from each fuel in each month
    missing_gf = (
        missing_gf
        / missing_gf.reset_index()
        .groupby(["plant_id_eia", "report_date"], dropna=False)
        .sum()
    )

    missing_gf = missing_gf.fillna(1)

    emission_factors = load_data.load_ghg_emission_factors()[
        ["energy_source_code", "co2_lb_per_mmbtu"]
    ]

    # merge in the emission factor
    missing_gf = missing_gf.reset_index().merge(
        emission_factors, how="left", on="energy_source_code"
    )

    # calculate weighted emission factor
    missing_gf["weighted_ef"] = (
        missing_gf["fuel_consumed_for_electricity_mmbtu"]
        * missing_gf["co2_lb_per_mmbtu"]
    )
    missing_gf = (
        missing_gf.groupby(["plant_id_eia", "report_date"], dropna=False)
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
    missing_co2["co2_mass_lb"] = (
        missing_co2["fuel_consumed_mmbtu"] * missing_co2["weighted_ef"]
    )

    # update in CEMS table
    cems.update(missing_co2[["co2_mass_lb"]])

    # update the co2 mass measurement code
    cems.loc[missing_co2.index, "co2_mass_measurement_code"] = "Imputed"

    return cems


def remove_cems_with_zero_monthly_data(cems):
    """
    Identifies months where zero generation, heat input, or emissions are reported 
    from each unit and removes associated hours from CEMS so that these can be filled using the eia923 data
    Inputs:
        cems: pandas dataframe of hourly cems data containing columns "plant_id_eia", "unitid" and "report_date"
    Returns:
        cems df with hourly observations for months when no emissions reported removed
    """
    # calculate the totals reported in each month
    cems_with_zero_monthly_emissions = cems.groupby(
        ["plant_id_eia", "unitid", "report_date"], dropna=False
    ).sum()[
        [
            "co2_mass_lb",
            "nox_mass_lb",
            "so2_mass_lb",
            "gross_generation_mwh",
            "fuel_consumed_mmbtu",
        ]
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
            ["plant_id_eia", "unitid", "report_date", "missing_data_flag"]
        ],
        how="left",
        on=["plant_id_eia", "unitid", "report_date"],
    )
    # remove any observations with the missing data flag
    print(
        f"   Removing {len(cems[cems['missing_data_flag'] == 'remove'])} observations from cems for unit-months where no data reported"
    )
    cems = cems[cems["missing_data_flag"] != "remove"]
    # drop the missing data flag column
    cems = cems.drop(columns="missing_data_flag")

    return cems


def calculate_electric_fuel_consumption_for_cems(cems, drop_interim_columns=True):
    """
    Calculates the portion of fuel consumption and CO2 emissions for electricity for each hour in CEMS.
    """
    # factors to convert to MMBTU
    mwh_to_mmbtu = 3.412142
    klb_to_mmbtu = 1.194  # NOTE: this might differ for each plant

    # calculate total heat output
    cems["heat_output_mmbtu"] = (cems["gross_generation_mwh"] * mwh_to_mmbtu) + (
        cems["steam_load_1000_lb"] * klb_to_mmbtu
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

    if drop_interim_columns:
        cems = cems.drop(columns=["heat_output_mmbtu", "frac_electricity"])

    return cems


def identify_hourly_data_source(eia923_allocated, cems, year):
    """Identifies whether there is hourly CEMS data available for each subplant-month.
    Possible categories:
        1. `cems`: For subplant-months for which we have hourly CEMS data for all CEMS units that make up that subplant,
            we will use the hourly values reported in CEMS. (Add a validation check for the net generation and fuel consumption totals)
        2. `partial_cems`: For subplant-months for which we have hourly CEMS data 
            for only some of the CEMS units that make up a subplant, we will use the reported 
            EIA-923 values to scale the partial hourly CEMS data from the other units to match the total value for the entire subplant. This will also calculate a partial subplant scaling factor for each data column (e.g. net generation, fuel consumption) by comparing the total monthly CEMS data to the monthly EIA-923 data.
        3. `eia`: for subplant-months for which no hourly data is reported in CEMS, 
            we will attempt to use EIA-930 data to assign an hourly profile to the monthly EIA-923 data
    Inputs:
        eia923_allocated:
        cems:
        year:
    Returns:
        eia923_allocated with new column `hourly_data_source`
    """

    # aggregate cems data to plant-unit-month
    cems_monthly = (
        cems.groupby(
            ["plant_id_eia", "subplant_id", "unitid", "report_date"], dropna=False
        )
        .sum()[
            [
                "net_generation_mwh",
                "fuel_consumed_mmbtu",
                "fuel_consumed_for_electricity_mmbtu",
                "co2_mass_lb",
            ]
        ]
        .reset_index()
    )

    all_data = eia923_allocated.copy()

    # create a binary column indicating whether data was reported in 923
    columns_to_test = [
        "net_generation_mwh",
        "fuel_consumed_mmbtu",
        "fuel_consumed_for_electricity_mmbtu",
        "co2_mass_lb",
        "co2_mass_lb_adjusted",
    ]
    all_data = all_data.assign(
        reported_eia923=lambda x: np.where(
            x[columns_to_test].notnull().all(axis=1), 1, 0
        )
    )

    # load the subplant crosswalk and identify unique unitids in each subplant
    units_in_subplant = pd.read_csv(
        f"../data/outputs/{year}/subplant_crosswalk.csv",
        dtype=get_dtypes(),
        parse_dates=["current_planned_operating_date", "retirement_date"],
    )[["plant_id_eia", "unitid", "subplant_id", "retirement_date"]].drop_duplicates()

    # remove units that retired before the current year
    units_in_subplant = units_in_subplant[
        ~(units_in_subplant["retirement_date"].dt.year < year)
    ]

    # get a count of the number of CEMS units in each subplant
    units_in_subplant = (
        units_in_subplant.groupby(["plant_id_eia", "subplant_id"], dropna=False)
        .count()["unitid"]
        .reset_index()
        .rename(columns={"unitid": "units_in_subplant"})
    )

    # create a dataframe that counts the number of units reported in CEMS in each subplant-month
    cems_units_reported = (
        cems_monthly.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .count()["unitid"]
        .reset_index()
        .rename(columns={"unitid": "subplant_units_reported"})
    )

    # merge in the total number of units that exist in each subplant
    # this will allow us to compare where a subplant-month is missing data from one or more units
    cems_units_reported = cems_units_reported.merge(
        units_in_subplant, how="left", on=["plant_id_eia", "subplant_id"]
    )

    # identify which subplant-months have complete or partial cems data
    cems_units_reported = cems_units_reported.assign(
        hourly_data_source=lambda x: np.where(
            (x.subplant_units_reported < x.units_in_subplant), "partial_cems", "cems"
        )
    )

    # merge in the data source column from CEMS
    all_data = all_data.merge(
        cems_units_reported[
            ["plant_id_eia", "subplant_id", "report_date", "hourly_data_source"]
        ],
        how="left",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="m:1",
    )

    # for the remaining plants, identify the hourly data source as EIA
    all_data["hourly_data_source"] = all_data["hourly_data_source"].fillna("eia")

    # remove any generator-months for which there is no data reported in either data source
    all_data = all_data[
        ~(
            (all_data["reported_eia923"] == 0)
            & (all_data["hourly_data_source"] == "eia")
        )
    ]

    all_data = all_data.drop(columns=["reported_eia923"])

    return all_data


def convert_gross_to_net_generation(cems, eia923_allocated, plant_attributes):
    """
    Converts hourly gross generation in CEMS to hourly net generation by calculating a gross to net generation ratio
    Inputs:

    Returns: 
        cems df with an added column for net_generation_mwh and a column indicated the method used to calculate net generation
    """

    gtn_conversions = calculate_gross_to_net_conversion_factors(
        cems, eia923_allocated, plant_attributes
    )

    factors_to_use = gtn_conversions[
        [
            "plant_id_eia",
            "subplant_id",
            "report_date",
            "hourly_shift_mw_monthly",
            "hourly_shift_mw_annual",
            "annual_plant_ratio",
            "annual_fuel_ratio",
        ]
    ]

    # merge the conversion factors we want to use into the cems data
    cems = cems.merge(
        factors_to_use, how="left", on=["plant_id_eia", "subplant_id", "report_date"]
    )
    # count the number of units in each subplant
    units_in_subplant = cems[
        ["plant_id_eia", "subplant_id", "report_date", "unitid"]
    ].drop_duplicates()
    units_in_subplant = (
        units_in_subplant.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .count()
        .reset_index()
        .rename(columns={"unitid": "units_in_subplant"})
    )
    cems = cems.merge(
        units_in_subplant, how="left", on=["plant_id_eia", "subplant_id", "report_date"]
    )

    cems["gtn_method"] = "monthly_shift_factor"
    # calculate net generation using the monthly shift factors where available
    cems["net_generation_mwh"] = cems["gross_generation_mwh"] + (
        cems["hourly_shift_mw_monthly"] / cems["units_in_subplant"]
    )
    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "annual_shift_factor"
    # next use the annual shift factor where available
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"]
        + (cems["hourly_shift_mw_annual"] / cems["units_in_subplant"])
    )
    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "annual_plant_ratio"
    # next use the annual plant ratio
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["annual_plant_ratio"]
    )
    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "annual_fuel_ratio"
    # next use the annual fuel ratio
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"] * cems["annual_fuel_ratio"]
    )
    cems.loc[cems["net_generation_mwh"].isna(), "gtn_method"] = "gross_as_net"
    # if nothing else is abailable, use the gross generation value
    cems["net_generation_mwh"] = cems["net_generation_mwh"].fillna(
        cems["gross_generation_mwh"]
    )

    # drop intermediate columns
    cems = cems.drop(
        columns=[
            "hourly_shift_mw_monthly",
            "hourly_shift_mw_annual",
            "annual_plant_ratio",
            "annual_fuel_ratio",
            "units_in_subplant",
        ]
    )

    return cems, gtn_conversions


def calculate_gross_to_net_conversion_factors(cems, eia923_allocated, plant_attributes):
    """
    Calculates gross to net ratios and shift factors
    """
    gross_gen_data = (
        cems.groupby(
            ["plant_id_eia", "subplant_id", "report_date", "datetime_utc"], dropna=False
        )
        .sum()["gross_generation_mwh"]
        .reset_index()
    )
    gross_gen_data = (
        gross_gen_data.groupby(
            ["plant_id_eia", "subplant_id", "report_date"], dropna=False
        )
        .agg({"datetime_utc": "count", "gross_generation_mwh": "sum"})
        .reset_index()
        .rename(columns={"datetime_utc": "hours_in_month"})
    )
    net_gen_data = (
        eia923_allocated.dropna(subset=["net_generation_mwh"])
        .groupby(["plant_id_eia", "subplant_id", "report_date"], dropna=False)
        .sum()["net_generation_mwh"]
        .reset_index()
    )

    # combine monthly gross and net generation data where we have data for both
    gtn_conversions = gross_gen_data.merge(
        net_gen_data,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        indicator="source",
    )

    # drop data that only exists in EIA, but not in CEMS, since there is no gross generation data to calculate NG for
    gtn_conversions = gtn_conversions[~(gtn_conversions["source"] == "right_only")]

    # calculate other groupings at the plant and annual levels
    annual_subplant_ratio = (
        gtn_conversions.dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia", "subplant_id"], dropna=False)
        .sum()[["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]]
        .reset_index()
    )
    monthly_plant_ratio = (
        gtn_conversions.groupby(["plant_id_eia", "report_date"], dropna=False)
        .sum()[["gross_generation_mwh", "net_generation_mwh"]]
        .reset_index()
    )
    annual_plant_ratio = (
        gtn_conversions.dropna(subset=["gross_generation_mwh", "net_generation_mwh"])
        .groupby(["plant_id_eia"], dropna=False)
        .sum()[["gross_generation_mwh", "net_generation_mwh"]]
        .reset_index()
    )

    # calculate the ratios at each aggregation level
    # fill missing values (due to divide by zero) with zero
    # replace infinite values with missing
    gtn_conversions["monthly_subplant_ratio"] = (
        (
            gtn_conversions["net_generation_mwh"]
            / gtn_conversions["gross_generation_mwh"]
        )
        .fillna(0)
        .replace([np.inf, -np.inf], np.nan)
    )
    annual_subplant_ratio["annual_subplant_ratio"] = (
        (
            annual_subplant_ratio["net_generation_mwh"]
            / annual_subplant_ratio["gross_generation_mwh"]
        )
        .fillna(0)
        .replace([np.inf, -np.inf], np.nan)
    )
    monthly_plant_ratio["monthly_plant_ratio"] = (
        (
            monthly_plant_ratio["net_generation_mwh"]
            / monthly_plant_ratio["gross_generation_mwh"]
        )
        .fillna(0)
        .replace([np.inf, -np.inf], np.nan)
    )
    annual_plant_ratio["annual_plant_ratio"] = (
        (
            annual_plant_ratio["net_generation_mwh"]
            / annual_plant_ratio["gross_generation_mwh"]
        )
        .fillna(0)
        .replace([np.inf, -np.inf], np.nan)
    )

    # calculate a monthly and annual shift factor
    gtn_conversions["hourly_shift_mw_monthly"] = (
        gtn_conversions["net_generation_mwh"] - gtn_conversions["gross_generation_mwh"]
    ) / (gtn_conversions["hours_in_month"])
    annual_subplant_ratio["hourly_shift_mw_annual"] = (
        annual_subplant_ratio["net_generation_mwh"]
        - annual_subplant_ratio["gross_generation_mwh"]
    ) / (annual_subplant_ratio["hours_in_month"])

    # drop the gross and net generation data from the dataframes at teh other aggregation levels
    annual_subplant_ratio = annual_subplant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh", "hours_in_month"]
    )
    monthly_plant_ratio = monthly_plant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh"]
    )
    annual_plant_ratio = annual_plant_ratio.drop(
        columns=["gross_generation_mwh", "net_generation_mwh"]
    )

    # merge the various ratios back into a single dataframe
    gtn_conversions = gtn_conversions.merge(
        annual_subplant_ratio, how="left", on=["plant_id_eia", "subplant_id"]
    )
    gtn_conversions = gtn_conversions.merge(
        monthly_plant_ratio, how="left", on=["plant_id_eia", "report_date"]
    )
    gtn_conversions = gtn_conversions.merge(
        annual_plant_ratio, how="left", on=["plant_id_eia"]
    )

    # where gross or net generation data was missing in a month, change the monthly ratios to missing
    gtn_conversions.loc[
        gtn_conversions[["gross_generation_mwh", "net_generation_mwh"]]
        .isna()
        .any(axis=1),
        ["monthly_subplant_ratio", "monthly_plant_ratio"],
    ] = np.NaN

    # calculate the mean ratio for all plants of a single fuel type
    annual_fuel_ratio = (
        annual_plant_ratio.merge(
            plant_attributes[["plant_id_eia", "plant_primary_fuel"]],
            how="left",
            on="plant_id_eia",
        )
        .groupby("plant_primary_fuel")
        .mean()["annual_plant_ratio"]
        .reset_index()
        .rename(columns={"annual_plant_ratio": "annual_fuel_ratio"})
    )

    # merge the plant primary fuel and the fuel ratios into the conversion table
    gtn_conversions = gtn_conversions.merge(
        plant_attributes[["plant_id_eia", "plant_primary_fuel"]],
        how="left",
        on="plant_id_eia",
    )
    gtn_conversions = gtn_conversions.merge(
        annual_fuel_ratio, how="left", on="plant_primary_fuel"
    )

    return gtn_conversions


def convert_gtn_using_ratio(gtn_conversion, level):
    values_to_fill = gtn_conversion[
        gtn_conversion["gtn_ratio"].isna()
        & ~gtn_conversion[f"gtn_ratio_{level}"].isna()
    ].index
    gtn_conversion.loc[values_to_fill, "gtn_ratio"] = gtn_conversion.loc[
        values_to_fill, f"gtn_ratio_{level}"
    ]
    gtn_conversion.loc[values_to_fill, "gtn_constant"] = 0
    gtn_conversion.loc[values_to_fill, "gtn_method"] = f"{level}_ratio"

    return gtn_conversion


def convert_gtn_using_regression(gtn_conversion, level):
    values_to_fill = gtn_conversion[
        gtn_conversion["gtn_ratio"].isna() & ~gtn_conversion[f"slope_{level}"].isna()
    ].index
    gtn_conversion.loc[values_to_fill, "gtn_ratio"] = gtn_conversion.loc[
        values_to_fill, f"slope_{level}"
    ]
    gtn_conversion.loc[values_to_fill, "gtn_constant"] = gtn_conversion.loc[
        values_to_fill, f"intercept_{level}"
    ]
    gtn_conversion.loc[values_to_fill, "gtn_method"] = f"{level}_regression"

    return gtn_conversion


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
        cems.groupby(["plant_id_eia", "report_date"], dropna=False)
        .sum()[
            [
                "gross_generation_mwh",
                "steam_load_1000_lb",
                "fuel_consumed_mmbtu",
                "co2_mass_lb",
            ]
        ]
        .reset_index()
    )
    eia_monthly = (
        gen_fuel_allocated.groupby(["plant_id_eia", "report_date"], dropna=False)
        .sum()
        .reset_index()
    )

    # identify all plant months that report no generation but have heat input and have no steam production
    missing_generation = cems_monthly[
        (cems_monthly["gross_generation_mwh"] == 0)
        & (cems_monthly["fuel_consumed_mmbtu"] > 0)
        & (cems_monthly["steam_load_1000_lb"] == 0)
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
    # missing_generation['steam_heat_rate_mmbtu_per_klb'] = (missing_generation['fuel_consumed_mmbtu'] - missing_generation['fuel_consumed_for_electricity_mmbtu']) / missing_generation['steam_load_1000_lb'] * missing_generation['heat_scaling_factor']
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
    # cems_ng_imputation['heat_content_for_electricity_mmbtu'] = cems_ng_imputation['fuel_consumed_mmbtu'] - (cems_ng_imputation['steam_load_1000_lb'].fillna(0) * cems_ng_imputation['steam_heat_rate_mmbtu_per_klb'])

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


def filter_unique_cems_data(cems, partial_cems):
    """Removes subplant-months from cems that also appear in partial_cems."""
    # filter the cems data to remove data that was scaled in partial_cems
    partial_cems_subplant_months = partial_cems[
        ["plant_id_eia", "subplant_id", "report_date"]
    ].drop_duplicates()
    filtered_cems = cems.merge(
        partial_cems_subplant_months,
        how="outer",
        on=["plant_id_eia", "subplant_id", "report_date"],
        indicator="source",
    )

    filtered_cems = filtered_cems[filtered_cems["source"] == "left_only"].drop(
        columns=["source"]
    )

    return filtered_cems


def aggregate_plant_data_to_ba_fuel(combined_plant_data, plant_frame):
    data_columns = [
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

    ba_fuel_data = combined_plant_data.merge(
        plant_frame, how="left", on=["plant_id_eia"]
    )
    ba_fuel_data = (
        ba_fuel_data.groupby(
            ["ba_code", "fuel_category", "datetime_utc", "report_date"], dropna=False
        )[data_columns]
        .sum()
        .reset_index()
    )
    return ba_fuel_data


def combine_plant_data(cems, partial_cems, shaped_eia_data):
    """
    Combines final hourly subplant data from each source into a single dataframe.
    Inputs:
        Pandas dataframes of shaped or original hourly data

    Note: returns dask dataframe (not used before this point in pipeline) because of data size

    """
    # Convert to dask because we are about to make a GIANT dataframe
    # 2,900,000 rows/partition leads to approx 1GB chunk size
    # cems = dd.from_pandas(cems, npartitions=20)
    # partial_cems = dd.from_pandas(partial_cems, npartitions=20)
    # shaped_eia_data = dd.from_pandas(shaped_eia_data, npartitions=20)

    KEY_COLUMNS = [
        "plant_id_eia",
        "datetime_utc",
        "report_date",
    ]

    DATA_COLUMNS = [
        "gross_generation_mwh",
        "net_generation_mwh",
        "steam_load_1000_lb",
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

    ALL_COLUMNS = KEY_COLUMNS + DATA_COLUMNS

    # group data by plant-hour and filter columns
    cems = (
        cems.groupby(KEY_COLUMNS, dropna=False,)
        .sum()
        .reset_index()[[col for col in cems.columns if col in ALL_COLUMNS]]
    )
    partial_cems = (
        partial_cems.groupby(KEY_COLUMNS, dropna=False,)
        .sum()
        .reset_index()[[col for col in partial_cems.columns if col in ALL_COLUMNS]]
    )
    shaped_eia_data = shaped_eia_data[
        [col for col in shaped_eia_data.columns if col in ALL_COLUMNS]
    ]

    # concat together
    combined_plant_data = pd.concat(
        [cems, partial_cems, shaped_eia_data], axis=0, ignore_index=True, copy=False
    )

    # groupby plant
    combined_plant_data = (
        combined_plant_data.groupby(KEY_COLUMNS, dropna=False).sum().reset_index()
    )

    combined_plant_data[DATA_COLUMNS] = combined_plant_data[DATA_COLUMNS].round(2)

    # re-order the columns
    combined_plant_data = combined_plant_data[ALL_COLUMNS]

    return combined_plant_data


def create_plant_attributes_table(cems, eia923_allocated, year, primary_fuel_table):

    # create a table with the unique plantids from both dataframes
    eia_plants = eia923_allocated[
        ["plant_id_eia", "plant_primary_fuel"]
    ].drop_duplicates()
    cems_plants = cems[["plant_id_eia"]].drop_duplicates()

    # merge primary fuel into cems
    cems_plants = cems_plants.merge(
        primary_fuel_table.drop_duplicates(subset="plant_id_eia")[
            ["plant_id_eia", "plant_primary_fuel"]
        ],
        how="left",
        on="plant_id_eia",
        validate="1:1",
    )

    # identify any CEMS-only plants that are missing a primary fuel assignment
    plants_missing_primary_fuel = list(
        cems_plants.loc[cems_plants["plant_primary_fuel"].isna(), "plant_id_eia"]
    )

    # calculate primary fuel for each of these missing plants
    cems_primary_fuel = cems.loc[
        cems["plant_id_eia"].isin(plants_missing_primary_fuel),
        ["plant_id_eia", "energy_source_code", "fuel_consumed_mmbtu"],
    ]

    # calculate the total fuel consumption by fuel type and keep the fuel code with the largest fuel consumption
    cems_primary_fuel = (
        cems_primary_fuel.groupby(["plant_id_eia", "energy_source_code"], dropna=False)
        .sum()
        .reset_index()
        .sort_values(by="fuel_consumed_mmbtu", ascending=False)
        .drop_duplicates(subset="plant_id_eia", keep="first")
        .drop(columns="fuel_consumed_mmbtu")
    )

    # merge the cems primary fuel back in and use it to fill any missing fuel codes
    cems_plants = cems_plants.merge(
        cems_primary_fuel, how="left", on="plant_id_eia", validate="1:1",
    )
    cems_plants["plant_primary_fuel"] = cems_plants["plant_primary_fuel"].fillna(
        cems_plants["energy_source_code"]
    )
    cems_plants = cems_plants.drop(columns="energy_source_code")

    # concat the two lists together
    plant_attributes = eia_plants.merge(
        cems_plants,
        how="outer",
        on="plant_id_eia",
        indicator="data_availability",
        suffixes=(None, "_cems"),
    )
    plant_attributes["plant_primary_fuel"] = plant_attributes[
        "plant_primary_fuel"
    ].fillna(plant_attributes["plant_primary_fuel_cems"])
    plant_attributes = plant_attributes.drop(columns=["plant_primary_fuel_cems"])
    plant_attributes["data_availability"] = plant_attributes[
        "data_availability"
    ].replace(
        {"left_only": "eia_only", "right_only": "cems_only", "both": "cems_and_eia"}
    )

    # assign a BA code and state code to each plant
    plant_attributes = assign_ba_code_to_plant(plant_attributes, year)

    # add a flag about whether the plant is distribution connected
    plant_attributes = identify_distribution_connected_plants(
        plant_attributes, year, voltage_threshold_kv=60
    )

    # assign a fuel category to each plant based on what is most likely to match with the category used in EIA-930
    plant_attributes = assign_fuel_category_to_ESC(
        df=plant_attributes, esc_column="plant_primary_fuel",
    )

    # add tz info
    plant_attributes = add_plant_local_timezone(plant_attributes, year)

    plant_attributes = apply_dtypes(plant_attributes)

    return plant_attributes


def assign_ba_code_to_plant(df, year):
    """
    Assigns a balancing authority code and state to each plant based on the plant id
    Inputs:
        df: a pandas dataframe containing a 'plant_id_eia' column
        year: four digit year number for the data
    Returns:
        df with a new column for 'ba_code' and 'state'
    """

    plant_ba = create_plant_ba_table(year)[
        ["plant_id_eia", "ba_code", "ba_code_physical", "state"]
    ]

    # merge the ba code into the dataframe
    df = df.merge(plant_ba, how="left", on="plant_id_eia",)

    if len(df[df["ba_code"].isna()]) > 0:
        print("   Warning: the following plants are missing ba_code:")
        print(df[df["ba_code"].isna()])

    # replace missing ba codes with NA
    df["ba_code"] = df["ba_code"].fillna("NA")
    df["ba_code_physical"] = df["ba_code_physical"].fillna("NA")

    return df


def create_plant_ba_table(year):
    """
    Creates a table assigning a ba_code and physical ba code to each plant id.
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

    # add plants from the plants_entity table in case any are missing from EIA-860
    plants_entity_ba = load_data.load_pudl_table("plants_entity_eia")[
        ["plant_id_eia", "balancing_authority_code_eia", "state"]
    ]
    plant_ba = plant_ba.merge(
        plants_entity_ba, how="outer", on="plant_id_eia", suffixes=(None, "_entity")
    )
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(plant_ba["balancing_authority_code_eia_entity"])
    plant_ba["state"] = plant_ba["state"].fillna(plant_ba["state_entity"].astype(str))
    plant_ba["balancing_authority_code_eia"] = plant_ba[
        "balancing_authority_code_eia"
    ].fillna(value=np.NaN)

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

    # for plants without a BA code assign the miscellaneous BA code based on the state
    plant_ba["ba_code"] = plant_ba["ba_code"].fillna(plant_ba["state"] + "MS")

    # add a physical ba code based on the owner of the transmission system
    plant_ba["ba_code_physical"] = plant_ba["ba_code"]
    plant_ba["ba_code_physical"].update(
        plant_ba["transmission_distribution_owner_name"].map(utility_as_ba_code)
    )

    # update based on mapping table when ambiguous
    physical_ba = pd.read_csv("../data/manual/physical_ba.csv", dtype=get_dtypes())
    plant_ba = plant_ba.merge(
        physical_ba,
        how="left",
        on=["ba_code", "transmission_distribution_owner_name"],
        suffixes=("", "_map"),
    )
    plant_ba["ba_code_physical"].update(plant_ba["ba_code_physical_map"])

    return plant_ba


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
    df,
    fuel_category_names=["fuel_category", "fuel_category_eia930"],
    esc_column="energy_source_code",
):
    """
    Assigns a fuel category to each energy source code in a dataframe.
    Args:
        df: pandas dataframe with column name that matches fuel_category_name and contains energy source codes
        fuel_category_name: list of the columns in energy_source_groups.csv that contains the desired category mapping
        esc_column: name of the column in df that contains the energy source codes to assign a category to
    Returns:
        df with additional column for fuel category
    """
    # load the fuel category table
    energy_source_groups = pd.read_csv(
        "../data/manual/energy_source_groups.csv", dtype=get_dtypes()
    )[["energy_source_code"] + fuel_category_names].rename(
        columns={"energy_source_code": esc_column}
    )
    # assign a fuel category to the monthly eia data
    df = df.merge(
        energy_source_groups[[esc_column] + fuel_category_names],
        how="left",
        on=esc_column,
    )

    return df


def add_plant_local_timezone(df, year):
    pudl_out = load_data.initialize_pudl_out(year=year)
    plant_tz = pudl_out.plants_eia860()[["plant_id_eia", "timezone"]]
    df = df.merge(plant_tz, how="left", on=["plant_id_eia"])

    return df


def aggregate_cems_to_subplant(cems):

    GROUPBY_COLUMNS = ["plant_id_eia", "subplant_id", "datetime_utc", "report_date"]

    DATA_COLUMNS = [
        "gross_generation_mwh",
        "net_generation_mwh",
        "steam_load_1000_lb",
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

    gtn_methods = cems[
        ["plant_id_eia", "subplant_id", "report_date", "gtn_method"]
    ].drop_duplicates()

    cems = cems.groupby(GROUPBY_COLUMNS, dropna=False).sum()[DATA_COLUMNS].reset_index()

    cems = cems.merge(
        gtn_methods,
        how="left",
        on=["plant_id_eia", "subplant_id", "report_date"],
        validate="m:1",
    )

    cems = apply_dtypes(cems)

    return cems

