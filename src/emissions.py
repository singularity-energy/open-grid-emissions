import pandas as pd
import numpy as np


import load_data
import validation
from column_checks import get_dtypes
from filepaths import *

from pudl.analysis.allocate_net_gen import distribute_annually_reported_data_to_months

CLEAN_FUELS = ["SUN", "MWH", "WND", "WAT", "WH", "PUR", "NUC"]


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
    df = df.merge(emission_factors, how="left", on="energy_source_code", validate="m:1")

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
            validate="m:1",
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
        df = df.merge(geothermal_efs, how="left", on=["plant_id_eia"], validate="m:1")

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
        f"{manual_folder()}geothermal_emission_factors.csv",
        dtype=get_dtypes(),
    ).loc[
        :, ["geotype_code", "co2_lb_per_mmbtu", "nox_lb_per_mmbtu", "so2_lb_per_mmbtu"]
    ]

    geothermal_geotypes = identify_geothermal_generator_geotype(year)

    # merge in the emission factor
    geo_efs = geothermal_geotypes.merge(
        geothermal_efs, how="left", on="geotype_code", validate="m:1"
    )

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
        total_plant_capacity, how="left", on="plant_id_eia", validate="m:1"
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


def adjust_fuel_and_emissions_for_CHP(df):
    """Allocates total emissions for electricity generation."""

    # calculate the electric allocation factor
    df = calculate_electric_allocation_factor(df)

    # overwrite the fuel consumed for electricity value using the allocation factor
    df["fuel_consumed_for_electricity_mmbtu"] = (
        df["fuel_consumed_mmbtu"] * df["electric_allocation_factor"]
    )

    if "co2_mass_lb" in df.columns:
        df["co2_mass_lb_for_electricity"] = (
            df["co2_mass_lb"] * df["electric_allocation_factor"]
        )
    if "co2_mass_lb_adjusted" in df.columns:
        df["co2_mass_lb_for_electricity_adjusted"] = (
            df["co2_mass_lb_adjusted"] * df["electric_allocation_factor"]
        )

    if "ch4_mass_lb" in df.columns:
        df["ch4_mass_lb_for_electricity"] = (
            df["ch4_mass_lb"] * df["electric_allocation_factor"]
        )
    if "ch4_mass_lb_adjusted" in df.columns:
        df["ch4_mass_lb_for_electricity_adjusted"] = (
            df["ch4_mass_lb_adjusted"] * df["electric_allocation_factor"]
        )

    if "n2o_mass_lb" in df.columns:
        df["n2o_mass_lb_for_electricity"] = (
            df["n2o_mass_lb"] * df["electric_allocation_factor"]
        )
    if "n2o_mass_lb_adjusted" in df.columns:
        df["n2o_mass_lb_for_electricity_adjusted"] = (
            df["n2o_mass_lb_adjusted"] * df["electric_allocation_factor"]
        )

    if "nox_mass_lb" in df.columns:
        df["nox_mass_lb_for_electricity"] = (
            df["nox_mass_lb"] * df["electric_allocation_factor"]
        )
    if "nox_mass_lb_adjusted" in df.columns:
        df["nox_mass_lb_for_electricity_adjusted"] = (
            df["nox_mass_lb_adjusted"] * df["electric_allocation_factor"]
        )

    if "so2_mass_lb" in df.columns:
        df["so2_mass_lb_for_electricity"] = (
            df["so2_mass_lb"] * df["electric_allocation_factor"]
        )
    if "so2_mass_lb_adjusted" in df.columns:
        df["so2_mass_lb_for_electricity_adjusted"] = (
            df["so2_mass_lb_adjusted"] * df["electric_allocation_factor"]
        )

    df = df.drop(columns=["electric_allocation_factor"])

    validation.test_chp_allocation(df)

    return df


def calculate_electric_allocation_factor(df):
    """
    Calculates an electric allocation factor for CHP plants based on the formula in the eGRID2020 technical guide.

    Requires a dataframe with the following columns: net_generation_mwh, fuel_consumed_mmbtu, fuel_consumed_for_electricity_mmbtu
    """

    mwh_to_mmbtu = 3.412142

    # calculate the useful thermal output
    # 0.8 is an assumed efficiency factor used by eGRID
    df["useful_thermal_output"] = 0.8 * (
        df["fuel_consumed_mmbtu"] - df["fuel_consumed_for_electricity_mmbtu"]
    )

    # convert generation to mmbtu
    df["generation_mmbtu"] = df["net_generation_mwh"] * mwh_to_mmbtu

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

    # CO2: adjust emissions for co2 for all biomass generators
    if "co2_mass_lb" in df.columns:
        df["co2_mass_lb_adjusted"] = df["co2_mass_lb"]
        df.loc[df["energy_source_code"].isin(biomass_fuels), "co2_mass_lb_adjusted"] = 0

    # CH4: for landfill gas (LFG), all other emissions are set to zero
    # this assumes that the gas would have been flared anyway if not used for electricity generation
    if "ch4_mass_lb" in df.columns:
        df["ch4_mass_lb_adjusted"] = df["ch4_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "ch4_mass_lb_adjusted"] = 0

    # N2O: LFG plants set to zero
    if "n2o_mass_lb" in df.columns:
        df["n2o_mass_lb_adjusted"] = df["n2o_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "n2o_mass_lb_adjusted"] = 0

    # NOX: assigned an adjusted value
    # this value is based on using NOx emissions from flaring as a baseline, and subtracting this from the actual emissions
    # to prevent negative emissions, we set the value = 0 if negative
    if "nox_mass_lb" in df.columns:
        df["nox_mass_lb_adjusted"] = df["nox_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "nox_mass_lb_adjusted"] = df.loc[
            df["energy_source_code"] == "LFG", "nox_mass_lb_adjusted"
        ] - (df.loc[df["energy_source_code"] == "LFG", "fuel_consumed_mmbtu"] * 0.078)
        df.loc[df["nox_mass_lb_adjusted"] < 0, "nox_mass_lb_adjusted"] = 0

    # SO2: LFG plants set to zero
    if "so2_mass_lb" in df.columns:
        df["so2_mass_lb_adjusted"] = df["so2_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "so2_mass_lb_adjusted"] = 0

    return df


def calculate_co2e_mass(df, year, gwp_horizon=100, ar5_climate_carbon_feedback=True):
    """
    Calculate CO2-equivalent emissions from CO2, CH4, and N2O. This is done
    by choosing one of the IPCC's emission factors for CH4 and N2O.


    If the `fuel_consumed_for_electricity_units` column is available, we also
    compute the adjusted emissions.
    """
    df_gwp = load_data.load_ipcc_gwp()

    # use the most recent AR that was available in the given year
    ipcc_version = df_gwp.loc[df_gwp["year"] <= year, "year"].max()

    # use the most recent AR that was available in the given year
    most_recent_AR_year = df_gwp.loc[df_gwp["year"] <= year, "year"].max()

    # identify the AR
    ipcc_version = list(
        df_gwp.loc[df_gwp["year"] == most_recent_AR_year, "ipcc_version"].unique()
    )
    if len(ipcc_version) > 1:
        if "AR5" in ipcc_version and ar5_climate_carbon_feedback:
            ipcc_version = "AR5_cc"
        elif "AR5" in ipcc_version and not ar5_climate_carbon_feedback:
            ipcc_version = "AR5_cc"
    else:
        ipcc_version = ipcc_version[0]

    gwp_to_use = df_gwp[df_gwp.ipcc_version == ipcc_version]

    ch4_gwp = gwp_to_use.loc[
        (gwp_to_use.gwp_horizon == gwp_horizon) & (gwp_to_use.gas == "ch4"),
        "gwp",
    ].item()
    n2o_gwp = gwp_to_use.loc[
        (gwp_to_use.gwp_horizon == gwp_horizon) & (gwp_to_use.gas == "n2o"),
        "gwp",
    ].item()

    df["co2e_mass_lb"] = (
        df["co2_mass_lb"] + ch4_gwp * df["ch4_mass_lb"] + n2o_gwp * df["n2o_mass_lb"]
    )

    if "co2_mass_lb_adjusted" in df:
        df["co2e_mass_lb_adjusted"] = (
            df["co2_mass_lb_adjusted"]
            + ch4_gwp * df["ch4_mass_lb_adjusted"]
            + n2o_gwp * df["n2o_mass_lb_adjusted"]
        )
    if "co2_mass_lb_for_electricity" in df:
        df["co2e_mass_lb_for_electricity"] = (
            df["co2_mass_lb_for_electricity"]
            + ch4_gwp * df["ch4_mass_lb_for_electricity"]
            + n2o_gwp * df["n2o_mass_lb_for_electricity"]
        )
    if "co2_mass_lb_for_electricity_adjusted" in df:
        df["co2e_mass_lb_for_electricity_adjusted"] = (
            df["co2_mass_lb_for_electricity_adjusted"]
            + ch4_gwp * df["ch4_mass_lb_for_electricity_adjusted"]
            + n2o_gwp * df["n2o_mass_lb_for_electricity_adjusted"]
        )

    return df


def calculate_nox_from_fuel_consumption(
    gen_fuel_allocated: pd.DataFrame, pudl_out, year
) -> pd.DataFrame:
    """
    Calculate NOx emissions from fuel consumption data.

    Apply emission rates in the following order based on availability:
       1. Month- and generator-specific uncontrolled NOx emission factors based on boiler design parameters of all associated boilers and heat content of consumed fuels
       2. Season- and generator-specific controlled NOx emissions factors based on boiler-specific NOx control equipment in operation.

    """

    # calculate uncontrolled nox emission factors based on boiler design parameters
    uncontrolled_nox_factors = calculate_generator_nox_ef_per_unit_from_boiler_type(
        gen_fuel_allocated, year, pudl_out
    )
    uncontrolled_nox_factors = convert_ef_to_lb_per_mmbtu(
        uncontrolled_nox_factors, pudl_out, "nox"
    )

    # merge nox efs into the gen_fuel_allocated
    gen_fuel_allocated = gen_fuel_allocated.merge(
        uncontrolled_nox_factors,
        how="left",
        on=[
            "report_date",
            "plant_id_eia",
            "energy_source_code",
            "prime_mover_code",
            "generator_id",
        ],
        validate="m:1",
    )
    # raise a warning if we are missing emission factors for any non-zero fuel consumption from non-clean fuels
    missing_ef = gen_fuel_allocated[
        gen_fuel_allocated["nox_ef_lb_per_mmbtu"].isna()
        & ~gen_fuel_allocated["fuel_consumed_mmbtu"].isna()
        & ~gen_fuel_allocated["energy_source_code"].isin(CLEAN_FUELS)
    ]
    if len(missing_ef) > 0:
        print(
            missing_ef[
                [
                    "report_date",
                    "plant_id_eia",
                    "energy_source_code",
                    "prime_mover_code",
                    "generator_id",
                ]
            ].drop_duplicates()
        )
        raise UserWarning("NOx emission factors are missing for the above records")
    gen_fuel_allocated["nox_mass_lb"] = (
        gen_fuel_allocated["fuel_consumed_mmbtu"]
        * gen_fuel_allocated["nox_ef_lb_per_mmbtu"]
    )

    # calculate the controlled nox rates
    controlled_nox_factors = calculate_generator_specific_controlled_nox_rates(year)

    # merge the controlled nox rates into gen_fuel_allocated
    gen_fuel_allocated = gen_fuel_allocated.merge(
        controlled_nox_factors,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )
    # calculate the controlled nox emissions based on the month
    gen_fuel_allocated = gen_fuel_allocated.assign(
        controlled_nox_mass_lb=lambda x: np.where(
            ((x.report_date.dt.month >= 5) & (x.report_date.dt.month <= 9)),
            x.fuel_consumed_mmbtu * x.controlled_ozone_season_nox_ef_lb_per_mmbtu,
            x.fuel_consumed_mmbtu * x.controlled_non_ozone_season_nox_ef_lb_per_mmbtu,
        )
    )
    # if there were not season-specific rates, fill in using the annual rate if available
    gen_fuel_allocated["controlled_nox_mass_lb"] = gen_fuel_allocated[
        "controlled_nox_mass_lb"
    ].fillna(
        gen_fuel_allocated["fuel_consumed_mmbtu"]
        * gen_fuel_allocated["controlled_annual_nox_ef_lb_per_mmbtu"]
    )
    # update the emision total using the controlled mass if available
    gen_fuel_allocated["nox_mass_lb"].update(
        gen_fuel_allocated["controlled_nox_mass_lb"]
    )

    # remove intermediate columns
    gen_fuel_allocated = gen_fuel_allocated.drop(
        columns=[
            "nox_ef_lb_per_mmbtu",
            "controlled_annual_nox_ef_lb_per_mmbtu",
            "controlled_ozone_season_nox_ef_lb_per_mmbtu",
            "controlled_non_ozone_season_nox_ef_lb_per_mmbtu",
            "controlled_nox_mass_lb",
        ]
    )

    return gen_fuel_allocated


def calculate_generator_nox_ef_per_unit_from_boiler_type(
    gen_fuel_allocated, year, pudl_out
):
    """Calculates a generator-specific Nox emission factor per unit fuel based on boiler firing type

    If a generator has multiple boilers, average the emission factor of all boilers
    """

    # get a dataframe with all unique generator-pm-esc combinations for emitting energy source types with data reported
    gen_keys_for_nox = gen_fuel_allocated.copy()[
        ~gen_fuel_allocated["fuel_consumed_mmbtu"].isna()
    ]
    gen_keys_for_nox = gen_keys_for_nox[
        [
            "report_date",
            "plant_id_eia",
            "generator_id",
            "prime_mover_code",
            "energy_source_code",
        ]
    ].drop_duplicates()
    gen_keys_for_nox = gen_keys_for_nox[
        ~gen_keys_for_nox["energy_source_code"].isin(CLEAN_FUELS)
    ]

    # load emission factors
    nox_emission_factors = load_data.load_nox_emission_factors()
    # remove duplicate factors
    nox_emission_factors = nox_emission_factors.drop_duplicates(
        subset=[
            "prime_mover_code",
            "energy_source_code",
            "boiler_bottom_type",
            "boiler_firing_type",
        ]
    )

    # load the boiler firing type info
    boiler_firing_type = load_boiler_firing_type(year)

    # identify the boiler firing type for each generator
    boiler_generator_assn = pudl_out.bga_eia860()
    # associate each boiler record with generator_id s
    boiler_firing_type = boiler_firing_type.merge(
        boiler_generator_assn[["plant_id_eia", "boiler_id", "generator_id"]],
        how="left",
        on=["plant_id_eia", "boiler_id"],
    )

    # merge the gen keys with the boiler firing types
    gen_nox_factors = gen_keys_for_nox.merge(
        boiler_firing_type,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:m",
    )

    # merge in the emission factors for spedcific boiler types
    gen_nox_factors = gen_nox_factors.merge(
        nox_emission_factors,
        how="left",
        on=[
            "prime_mover_code",
            "energy_source_code",
            "boiler_bottom_type",
            "boiler_firing_type",
        ],
        validate="m:1",
    )

    # identify missing emission factors and replace with zeros
    missing_nox_efs = (
        gen_nox_factors.loc[
            gen_nox_factors["emission_factor"].isna(),
            [
                "prime_mover_code",
                "energy_source_code",
                "boiler_bottom_type",
                "boiler_firing_type",
            ],
        ]
        .drop_duplicates()
        .sort_values(
            by=[
                "energy_source_code",
                "prime_mover_code",
                "boiler_firing_type",
                "boiler_bottom_type",
            ]
        )
    )
    if len(missing_nox_efs) > 0:
        print(" ")
        print(
            "Warning: NOx emission factors are missing for the following boiler types. An emission factor of zero will be used for these boilers."
        )
        print(missing_nox_efs)
        print(" ")
    gen_nox_factors["emission_factor"] = gen_nox_factors["emission_factor"].fillna(0)

    # average the emission factors for all boilers associated with each generator
    gen_nox_factors = (
        gen_nox_factors.groupby(
            [
                "report_date",
                "plant_id_eia",
                "generator_id",
                "prime_mover_code",
                "energy_source_code",
                "emission_factor_denominator",
            ],
            dropna=False,
        )["emission_factor"]
        .mean()
        .reset_index()
    )

    return gen_nox_factors


def load_boiler_firing_type(year):

    boiler_design_parameters_eia860 = load_data.load_boiler_design_parameters_eia860(
        year
    )

    firing_type_description = {
        "CB": "CELLBURNER",
        "CY": "CYCLONE",
        "DB": "DUCTBURNER",
        "FB": "FLUIDIZED",
        "SS": "STOKER",
        "TF": "TANGENTIAL",
        "VF": "VERTICAL",
        "WF": "WALL",
        "OT": "OTHER",
    }

    # only keep boilers that are operational
    boiler_firing_type = boiler_design_parameters_eia860.copy()[
        boiler_design_parameters_eia860["operational_status"] == "OP"
    ]

    boiler_firing_type["boiler_firing_type"] = boiler_firing_type["firing_type_1"].map(
        firing_type_description
    )

    boiler_firing_type["boiler_bottom_type"] = boiler_firing_type[
        "boiler_bottom_type"
    ].replace({"D": "DRY", "W": "WET"})

    boiler_firing_type = boiler_firing_type[
        ["plant_id_eia", "boiler_id", "boiler_bottom_type", "boiler_firing_type"]
    ].dropna(subset=["boiler_bottom_type", "boiler_firing_type"], thresh=1)

    return boiler_firing_type


def convert_ef_to_lb_per_mmbtu(gen_emission_factors, pudl_out, pollutant):

    # get the reported fuel heat content values from EIA-923
    (
        plant_specific_fuel_heat_content,
        national_avg_fuel_heat_content,
        annual_avg_fuel_heat_content,
    ) = return_monthly_plant_fuel_heat_content(pudl_out)

    # merge in plant pm specific fuel heat content values
    gen_emission_factors = gen_emission_factors.merge(
        plant_specific_fuel_heat_content,
        how="left",
        on=["report_date", "plant_id_eia", "energy_source_code", "prime_mover_code"],
        validate="m:1",
    )
    # merge in national monthly average fuel heat content values and use to fill missing values
    gen_emission_factors = gen_emission_factors.merge(
        national_avg_fuel_heat_content,
        how="left",
        on=[
            "report_date",
            "energy_source_code",
        ],
        validate="m:1",
        suffixes=(None, "_national"),
    )
    gen_emission_factors["fuel_mmbtu_per_unit"] = gen_emission_factors[
        "fuel_mmbtu_per_unit"
    ].fillna(gen_emission_factors["fuel_mmbtu_per_unit_national"])
    # merge in annual average fuel heat content and use to fill missing values
    gen_emission_factors = gen_emission_factors.merge(
        annual_avg_fuel_heat_content,
        how="left",
        on=[
            "energy_source_code",
        ],
        validate="m:1",
        suffixes=(None, "_annual"),
    )
    gen_emission_factors["fuel_mmbtu_per_unit"] = gen_emission_factors[
        "fuel_mmbtu_per_unit"
    ].fillna(gen_emission_factors["fuel_mmbtu_per_unit_annual"])

    # check to make sure that we have fuel content data for all plants that need it
    missing_fuel_content = gen_emission_factors[
        gen_emission_factors["fuel_mmbtu_per_unit"].isna()
        & (gen_emission_factors["emission_factor_denominator"] != "mmbtu")
    ]
    if len(missing_fuel_content) > 0:
        raise UserWarning(
            f"Warning: The heat content for the following fuels is missing and NOx emissions will not be calculated for this fuel:{list(missing_fuel_content.energy_source_code.unique())}"
        )

    # convert emission factors from lb per unit to lb per mmbtu if the factor is not already in units of lb/mmbtu
    if pollutant == "nox":
        gen_emission_factors = gen_emission_factors.assign(
            nox_ef_lb_per_mmbtu=lambda x: np.where(
                (x.emission_factor_denominator == "mmbtu"),
                x.emission_factor,
                (x.emission_factor / x.fuel_mmbtu_per_unit),
            )
        )
    elif pollutant == "so2":
        gen_emission_factors = gen_emission_factors.assign(
            so2_ef_lb_per_mmbtu=lambda x: np.where(
                (x.emission_factor_denominator == "mmbtu"),
                x.emission_factor,
                (x.emission_factor / x.fuel_mmbtu_per_unit),
            )
        )
    else:
        raise UserWarning(
            f"arg `pollutant` must be 'nox' or 'so2'. You specified '{pollutant}'"
        )

    # drop intermediate columns
    gen_emission_factors = gen_emission_factors.drop(
        columns=[
            "emission_factor_denominator",
            "emission_factor",
            "fuel_mmbtu_per_unit",
            "fuel_mmbtu_per_unit_national",
            "fuel_mmbtu_per_unit_annual",
        ]
    )

    return gen_emission_factors


def return_monthly_plant_fuel_heat_content(pudl_out):
    # load information about the monthly heat input of fuels
    plant_specific_fuel_heat_content = (
        pudl_out.gf_eia923()
        .loc[
            :,
            [
                "plant_id_eia",
                "energy_source_code",
                "prime_mover_code",
                "report_date",
                "fuel_mmbtu_per_unit",
            ],
        ]
        .pipe(
            distribute_annually_reported_data_to_months,
            key_columns=["plant_id_eia", "energy_source_code", "prime_mover_code"],
            data_column_name="fuel_mmbtu_per_unit",
        )
    )
    plant_specific_fuel_heat_content = plant_specific_fuel_heat_content[
        ~plant_specific_fuel_heat_content["energy_source_code"].isin(CLEAN_FUELS)
    ]

    # replace zero heat content with missing values
    plant_specific_fuel_heat_content[
        "fuel_mmbtu_per_unit"
    ] = plant_specific_fuel_heat_content["fuel_mmbtu_per_unit"].replace(0, np.NaN)
    plant_specific_fuel_heat_content[
        "fuel_mmbtu_per_unit"
    ] = plant_specific_fuel_heat_content["fuel_mmbtu_per_unit"].replace(np.inf, np.NaN)

    # calculate the average monthly heat content for a fuel
    national_avg_fuel_heat_content = (
        plant_specific_fuel_heat_content.drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "report_date"], dropna=False)
        .mean()
        .reset_index()
    )

    annual_avg_fuel_heat_content = (
        national_avg_fuel_heat_content.groupby(["energy_source_code"], dropna=False)
        .mean()
        .reset_index()
    )

    # change the report date columns back to datetimes
    plant_specific_fuel_heat_content["report_date"] = pd.to_datetime(
        plant_specific_fuel_heat_content["report_date"]
    )
    national_avg_fuel_heat_content["report_date"] = pd.to_datetime(
        national_avg_fuel_heat_content["report_date"]
    )

    return (
        plant_specific_fuel_heat_content,
        national_avg_fuel_heat_content,
        annual_avg_fuel_heat_content,
    )


def calculate_unit_specific_controlled_nox_rates(year):
    nox_rates = load_controlled_nox_emission_rates(year)
    nox_rates = calculate_non_ozone_season_nox_rate(nox_rates)
    weighted_nox_rates = calculate_weighted_nox_rates(year, nox_rates, "unitid")

    return weighted_nox_rates


def calculate_boiler_specific_controlled_nox_rates(year):
    nox_rates = load_controlled_nox_emission_rates(year)
    nox_rates = calculate_non_ozone_season_nox_rate(nox_rates)
    weighted_nox_rates = calculate_weighted_nox_rates(year, nox_rates, "boiler_id")

    return weighted_nox_rates


def calculate_generator_specific_controlled_nox_rates(year):
    nox_rates = load_controlled_nox_emission_rates(year)
    nox_rates = calculate_non_ozone_season_nox_rate(nox_rates)
    weighted_nox_rates = calculate_weighted_nox_rates(year, nox_rates, "generator_id")

    return weighted_nox_rates


def load_controlled_nox_emission_rates(year):
    # load the emissions control data
    emissions_controls_eia923 = load_data.load_emissions_controls_eia923(year)

    # create a dataframe that contains only operating NOx emission control data
    nox_rates = emissions_controls_eia923[
        ~emissions_controls_eia923["nox_control_id"].isna()
    ]
    nox_rates = nox_rates[nox_rates["operational_status"] == "OP"]
    nox_rates = nox_rates[
        [
            "plant_id_eia",
            "nox_control_id",
            "hours_in_service",
            "annual_nox_emission_rate_lb_per_mmbtu",
            "ozone_season_nox_emission_rate_lb_per_mmbtu",
        ]
    ].dropna(
        subset=[
            "annual_nox_emission_rate_lb_per_mmbtu",
            "ozone_season_nox_emission_rate_lb_per_mmbtu",
        ],
        thresh=1,
    )
    nox_rates = nox_rates.rename(
        columns={
            "annual_nox_emission_rate_lb_per_mmbtu": "controlled_annual_nox_ef_lb_per_mmbtu",
            "ozone_season_nox_emission_rate_lb_per_mmbtu": "controlled_ozone_season_nox_ef_lb_per_mmbtu",
        }
    )

    return nox_rates


def calculate_weighted_nox_rates(year, nox_rates, aggregation_level):
    """Aggregates nox rate data from nox_control_id to boiler_id, generator_id, or unitid"""
    # load the association tables
    boiler_nox_association_eia860 = load_data.load_boiler_nox_association_eia860(year)

    # merge boiler ids associated with each nox_control_id
    nox_rates = nox_rates.merge(
        boiler_nox_association_eia860[
            [
                "plant_id_eia",
                "nox_control_id",
                "boiler_id",
            ]
        ],
        how="left",
        on=["plant_id_eia", "nox_control_id"],
        validate="m:m",
    )

    # if there are any missing boiler_ids, fill using the nox_control_id, which is likely to match a boiler
    nox_rates["boiler_id"] = nox_rates["boiler_id"].fillna(nox_rates["nox_control_id"])

    if aggregation_level == "generator_id":
        boiler_generator_assn = load_data.initialize_pudl_out(year).bga_eia860()
        # associate a generator_id with each record
        nox_rates = nox_rates.merge(
            boiler_generator_assn[["plant_id_eia", "boiler_id", "generator_id"]],
            how="left",
            on=["plant_id_eia", "boiler_id"],
            validate="m:m",
        )
    elif aggregation_level == "unitid":
        epa_eia_crosswalk = load_data.load_epa_eia_crosswalk(year)
        boiler_to_unit_crosswalk = epa_eia_crosswalk[
            ["plant_id_eia", "unitid", "boiler_id"]
        ].drop_duplicates()
        # associate a unitid with each record
        nox_rates = nox_rates.merge(
            boiler_to_unit_crosswalk,
            how="left",
            on=["plant_id_eia", "boiler_id"],
            validate="m:m",
        )

    # calculate a weighted average for each boiler or generator
    weighted_nox_rates = calculate_weighted_averages(
        nox_rates,
        groupby_columns=["plant_id_eia", aggregation_level],
        data_cols=[
            "controlled_annual_nox_ef_lb_per_mmbtu",
            "controlled_ozone_season_nox_ef_lb_per_mmbtu",
            "controlled_non_ozone_season_nox_ef_lb_per_mmbtu",
        ],
        weight_col="hours_in_service",
    )
    return weighted_nox_rates


def calculate_weighted_averages(df, groupby_columns, data_cols, weight_col):
    """helper function for calculating weighted averages of one or more columns in a dataframe."""
    wa = df.copy()
    for data_col in data_cols:
        wa[f"{data_col}_data_times_weight"] = wa[data_col] * wa[weight_col]
        wa[f"{data_col}_weight_where_notnull"] = wa[weight_col] * pd.notnull(
            wa[data_col]
        )
    g = wa.groupby(groupby_columns, dropna=False)
    result = pd.DataFrame()
    for data_col in data_cols:
        result[data_col] = (
            g[f"{data_col}_data_times_weight"].sum()
            / g[f"{data_col}_weight_where_notnull"].sum()
        )
    result = result.reset_index()

    return result


def calculate_non_ozone_season_nox_rate(weighted_nox_rates):
    annual_col = "controlled_annual_nox_ef_lb_per_mmbtu"
    oz_col = "controlled_ozone_season_nox_ef_lb_per_mmbtu"
    non_oz_col = "controlled_non_ozone_season_nox_ef_lb_per_mmbtu"

    # ozone season is May - Sept (5 months).
    # To get the average emission rate for the 7 non-ozone season months, we assume similar operation across all months
    # annual_avg = [(5* oz_avg) + (7 * non_oz_avg)] / 12
    weighted_nox_rates[non_oz_col] = (
        (12 * weighted_nox_rates[annual_col]) - (5 * weighted_nox_rates[oz_col])
    ) / 7

    # if there are any rates that we calculate as negative, replace with zero
    weighted_nox_rates.loc[weighted_nox_rates[non_oz_col] < 0, non_oz_col] = 0

    return weighted_nox_rates


def calculate_so2_from_fuel_consumption(gen_fuel_allocated, pudl_out, year):
    """
    Calculate SO2 emissions from fuel consumption data.

    Apply emission rates in the following order based on availability:
       1. Month- and generator-specific uncontrolled SO2 emission factors based on boiler design parameters of all associated boilers and heat content of consumed fuels
       2. Adjust SO2 emissions using boiler-specific SO2 removal efficiencies if available. Otherwise, assume SO2 emissions are uncontrolled.

    """

    # load uncontrolled emission factors
    uncontrolled_so2_factors = calculate_generator_so2_ef_per_unit_from_boiler_type(
        gen_fuel_allocated, year, pudl_out
    )

    # adjust for sulfur content
    uncontrolled_so2_factors = adjust_so2_efs_for_fuel_sulfur_content(
        uncontrolled_so2_factors, pudl_out
    )

    # convert to mmbtu
    uncontrolled_so2_factors = convert_ef_to_lb_per_mmbtu(
        uncontrolled_so2_factors, pudl_out, "so2"
    )

    # merge so2 efs into the gen_fuel_allocated
    gen_fuel_allocated = gen_fuel_allocated.merge(
        uncontrolled_so2_factors,
        how="left",
        on=[
            "report_date",
            "plant_id_eia",
            "energy_source_code",
            "prime_mover_code",
            "generator_id",
        ],
        validate="m:1",
    )
    # raise a warning if we are missing emission factors for any non-zero fuel consumption from non-clean fuels
    missing_ef = gen_fuel_allocated[
        gen_fuel_allocated["so2_ef_lb_per_mmbtu"].isna()
        & ~gen_fuel_allocated["fuel_consumed_mmbtu"].isna()
        & ~gen_fuel_allocated["energy_source_code"].isin(CLEAN_FUELS)
    ]
    if len(missing_ef) > 0:
        print(
            missing_ef[
                [
                    "report_date",
                    "plant_id_eia",
                    "energy_source_code",
                    "prime_mover_code",
                    "generator_id",
                ]
            ].drop_duplicates()
        )
        raise UserWarning("SO2 emission factors are missing for the above records")
    gen_fuel_allocated["so2_mass_lb"] = (
        gen_fuel_allocated["fuel_consumed_mmbtu"]
        * gen_fuel_allocated["so2_ef_lb_per_mmbtu"]
    )

    # load so2 removal efficiency information
    so2_control_efficiency = load_so2_control_efficiencies(year)
    so2_control_efficiency = calculate_weighted_so2_control_efficiency(
        year, so2_control_efficiency, "generator_id"
    )

    # calculate the controlled so2 emissions
    gen_fuel_allocated = gen_fuel_allocated.merge(
        so2_control_efficiency,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )
    # assume all other generators have uncontrolled emissions
    gen_fuel_allocated["so2_removal_efficiency_annual"] = gen_fuel_allocated[
        "so2_removal_efficiency_annual"
    ].fillna(0)

    # calculate controlled so2 emissions
    gen_fuel_allocated["so2_mass_lb"] = gen_fuel_allocated["so2_mass_lb"] * (
        1 - gen_fuel_allocated["so2_removal_efficiency_annual"]
    )

    # remove intermediate columns
    gen_fuel_allocated = gen_fuel_allocated.drop(
        columns=[
            "so2_ef_lb_per_mmbtu",
            "so2_removal_efficiency_annual",
            "so2_removal_efficiency_at_full_load",
        ]
    )

    return gen_fuel_allocated


def calculate_generator_so2_ef_per_unit_from_boiler_type(
    gen_fuel_allocated, year, pudl_out
):
    """Calculates a generator-specific Nox emission factor per unit fuel based on boiler firing type

    If a generator has multiple boilers, average the emission factor of all boilers
    """

    # get a dataframe with all unique generator-pm-esc combinations for emitting energy source types with data reported
    gen_keys_for_so2 = gen_fuel_allocated.copy()[
        ~gen_fuel_allocated["fuel_consumed_mmbtu"].isna()
    ]
    gen_keys_for_so2 = gen_keys_for_so2[
        [
            "report_date",
            "plant_id_eia",
            "generator_id",
            "prime_mover_code",
            "energy_source_code",
        ]
    ].drop_duplicates()
    gen_keys_for_so2 = gen_keys_for_so2[
        ~gen_keys_for_so2["energy_source_code"].isin(CLEAN_FUELS)
    ]

    # load emission factors
    so2_emission_factors = load_data.load_so2_emission_factors()

    # load the boiler firing type info
    boiler_firing_type = load_boiler_firing_type(year)
    # drop the boiler bottom type data
    boiler_firing_type = boiler_firing_type.drop(columns="boiler_bottom_type")
    boiler_firing_type = boiler_firing_type.drop_duplicates()

    # identify the boiler firing type for each generator
    boiler_generator_assn = pudl_out.bga_eia860()
    # associate each boiler record with generator_id s
    boiler_firing_type = boiler_firing_type.merge(
        boiler_generator_assn[["plant_id_eia", "boiler_id", "generator_id"]],
        how="left",
        on=["plant_id_eia", "boiler_id"],
        validate="1:m",
    )

    # merge the gen keys with the boiler firing types
    gen_so2_factors = gen_keys_for_so2.merge(
        boiler_firing_type,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:m",
    )

    # merge in the emission factors for specific boiler types
    gen_so2_factors = gen_so2_factors.merge(
        so2_emission_factors,
        how="left",
        on=[
            "prime_mover_code",
            "energy_source_code",
            "boiler_firing_type",
        ],
        validate="m:1",
    )

    # fill in geotype-specific geothermal emission factors
    if gen_so2_factors["energy_source_code"].str.contains("GEO").any():
        gen_so2_factors = add_geothermal_emission_factors(
            gen_so2_factors,
            year,
            include_co2=False,
            include_nox=False,
            include_so2=True,
        )
        gen_so2_factors.loc[
            gen_so2_factors["energy_source_code"] == "GEO", "emission_factor"
        ] = gen_so2_factors.loc[
            gen_so2_factors["energy_source_code"] == "GEO", "so2_lb_per_mmbtu"
        ]
        gen_so2_factors.loc[
            gen_so2_factors["energy_source_code"] == "GEO", "emission_factor_numerator"
        ] = "lb"
        gen_so2_factors.loc[
            gen_so2_factors["energy_source_code"] == "GEO",
            "emission_factor_denominator",
        ] = "mmbtu"
        gen_so2_factors.loc[
            gen_so2_factors["energy_source_code"] == "GEO", "multiply_by_sulfur_content"
        ] = 0
        gen_so2_factors = gen_so2_factors.drop(columns="so2_lb_per_mmbtu")

    # identify missing emission factors and replace with zeros
    missing_so2_efs = (
        gen_so2_factors.loc[
            gen_so2_factors["emission_factor"].isna(),
            [
                "prime_mover_code",
                "energy_source_code",
                "boiler_firing_type",
            ],
        ]
        .drop_duplicates()
        .sort_values(
            by=[
                "energy_source_code",
                "prime_mover_code",
                "boiler_firing_type",
            ]
        )
    )
    if len(missing_so2_efs) > 0:
        print(" ")
        print(
            "Warning: SO2 emission factors are missing for the following boiler types. An emission factor of zero will be used for these boilers."
        )
        print(missing_so2_efs)
        print(" ")
    gen_so2_factors["emission_factor"] = gen_so2_factors["emission_factor"].fillna(0)
    gen_so2_factors["multiply_by_sulfur_content"] = gen_so2_factors[
        "multiply_by_sulfur_content"
    ].fillna(0)

    # average the emission factors for all boilers associated with each generator
    gen_so2_factors = (
        gen_so2_factors.groupby(
            [
                "report_date",
                "plant_id_eia",
                "generator_id",
                "prime_mover_code",
                "energy_source_code",
                "emission_factor_denominator",
                "multiply_by_sulfur_content",
            ],
            dropna=False,
        )["emission_factor"]
        .mean()
        .reset_index()
    )

    return gen_so2_factors


def return_monthly_plant_fuel_sulfur_content(pudl_out):
    plant_specific_fuel_sulfur_content = pudl_out.bf_eia923().loc[
        :,
        [
            "plant_id_eia",
            "boiler_id",
            "energy_source_code",
            "report_date",
            "fuel_consumed_units",
            "sulfur_content_pct",
        ],
    ]

    # divide by 100 to get decimal percent
    plant_specific_fuel_sulfur_content["sulfur_content_pct"] = (
        plant_specific_fuel_sulfur_content["sulfur_content_pct"] / 100
    )

    plant_specific_fuel_sulfur_content = plant_specific_fuel_sulfur_content[
        ~plant_specific_fuel_sulfur_content["energy_source_code"].isin(CLEAN_FUELS)
    ]

    # merge in the prime mover data
    plant_specific_fuel_sulfur_content = plant_specific_fuel_sulfur_content.merge(
        pd.read_sql("boilers_entity_eia", pudl_out.pudl_engine),
        how="left",
        on=["plant_id_eia", "boiler_id"],
        validate="m:1",
    )

    # calculate a weighted average for each  generator
    plant_specific_fuel_sulfur_content = calculate_weighted_averages(
        plant_specific_fuel_sulfur_content,
        groupby_columns=[
            "plant_id_eia",
            "energy_source_code",
            "prime_mover_code",
            "report_date",
        ],
        data_cols=[
            "sulfur_content_pct",
        ],
        weight_col="fuel_consumed_units",
    )

    # calculate the average monthly heat content for a fuel
    national_avg_fuel_sulfur_content = (
        plant_specific_fuel_sulfur_content.drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "report_date"], dropna=False)
        .mean()
        .reset_index()
    )

    annual_avg_fuel_sulfur_content = (
        national_avg_fuel_sulfur_content.groupby(["energy_source_code"], dropna=False)
        .mean()
        .reset_index()
    )

    # change the report date columns back to datetimes
    plant_specific_fuel_sulfur_content["report_date"] = pd.to_datetime(
        plant_specific_fuel_sulfur_content["report_date"]
    )
    national_avg_fuel_sulfur_content["report_date"] = pd.to_datetime(
        national_avg_fuel_sulfur_content["report_date"]
    )

    plant_specific_fuel_sulfur_content

    return (
        plant_specific_fuel_sulfur_content,
        national_avg_fuel_sulfur_content,
        annual_avg_fuel_sulfur_content,
    )


def adjust_so2_efs_for_fuel_sulfur_content(uncontrolled_so2_factors, pudl_out):
    # multiply factors by sulfur content
    (
        plant_specific_fuel_sulfur_content,
        national_avg_fuel_sulfur_content,
        annual_avg_fuel_sulfur_content,
    ) = return_monthly_plant_fuel_sulfur_content(pudl_out)

    # merge in plant pm specific fuel sulfur content values
    uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
        plant_specific_fuel_sulfur_content,
        how="left",
        on=["report_date", "plant_id_eia", "energy_source_code", "prime_mover_code"],
        validate="m:1",
    )
    # merge in national monthly average fuel sulfur content values and use to fill missing values
    uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
        national_avg_fuel_sulfur_content,
        how="left",
        on=[
            "report_date",
            "energy_source_code",
        ],
        validate="m:1",
        suffixes=(None, "_national"),
    )
    uncontrolled_so2_factors["sulfur_content_pct"] = uncontrolled_so2_factors[
        "sulfur_content_pct"
    ].fillna(uncontrolled_so2_factors["sulfur_content_pct_national"])
    # merge in annual average fuel sulfur content and use to fill missing values
    uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
        annual_avg_fuel_sulfur_content,
        how="left",
        on=[
            "energy_source_code",
        ],
        validate="m:1",
        suffixes=(None, "_annual"),
    )
    uncontrolled_so2_factors["sulfur_content_pct"] = uncontrolled_so2_factors[
        "sulfur_content_pct"
    ].fillna(uncontrolled_so2_factors["sulfur_content_pct_annual"])

    # check that we are not missing sulfur content for any generators that need it
    missing_sulfur_content = uncontrolled_so2_factors[
        uncontrolled_so2_factors["sulfur_content_pct"].isna()
        & (uncontrolled_so2_factors["multiply_by_sulfur_content"] == 1)
    ]
    if len(missing_sulfur_content) > 0:
        print(
            missing_sulfur_content[
                [
                    "plant_id_eia",
                    "generator_id",
                    "prime_mover_code",
                    "energy_source_code",
                ]
            ].drop_duplicates()
        )
        raise UserWarning(
            "Sulfur content data is missing for the above units that require it."
        )
    uncontrolled_so2_factors.loc[
        uncontrolled_so2_factors["sulfur_content_pct"].isna()
        & (uncontrolled_so2_factors["multiply_by_sulfur_content"] == 0),
        "sulfur_content_pct",
    ] = 0

    # perform the adjustment
    uncontrolled_so2_factors = uncontrolled_so2_factors.assign(
        emission_factor=lambda x: np.where(
            (x.multiply_by_sulfur_content == 1),
            (x.emission_factor * x.sulfur_content_pct),
            x.emission_factor,
        )
    )

    # drop intermediate columns
    uncontrolled_so2_factors = uncontrolled_so2_factors.drop(
        columns=[
            "multiply_by_sulfur_content",
            "sulfur_content_pct",
            "sulfur_content_pct_national",
            "sulfur_content_pct_annual",
        ]
    )

    return uncontrolled_so2_factors


def load_so2_control_efficiencies(year):
    # load the emissions control data
    emissions_controls_eia923 = load_data.load_emissions_controls_eia923(year)

    # create a dataframe that contains only operating NOx emission control data
    so2_efficiency = emissions_controls_eia923[
        ~emissions_controls_eia923["so2_control_id"].isna()
    ]
    so2_efficiency = so2_efficiency[so2_efficiency["operational_status"] == "OP"]
    so2_efficiency = so2_efficiency[
        [
            "plant_id_eia",
            "so2_control_id",
            "hours_in_service",
            "so2_removal_efficiency_annual",
            "so2_removal_efficiency_at_full_load",
        ]
    ].dropna(
        subset=[
            "so2_removal_efficiency_annual",
            "so2_removal_efficiency_at_full_load",
        ],
        thresh=1,
    )

    # validate that the efficiencies are between 0 and 1
    bad_efficiencies = so2_efficiency[
        (so2_efficiency["so2_removal_efficiency_annual"] < 0)
        | (so2_efficiency["so2_removal_efficiency_annual"] > 1)
    ]
    if len(bad_efficiencies) > 0:
        raise UserWarning(
            "Warning: certain loaded SO2 removal efficiencies are either negative or > 100%"
        )

    return so2_efficiency


def calculate_weighted_so2_control_efficiency(
    year, so2_control_efficiency, aggregation_level
):
    """Aggregates so2 rate data from so2_control_id to generator_id"""
    # load the association tables
    boiler_so2_association_eia860 = load_data.load_boiler_so2_association_eia860(year)

    # merge boiler ids associated with each so2_control_id
    so2_control_efficiency = so2_control_efficiency.merge(
        boiler_so2_association_eia860[
            [
                "plant_id_eia",
                "so2_control_id",
                "boiler_id",
            ]
        ],
        how="left",
        on=["plant_id_eia", "so2_control_id"],
        validate="m:m",
    )

    # if there are any missing boiler_ids, fill using the so2_control_id, which is likely to match a boiler
    so2_control_efficiency["boiler_id"] = so2_control_efficiency["boiler_id"].fillna(
        so2_control_efficiency["so2_control_id"]
    )

    if aggregation_level == "generator_id":
        boiler_generator_assn = load_data.initialize_pudl_out(year).bga_eia860()
        # associate a generator_id with each record
        so2_control_efficiency = so2_control_efficiency.merge(
            boiler_generator_assn[["plant_id_eia", "boiler_id", "generator_id"]],
            how="left",
            on=["plant_id_eia", "boiler_id"],
            validate="m:m",
        )

    # calculate a weighted average for each boiler or generator
    weighted_so2_control_efficiency = calculate_weighted_averages(
        so2_control_efficiency,
        groupby_columns=["plant_id_eia", aggregation_level],
        data_cols=[
            "so2_removal_efficiency_annual",
            "so2_removal_efficiency_at_full_load",
        ],
        weight_col="hours_in_service",
    )
    return weighted_so2_control_efficiency


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
        (cems["co2_mass_lb"] == 0) & (cems["fuel_consumed_mmbtu"] > 0),
        "co2_mass_lb",
    ] = np.NaN

    # check that all records with missing co2 have a non-missing energy source code
    missing_esc = cems[
        (cems["energy_source_code"].isna()) & (cems["co2_mass_lb"].isna())
    ]
    if len(missing_esc) > 0:
        print(" ")
        print(
            "Warning: the following units are missing co2 data and energy source codes. This may be because they burn multiple fuels."
        )
        print(missing_esc[["plant_id_eia", "unitid"]].drop_duplicates())
        print(" ")

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
        emission_factors, how="left", on="energy_source_code", validate="m:1"
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
        missing_gf, how="left", on=["plant_id_eia", "report_date"], validate="m:1"
    ).set_index(missing_index)

    # calculate missing co2 data
    missing_co2["co2_mass_lb"] = (
        missing_co2["fuel_consumed_mmbtu"] * missing_co2["weighted_ef"]
    )

    # update in CEMS table
    cems.update(missing_co2[["co2_mass_lb"]])

    # update the co2 mass measurement code
    cems.loc[missing_co2.index, "co2_mass_measurement_code"] = "Imputed"

    # check that there are no missing co2 values left
    still_missing_co2 = cems[cems["co2_mass_lb"].isna()]
    if len(still_missing_co2) > 0:
        raise UserWarning(
            "There are still misssing CO2 values remaining after filling missing CO2 values in CEMS"
        )

    return cems
