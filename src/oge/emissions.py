import gc
import pandas as pd
import numpy as np

import oge.load_data as load_data
import oge.validation as validation
from oge.column_checks import apply_dtypes, get_dtypes
from oge.filepaths import reference_table_folder
from oge.helpers import create_plant_ba_table
from oge.logging_util import get_logger
from oge.constants import (
    BIOMASS_FUELS,
    CLEAN_FUELS,
    ConversionFactors,
    chp_gross_thermal_output_efficiency,
    chp_useful_thermal_output_efficiency,
    nox_lb_per_mmbtu_flared_landfill_gas,
    current_early_release_year,
)

logger = get_logger(__name__)


def calculate_ghg_emissions_from_fuel_consumption(
    df: pd.DataFrame,
    year: int,
    include_co2: bool = True,
    include_ch4: bool = True,
    include_n2o: bool = True,
) -> pd.DataFrame:
    """Calculates GHG emissions from fuel consumption.

    Args:
        df (pd.DataFrame): data frame to operate on.
        year (int): a four-digit year.
        include_co2 (bool, optional): whether to include CO2 or not. Defaults to True.
        include_ch4 (bool, optional): whether to include CH4 or not. Defaults to True.
        include_n2o (bool, optional): whether to include N2O or not. Defaults to True.

    Returns:
        pd.DataFrame: original data frame with additional emission columns.
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

    # drop intermediate columns in-place to avoid memory allocation issues
    df.drop(columns=efs_to_use, inplace=True)

    return df


def add_geothermal_emission_factors(
    df: pd.DataFrame,
    year: int,
    include_co2: bool = True,
    include_nox: bool = True,
    include_so2: bool = True,
) -> pd.DataFrame:
    """Adds geothermal emission factors to a data frame.

    Args:
        df (pd.DataFrame): data frame to operate on.
        year (int): a four-digit year.
        include_co2 (bool, optional): whether to include CO2 or not. Defaults to True.
        include_nox (bool, optional): whether to include NOx or not. Defaults to True.
        include_so2 (bool, optional): whether to include SO2 or not. Defaults to True.

    Returns:
        pd.DataFrame: original data frame with additional geothermal emission columns.
    """
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

    # if there is a merge key for generator id, merge in the geothermal EFs on
    # generator id
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

    # drop intermediate columns in-place to avoid memory allocation issues
    cols_to_drop = [f"{e}_lb_per_mmbtu_geo" for e in emissions_to_calc]
    df.drop(columns=cols_to_drop, inplace=True, errors="ignore")

    return df


def calculate_geothermal_emission_factors(year: int) -> pd.DataFrame:
    """Updates the list of geothermal plants provided by EPA using EIA data. Calculates
    a weighted average EF for each plant-month based on the fraction of fuel consumed
    from each type of prime mover (steam, binary, flash).

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table enclosing the geothermal emission factors
    """
    # load geothermal efs
    geothermal_efs = pd.read_csv(
        reference_table_folder("geothermal_emission_factors.csv"),
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


def identify_geothermal_generator_geotype(year: int) -> pd.DataFrame:
    """Identifies whether each geothermal generator is binary, flash, or dry steam.


    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table indicating the type of each geothermal generator.
    """
    geothermal_geotype = load_data.load_pudl_table("core_eia860__scd_generators", year)
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

    # default steam turbines to flash steam b/c flash is more common than steam
    # according to EIA.
    # Source: https://www.eia.gov/energyexplained/geothermal/geothermal-power-plants.php
    geo_map = {"BT": "B", "ST": "F"}
    geothermal_geotype["geotype_code"] = geothermal_geotype["prime_mover_code"].map(
        geo_map
    )

    # According to the NREL Geothermal report:
    #  (https://www.nrel.gov/docs/fy21osti/78291.pdf)
    # The Geysers Complex in California contains the only "Dry Steam" (geotype S)
    # plants in the country. The Geysers has utility_id_eia = 7160
    geothermal_geotype.loc[
        (geothermal_geotype["utility_id_eia"] == 7160)
        & (geothermal_geotype["geotype_code"] == "F"),
        "geotype_code",
    ] = "S"

    # calculate what fraction of the plant's nameplate capcity each generator is
    # responsible for assume that missing power factor is 100%
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


def adjust_fuel_and_emissions_for_chp(df: pd.DataFrame) -> pd.DataFrame:
    """Allocates total emissions for electricity generation.

    Args:
        df (pd.DataFrame): dataframe containing fuel and emissions columns mapped to
            'subplant_id'. Must contain fuel_consumed_mmbtu, 'net_generation_mwh',
            and 'fuel_consumed_mmbtu_for_electricity' columns

    Returns:
        pd.DataFrame: the input dataframe with '_for_electricity' columns added for all
            fuel and emissions columns
    """

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

    # Drop in-place to avoid creating a copy of the large dataframe
    df.drop(columns=["electric_allocation_factor"], inplace=True)

    validation.test_chp_allocation(df)

    return df


def calculate_electric_allocation_factor(df: pd.DataFrame) -> pd.DataFrame:
    """Calculates an electric allocation factor for CHP subplants.

    The methodology is based on the formula in the eGRID2020 technical guide.

    We calculate the electric allocation factor for each subplant, rather than for each
    generator_id or emissions_unit_id_epa since fuel consumption and generation data are
    not always assigned to the same generator within a subplant.

    Calculates a unique EAF for each subplant-month-fuel.


    Args:
        df (pd.DataFrame): dataframe with 'fuel_consumed_mmbtu', 'net_generation_mwh',
            and fuel_consumed_mmbtu_for_electricity columns, mapped to 'subplant_id'

    Returns:
        pd.DataFrame: the input dataframe with an 'electric_allocation_factor' column
            added
    """
    groupby_cols = ["plant_id_eia", "subplant_id", "report_date"]
    # for EIA-923 we want to group on energy_source_code as well, but this does not
    # exist in CEMS
    if "energy_source_code" in df.columns:
        groupby_cols += ["energy_source_code"]

    df_subplant = (
        df.groupby(groupby_cols)[
            [
                "fuel_consumed_mmbtu",
                "fuel_consumed_for_electricity_mmbtu",
                "net_generation_mwh",
            ]
        ]
        .sum()
        .reset_index()
    )

    # calculate the gross thermal output
    # 0.8 is an assumed efficiency factor used by eGRID
    df_subplant["gross_thermal_output_for_heating_mmbtu"] = (
        chp_gross_thermal_output_efficiency
        * (
            df_subplant["fuel_consumed_mmbtu"]
            - df_subplant["fuel_consumed_for_electricity_mmbtu"]
        )
    )
    # calculate the useful thermal output
    # 0.75 is an assumed efficiency factor used by eGRID
    df_subplant["useful_thermal_output_mmbtu"] = (
        chp_useful_thermal_output_efficiency
        * df_subplant["gross_thermal_output_for_heating_mmbtu"]
    )

    # convert generation to mmbtu
    df_subplant["generation_mmbtu"] = (
        df_subplant["net_generation_mwh"] * ConversionFactors.mwh_to_mmbtu
    )

    # calculate the electric allocation factor
    df_subplant["electric_allocation_factor"] = df_subplant["generation_mmbtu"] / (
        df_subplant["generation_mmbtu"] + df_subplant["useful_thermal_output_mmbtu"]
    )

    # if the allocation factor < 0, set to zero
    # if the allocation factor > 1, set to one
    df_subplant["electric_allocation_factor"] = df_subplant[
        "electric_allocation_factor"
    ].clip(0, 1)

    # fill any missing factors with 1
    df_subplant["electric_allocation_factor"] = (
        df_subplant["electric_allocation_factor"].fillna(1).astype("float32")
    )

    # merge the allocation factors back into the original dataframe
    df = df.merge(
        df_subplant[groupby_cols + ["electric_allocation_factor"]],
        how="left",
        on=groupby_cols,
        validate="m:1",
    )

    return df


def adjust_emissions_for_biomass(df: pd.DataFrame) -> pd.DataFrame:
    """Adds emission in data frame where biomass emissions have been removed.

    Args:
        df (pd.DataFrame): original data frame.

    Returns:
        pd.DataFrame: original data frame with new adjusted emissions column that
            sets any biomass emissions to zero.
    """

    # create a column for adjusted biomass emissions, setting these emissions to zero
    # CO2: adjust emissions for co2 for all biomass generators
    if "co2_mass_lb" in df.columns:
        df["co2_mass_lb_adjusted"] = df["co2_mass_lb"]
        df.loc[df["energy_source_code"].isin(BIOMASS_FUELS), "co2_mass_lb_adjusted"] = 0

    # CH4: for landfill gas (LFG), all other emissions are set to zero
    # this assumes that the gas would have been flared anyway if not used for
    # electricity generation
    if "ch4_mass_lb" in df.columns:
        df["ch4_mass_lb_adjusted"] = df["ch4_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "ch4_mass_lb_adjusted"] = 0

    # N2O: LFG plants set to zero
    if "n2o_mass_lb" in df.columns:
        df["n2o_mass_lb_adjusted"] = df["n2o_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "n2o_mass_lb_adjusted"] = 0

    # NOx: assigned an adjusted value
    # this value is based on using NOx emissions from flaring as a baseline, and
    # subtracting this from the actual emissions to prevent negative emissions, we set
    # the value = 0 if negative
    if "nox_mass_lb" in df.columns:
        df["nox_mass_lb_adjusted"] = df["nox_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "nox_mass_lb_adjusted"] = df.loc[
            df["energy_source_code"] == "LFG", "nox_mass_lb_adjusted"
        ] - (
            df.loc[df["energy_source_code"] == "LFG", "fuel_consumed_mmbtu"]
            * nox_lb_per_mmbtu_flared_landfill_gas
        )
        df.loc[df["nox_mass_lb_adjusted"] < 0, "nox_mass_lb_adjusted"] = 0

    # SO2: LFG plants set to zero
    if "so2_mass_lb" in df.columns:
        df["so2_mass_lb_adjusted"] = df["so2_mass_lb"]
        df.loc[df["energy_source_code"] == "LFG", "so2_mass_lb_adjusted"] = 0

    return df


def calculate_co2e_mass(
    df: pd.DataFrame,
    version: int | str,
    gwp_horizon: int = 100,
    ar5_climate_carbon_feedback: bool = True,
) -> pd.DataFrame:
    """Calculate CO2-equivalent emissions from CO2, CH4, and N2O. This is done by
    choosing one of the IPCC's emission factors for CH4 and N2O.

    Note that the `version` parameter controls the version of the IPCC AR that will be
    chosen. If `version` is an int then it will be interpreted as a year and the most
    recent AR that was available that year will be picked. If `version` is a str then
    it will be interpreted as an IPCC Assessment Report (AR). Therefore, valid
    `version` are: int >= 1995 and any str equals to SAR, TAR, AR4, AR5, AR5_cc, or AR6.
    Note that when `version` is a str the `ar5_climate_carbon_feedback` has no effect
    since the user has direct access to AR5 or AR5_cc via the `version` parameter.


    If the `fuel_consumed_for_electricity_units` column is available, we also compute
    the adjusted emissions.

    Args:
        df (pd.DataFrame): data frame with CO2, CH4 and N2O emissions.
        version (int | str): version of the IPCC to consider.
        gwp_horizon (int, optional): the GWP horizon to consider. Defaults to 100.
        ar5_climate_carbon_feedback (bool, optional): When vversion is is a year, this
            parameter allows to select the climate carbon feedback of the fifth version
            of the assessment report. Defaults to True.

    Raises:
        ValueError: if `version` is invalid.

    Returns:
        pd.DataFrame: original data frame with the CO2e columns.
    """
    df_gwp = load_data.load_ipcc_gwp()

    if (
        isinstance(version, int) and version >= df_gwp["year"].min()
    ) or version == "latest":
        year = version if isinstance(version, int) else df_gwp["year"].max()

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
    elif version in df_gwp["ipcc_version"].unique():
        ipcc_version = version
    else:
        raise ValueError(
            f"version must be a year >= {df_gwp['year'].min()} or a valid AR"
        )

    gwp_to_use = df_gwp[df_gwp.ipcc_version == ipcc_version]

    ch4_gwp = gwp_to_use.loc[
        (gwp_to_use.gwp_horizon == gwp_horizon) & (gwp_to_use.gas == "ch4"),
        "gwp",
    ].item()
    n2o_gwp = gwp_to_use.loc[
        (gwp_to_use.gwp_horizon == gwp_horizon) & (gwp_to_use.gas == "n2o"),
        "gwp",
    ].item()

    # fill missing ch4 and n2o with zero so that the calculation works for geothermal
    # plants don't fill co2 so that if the data actually are misisng, a missing value
    # is also returned
    df["co2e_mass_lb"] = (
        df["co2_mass_lb"]
        + (ch4_gwp * df["ch4_mass_lb"].fillna(0))
        + (n2o_gwp * df["n2o_mass_lb"].fillna(0))
    )

    if "co2_mass_lb_adjusted" in df:
        df["co2e_mass_lb_adjusted"] = (
            df["co2_mass_lb_adjusted"]
            + (ch4_gwp * df["ch4_mass_lb_adjusted"].fillna(0))
            + (n2o_gwp * df["n2o_mass_lb_adjusted"].fillna(0))
        )
    if "co2_mass_lb_for_electricity" in df:
        df["co2e_mass_lb_for_electricity"] = (
            df["co2_mass_lb_for_electricity"]
            + (ch4_gwp * df["ch4_mass_lb_for_electricity"].fillna(0))
            + (n2o_gwp * df["n2o_mass_lb_for_electricity"].fillna(0))
        )
    if "co2_mass_lb_for_electricity_adjusted" in df:
        df["co2e_mass_lb_for_electricity_adjusted"] = (
            df["co2_mass_lb_for_electricity_adjusted"]
            + (ch4_gwp * df["ch4_mass_lb_for_electricity_adjusted"].fillna(0))
            + (n2o_gwp * df["n2o_mass_lb_for_electricity_adjusted"].fillna(0))
        )

    return df


def calculate_nox_from_fuel_consumption(
    gen_fuel_allocated: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Calculates NOx emissions from fuel consumption data.

    Note:
        Apply emission rates in the following order based on availability:
        1. Month- and generator-specific uncontrolled NOx emission factors based on
        boiler design parameters of all associated boilers and heat content of consumed
        fuels
        2. Season- and generator-specific controlled NOx emissions factors based on
        boiler-specific NOx control equipment in operation.

    Args:
        gen_fuel_allocated (pd.DataFrame): table enclosing generation and fuel
            consumption at generator level.
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: original data frame with NOx emissions.
    """

    # calculate uncontrolled NOx emission factors based on boiler design parameters
    uncontrolled_nox_factors = calculate_generator_nox_ef_per_unit_from_boiler_type(
        gen_fuel_allocated, year
    )
    # convert all EFs to a standardized unit based on fuel heat content
    uncontrolled_nox_factors = convert_ef_to_lb_per_mmbtu(
        uncontrolled_nox_factors, year, "nox"
    )

    # For generators associated with multiple boilers, average all NOx factors together
    # to have a single factor per generator
    # TODO: In the future, these factors should be weighed by the boiler fuel allocation
    # factors, and not just a straight average.
    uncontrolled_nox_factors = (
        uncontrolled_nox_factors.groupby(
            [
                "report_date",
                "plant_id_eia",
                "generator_id",
                "prime_mover_code",
                "energy_source_code",
            ],
            dropna=False,
        )["nox_ef_lb_per_mmbtu"]
        .mean()
        .reset_index()
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
    # raise a warning if we are missing emission factors for any non-zero fuel
    # consumption from non-clean fuels
    missing_ef = gen_fuel_allocated[
        gen_fuel_allocated["nox_ef_lb_per_mmbtu"].isna()
        & (gen_fuel_allocated["fuel_consumed_mmbtu"] > 0)
        & ~gen_fuel_allocated["energy_source_code"].isin(CLEAN_FUELS)
    ]

    if len(missing_ef) > 0:
        logger.warning("NOx emission factors are missing for the following records")
        logger.warning("Missing factors for FC prime movers are currently expected")
        logger.warning(
            "\n"
            + missing_ef[
                [
                    "report_date",
                    "plant_id_eia",
                    "energy_source_code",
                    "prime_mover_code",
                    "generator_id",
                ]
            ]
            .drop_duplicates()
            .merge(
                create_plant_ba_table(year)[["plant_id_eia", "ba_code"]],
                how="left",
                on="plant_id_eia",
                validate="m:1",
            )
            .to_string()
        )
    gen_fuel_allocated["nox_mass_lb"] = (
        gen_fuel_allocated["fuel_consumed_mmbtu"]
        * gen_fuel_allocated["nox_ef_lb_per_mmbtu"]
    )

    # calculate the controlled NOx rates
    controlled_nox_factors = calculate_generator_specific_controlled_nox_rates(year)

    # merge the controlled NOx rates into gen_fuel_allocated
    gen_fuel_allocated = gen_fuel_allocated.merge(
        controlled_nox_factors,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:1",
    )
    # calculate the controlled NOx emissions based on the month
    gen_fuel_allocated = gen_fuel_allocated.assign(
        controlled_nox_mass_lb=lambda x: np.where(
            ((x.report_date.dt.month >= 5) & (x.report_date.dt.month <= 9)),
            x.fuel_consumed_mmbtu * x.controlled_ozone_season_nox_ef_lb_per_mmbtu,
            x.fuel_consumed_mmbtu * x.controlled_non_ozone_season_nox_ef_lb_per_mmbtu,
        )
    )
    # if there were not season-specific rates, fill in using the annual rate if
    # available
    gen_fuel_allocated["controlled_nox_mass_lb"] = (
        gen_fuel_allocated["controlled_nox_mass_lb"]
        .astype(get_dtypes()["controlled_nox_mass_lb"])
        .fillna(
            gen_fuel_allocated["fuel_consumed_mmbtu"]
            * gen_fuel_allocated["controlled_annual_nox_ef_lb_per_mmbtu"]
        )
    )
    # update the emision total using the controlled mass if available
    gen_fuel_allocated.update(
        {"nox_mass_lb": gen_fuel_allocated["controlled_nox_mass_lb"]}
    )

    # remove intermediate columns in-place to avoid memory allocation issues
    gen_fuel_allocated.drop(
        columns=[
            "nox_ef_lb_per_mmbtu",
            "controlled_annual_nox_ef_lb_per_mmbtu",
            "controlled_ozone_season_nox_ef_lb_per_mmbtu",
            "controlled_non_ozone_season_nox_ef_lb_per_mmbtu",
            "controlled_nox_mass_lb",
        ],
        inplace=True,
        errors="ignore",
    )

    return gen_fuel_allocated


def calculate_generator_nox_ef_per_unit_from_boiler_type(
    gen_fuel_allocated: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Calculates a boiler-specific NOx emission factor per unit fuel based on boiler
    firing type, and associates with each generator.

    Args:
        gen_fuel_allocated (pd.DataFrame): table enclosing generation and fuel
            consumption at generator level.
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: original data frame with additional information regarding the
            generator's boiler type design and the associated uncotrolled NOx emission
            factor.
    """

    # get a dataframe with all unique generator-pm-esc combinations for emitting energy
    # source types with data reported
    gen_keys_for_nox = gen_fuel_allocated[
        (gen_fuel_allocated["fuel_consumed_mmbtu"] > 0)
    ].copy()
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
            "wet_dry_bottom",
            "boiler_firing_type",
        ]
    )

    # load the boiler firing type info
    boiler_firing_type = load_boiler_firing_type(year)

    # identify the boiler firing type for each generator
    boiler_generator_assn = load_data.load_pudl_table(
        "core_eia860__assn_boiler_generator",
        year,
        columns=["plant_id_eia", "boiler_id", "generator_id"],
    )
    # associate each boiler record with generator_id s
    boiler_firing_type = boiler_firing_type.merge(
        boiler_generator_assn,
        how="left",
        on=["plant_id_eia", "boiler_id"],
        validate="1:m",
    )

    # merge the gen keys with the boiler firing types
    gen_nox_factors = gen_keys_for_nox.merge(
        boiler_firing_type,
        how="left",
        on=["plant_id_eia", "generator_id"],
        validate="m:m",
    )

    gen_nox_factors["wet_dry_bottom"] = gen_nox_factors["wet_dry_bottom"].fillna("none")
    gen_nox_factors["boiler_firing_type"] = gen_nox_factors[
        "boiler_firing_type"
    ].fillna("none")

    # merge in the emission factors for specific boiler types
    gen_nox_factors = gen_nox_factors.merge(
        nox_emission_factors,
        how="left",
        on=[
            "prime_mover_code",
            "energy_source_code",
            "wet_dry_bottom",
            "boiler_firing_type",
        ],
        validate="m:1",
    )

    # fill in geotype-specific geothermal emission factors
    if gen_nox_factors["energy_source_code"].str.contains("GEO").any():
        gen_nox_factors = add_geothermal_emission_factors(
            gen_nox_factors,
            year,
            include_co2=False,
            include_nox=True,
            include_so2=False,
        )
        gen_nox_factors.loc[
            gen_nox_factors["energy_source_code"] == "GEO", "emission_factor"
        ] = gen_nox_factors.loc[
            gen_nox_factors["energy_source_code"] == "GEO", "nox_lb_per_mmbtu"
        ]
        gen_nox_factors.loc[
            gen_nox_factors["energy_source_code"] == "GEO", "emission_factor_numerator"
        ] = "lb"
        gen_nox_factors.loc[
            gen_nox_factors["energy_source_code"] == "GEO",
            "emission_factor_denominator",
        ] = "mmbtu"
        gen_nox_factors = gen_nox_factors.drop(columns="nox_lb_per_mmbtu")

    # identify missing emission factors and replace with PM-fuel factors if available
    missing_nox_efs = (
        gen_nox_factors.loc[
            gen_nox_factors["emission_factor"].isna(),
            [
                "prime_mover_code",
                "energy_source_code",
                "wet_dry_bottom",
                "boiler_firing_type",
            ],
        ]
        .drop_duplicates()
        .sort_values(
            by=[
                "energy_source_code",
                "prime_mover_code",
                "boiler_firing_type",
                "wet_dry_bottom",
            ]
        )
    )
    if len(missing_nox_efs) > 0:
        logger.warning(
            "NOx emission factors are missing for the following boiler types. "
            "A prime mover-fuel level factor will be used if available."
        )
        logger.warning("Missing factors for FC prime movers are currently expected")
        logger.warning("\n" + missing_nox_efs.to_string())
    gen_nox_factors = fill_missing_factors_based_on_pm_fuel(
        nox_emission_factors, gen_nox_factors
    )

    # identify missing emission factors and replace with zeros
    missing_nox_efs = (
        gen_nox_factors.loc[
            gen_nox_factors["emission_factor"].isna(),
            [
                "prime_mover_code",
                "energy_source_code",
                "wet_dry_bottom",
                "boiler_firing_type",
            ],
        ]
        .drop_duplicates()
        .sort_values(
            by=[
                "energy_source_code",
                "prime_mover_code",
                "boiler_firing_type",
                "wet_dry_bottom",
            ]
        )
    )

    if len(missing_nox_efs) > 0:
        logger.warning(
            """
            After filling with PM-fuel factors, NOx emission factors are still missing 
            for the following boiler types. An emission factor of zero will be used for 
            these boilers. Missing factors for FC prime movers are currently expected.
            """
        )
        logger.warning("\n" + missing_nox_efs.to_string())

    gen_nox_factors["emission_factor"] = gen_nox_factors["emission_factor"].fillna(0)

    return gen_nox_factors


def load_boiler_firing_type(year: int) -> pd.DataFrame:
    """Loads boiler firing type.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table enclosing the boiler design parameters.
    """
    boiler_firing_type = load_data.load_pudl_table(
        "core_eia860__scd_boilers",
        year,
        columns=[
            "plant_id_eia",
            "boiler_id",
            "firing_type_1",
            "wet_dry_bottom",
        ],
    )

    firing_types_eia = load_data.load_pudl_table("core_eia__codes_firing_types")

    boiler_firing_type["boiler_firing_type"] = (
        boiler_firing_type["firing_type_1"]
        .map(dict(zip(firing_types_eia["code"], firing_types_eia["label"])))
        .fillna("none")
        .str.lower()
    )

    boiler_firing_type["wet_dry_bottom"] = (
        boiler_firing_type["wet_dry_bottom"]
        .replace({"D": "dry", "W": "wet"})
        .fillna("none")
        .str.lower()
    )

    boiler_firing_type = boiler_firing_type[
        ["plant_id_eia", "boiler_id", "wet_dry_bottom", "boiler_firing_type"]
    ].dropna(subset=["wet_dry_bottom", "boiler_firing_type"], thresh=1)

    return boiler_firing_type


def fill_missing_factors_based_on_pm_fuel(
    emission_factors: pd.DataFrame, gen_factors: pd.DataFrame
) -> pd.DataFrame:
    """Fills emission factors based on the PM-fuel minimum factor if boiler firing type
    information is not available. The minimum factor for a PM-fuel grouping is used to
    be conservative, and is consistent with the eGRID methodology,

    Args:
        emission_factors (pd.DataFrame): data frame containing boiler-firing type
            specific NOx or SO2 emission factors.
        gen_factors (pd.DataFrame): generation and fuel consumption data frame with
            boiler-specific factors.

    Returns:
        pd.DataFrame: original data frame with missing factors filled out with
            available PM-fuel.
    """
    # identify the most conservative  emission factor for each PM-fuel combo
    pm_fuel_factors = (
        emission_factors[emission_factors["emission_factor_denominator"] != "mmbtu"]
        .groupby(
            [
                "energy_source_code",
                "prime_mover_code",
                "emission_factor_numerator",
                "emission_factor_denominator",
            ]
        )["emission_factor"]
        .min()
        .reset_index()
    )

    # merge the pm_fuel factors in and use them to fill any missing boiler-specific
    # factors
    gen_factors = gen_factors.merge(
        pm_fuel_factors,
        how="left",
        on=["energy_source_code", "prime_mover_code"],
        validate="m:1",
        suffixes=(None, "_pm_fuel"),
    )
    for col in [
        "emission_factor",
        "emission_factor_numerator",
        "emission_factor_denominator",
    ]:
        gen_factors[col] = gen_factors[col].fillna(gen_factors[f"{col}_pm_fuel"])
        gen_factors = gen_factors.drop(columns=f"{col}_pm_fuel")

    return gen_factors


def convert_ef_to_lb_per_mmbtu(gen_emission_factors, year, pollutant):
    # get the reported fuel heat content values from EIA-923
    (
        plant_specific_fuel_heat_content,
        national_avg_fuel_heat_content,
        annual_avg_fuel_heat_content,
    ) = return_monthly_plant_fuel_heat_content(year)

    # merge in plant pm specific fuel heat content values
    gen_emission_factors = gen_emission_factors.merge(
        plant_specific_fuel_heat_content,
        how="left",
        on=["report_date", "plant_id_eia", "energy_source_code", "prime_mover_code"],
        validate="m:1",
    )
    # merge in national monthly average fuel heat content values and use to fill
    # missing values
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
        logger.warning(
            "The heat content for the following fuels is missing and NOx emissions "
            "will not be calculated for these fuel: "
            f"{list(missing_fuel_content.energy_source_code.unique())}"
        )

    # convert emission factors from lb per unit to lb per mmbtu if the factor is not
    # already in units of lb/mmbtu
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


def return_monthly_plant_fuel_heat_content(
    year: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculates the heat content of fuels at different spatial/temporal resolution.

    Args:
        year (int): a four-digit year.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: the plant/monthly heat content
            of fuels, the national/monthly average heat content of fuels and the
            national/annual average of fuels.
    """
    # load information about the monthly heat input of fuels
    plant_specific_fuel_heat_content = load_data.load_pudl_table(
        "out_eia923__generation_fuel_combined",
        year,
        columns=[
            "plant_id_eia",
            "energy_source_code",
            "prime_mover_code",
            "report_date",
            "fuel_mmbtu_per_unit",
        ],
    )
    plant_specific_fuel_heat_content = plant_specific_fuel_heat_content[
        ~plant_specific_fuel_heat_content["energy_source_code"].isin(CLEAN_FUELS)
    ]

    # replace zero heat content with missing values
    plant_specific_fuel_heat_content["fuel_mmbtu_per_unit"] = (
        plant_specific_fuel_heat_content["fuel_mmbtu_per_unit"].replace(0, np.NaN)
    )
    plant_specific_fuel_heat_content["fuel_mmbtu_per_unit"] = (
        plant_specific_fuel_heat_content["fuel_mmbtu_per_unit"].replace(np.inf, np.NaN)
    )

    # calculate the average monthly heat content for a fuel
    national_avg_fuel_heat_content = (
        plant_specific_fuel_heat_content.drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "report_date"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    annual_avg_fuel_heat_content = (
        national_avg_fuel_heat_content.groupby(["energy_source_code"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    return (
        apply_dtypes(plant_specific_fuel_heat_content),
        apply_dtypes(national_avg_fuel_heat_content),
        apply_dtypes(annual_avg_fuel_heat_content),
    )


def calculate_unit_specific_controlled_nox_rates(year: int) -> pd.DataFrame:
    """Calculates controlled NOx emission rates at unit level.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table enclosing controlled NOx emission rates at unit level.
    """
    nox_rates = load_controlled_nox_emission_rates(year)
    nox_rates = calculate_non_ozone_season_nox_rate(nox_rates)
    weighted_nox_rates = calculate_weighted_nox_rates(
        year, nox_rates, "emissions_unit_id_epa"
    )

    return weighted_nox_rates


def calculate_boiler_specific_controlled_nox_rates(year: int) -> pd.DataFrame:
    """Calculates controlled NOx emission rates at boiler level.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table enclosing controlled NOx emission rates at boiler level.
    """
    nox_rates = load_controlled_nox_emission_rates(year)
    nox_rates = calculate_non_ozone_season_nox_rate(nox_rates)
    weighted_nox_rates = calculate_weighted_nox_rates(year, nox_rates, "boiler_id")

    return weighted_nox_rates


def calculate_generator_specific_controlled_nox_rates(year: int) -> pd.DataFrame:
    """Calculates controlled NOx emission rates at generator level.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table enclosing controlled NOx emission rates at generator level.
    """
    nox_rates = load_controlled_nox_emission_rates(year)
    nox_rates = calculate_non_ozone_season_nox_rate(nox_rates)
    weighted_nox_rates = calculate_weighted_nox_rates(year, nox_rates, "generator_id")

    return weighted_nox_rates


def load_controlled_nox_emission_rates(year: int) -> pd.DataFrame:
    """Loads contrilled NOx emission rates data.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table enclosing controlled annual NOx emission rates and the
            controlled NOx emission rates for the ozone season.
    """
    # load the emissions control data
    emissions_controls_eia923 = load_data.load_emissions_controls_eia923(year)

    # create a dataframe that contains only NOx emission data for operating control
    # equipment
    nox_rates = emissions_controls_eia923.dropna(
        axis="index",
        how="all",
        subset=[
            "annual_nox_emission_rate_lb_per_mmbtu",
            "ozone_season_nox_emission_rate_lb_per_mmbtu",
        ],
    )
    nox_rates = nox_rates[nox_rates["operational_status"] == "OP"]
    nox_rates = nox_rates[
        [
            "plant_id_eia",
            "nox_control_id_eia",
            "particulate_control_id_eia",
            "so2_control_id_eia",
            "mercury_control_id_eia",
            "hours_in_service",
            "annual_nox_emission_rate_lb_per_mmbtu",
            "ozone_season_nox_emission_rate_lb_per_mmbtu",
        ]
    ]
    nox_rates = nox_rates.rename(
        columns={
            "annual_nox_emission_rate_lb_per_mmbtu": "controlled_annual_nox_ef_lb_per_mmbtu",
            "ozone_season_nox_emission_rate_lb_per_mmbtu": "controlled_ozone_season_nox_ef_lb_per_mmbtu",
        }
    )

    return nox_rates


def calculate_weighted_nox_rates(
    year: int, nox_rates: pd.DataFrame, aggregation_level: str
) -> pd.DataFrame:
    """Aggregates NOx rates data for a given aggregation level.

    Args:
        year (int): a four-digit year
        nox_rates (pd.DataFrame): table enclosing the controlled NOx emission rates.
        aggregation_level (str): aggregation level to consider. Can be either
        'generator_id' or 'emissions_unit_id_epa'.

    Returns:
        pd.DataFrame: table enclosing NOx rates for a given aggregation level weighted
            by the number of hours the control equipment has been in service.
    """

    nox_rates = associate_control_ids_with_boiler_id(
        nox_rates,
        year,
        pollutant="nox",
    )

    if aggregation_level == "generator_id":
        boiler_generator_assn = load_data.load_pudl_table(
            "core_eia860__assn_boiler_generator",
            year,
            columns=["plant_id_eia", "boiler_id", "generator_id"],
        )
        # associate a generator_id with each record
        nox_rates = nox_rates.merge(
            boiler_generator_assn,
            how="left",
            on=["plant_id_eia", "boiler_id"],
            validate="m:m",
        )
    elif aggregation_level == "emissions_unit_id_epa":
        # load subplant crosswalk and map boiler-generator associations
        unit_boiler_assn = load_data.load_unit_to_boiler_associations(year)
        # associate a emissions_unit_id_epa with each record
        nox_rates = nox_rates.merge(
            unit_boiler_assn,
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


def associate_control_ids_with_boiler_id(
    df: pd.DataFrame, year: int, pollutant: str
) -> pd.DataFrame:
    """Maps emission control id with boiler id.

    Because some emission control data is missing a corresponding control id for some
    pollutant, this function attempts to use all control ids to associate each
    obervation with a boiler id, in a specific order.

    Args:
        df (pd.DataFrame): original data frame with boiler id.
        year (int): a four-digit year.
        pollutant (str): pollutant to consider. Either 'particulate', 'so2', 'nox' or
        'mercury'

    Returns:
        pd.DataFrame: original data frame with additional column indicating the
            control id for pollutant.
    """
    all_pollutants = ["particulate", "so2", "nox", "mercury"]
    # reorder the pollutant list to make the primary pollutant first
    all_pollutants.remove(pollutant)
    pollutant_order = [pollutant] + all_pollutants

    boiler_association_eia860 = load_data.load_pudl_table(
        "core_eia860__assn_yearly_boiler_emissions_control_equipment", year
    )

    counter = 1

    for pol in pollutant_order:
        pol_control_id_assn = boiler_association_eia860[
            boiler_association_eia860["emission_control_id_type"] == pol
        ].rename(columns={"emission_control_id_eia": f"{pol}_control_id_eia"})

        # if this is the first time through the loop
        if counter == 1:
            # load the association table and rename the boiler id column to specify
            # which table it came from

            df = df.merge(
                pol_control_id_assn[
                    [
                        "plant_id_eia",
                        f"{pol}_control_id_eia",
                        "boiler_id",
                    ]
                ],
                how="left",
                on=["plant_id_eia", f"{pol}_control_id_eia"],
                validate="m:m",
            )
            counter += 1
        else:
            # split out the data that is still missing a boiler_id
            missing_boiler_id = df[
                df["boiler_id"].isna() & ~df[f"{pol}_control_id_eia"].isna()
            ].drop(columns=["boiler_id"])
            # remove this data from the original dataframe
            df = df[~(df["boiler_id"].isna() & ~df[f"{pol}_control_id_eia"].isna())]

            missing_boiler_id = missing_boiler_id.merge(
                pol_control_id_assn[
                    [
                        "plant_id_eia",
                        f"{pol}_control_id_eia",
                        "boiler_id",
                    ]
                ],
                how="left",
                on=["plant_id_eia", f"{pol}_control_id_eia"],
                validate="m:m",
            )
            # add this data back to the original dataframe
            df = pd.concat([df, missing_boiler_id], axis=0, copy=False)

    # if there are any missing boiler_ids, fill using the nox_control_id, which is
    # likely to match a boiler
    df["boiler_id"] = df["boiler_id"].fillna(df[f"{pollutant}_control_id_eia"])

    return df


def calculate_weighted_averages(
    df: pd.DataFrame, groupby_columns: list[str], data_cols: list[str], weight_col: str
) -> pd.DataFrame:
    """Calculates weighted averages of one or more columns in a dataframe.

    Args:
        df (pd.DataFrame): original data frame.
        groupby_columns (list[str]): column(s) used to group the data frame before
            calculating the weighted average
        data_cols (list[str]): column(s) for which the weighted average will be
            calculated.
        weight_col (str): column to use as weight.

    Returns:
        pd.DataFrame: data frame with `groupby_columns` and weighted average column(s).
    """
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
    if result.index.empty:
        for col in groupby_columns:
            result[col] = ""
    else:
        result = result.reset_index()

    return result


def calculate_non_ozone_season_nox_rate(
    weighted_nox_rates: pd.DataFrame,
) -> pd.DataFrame:
    """Calculates NOx emission rates for non-ozone season.

    Args:
        weighted_nox_rates (pd.DataFrame): original data frame with the annual and
            ozone-season Nox emission rates.

    Returns:
        pd.DataFrame: original data frame with an additional non-ozone season NOx
            emission rates column.
    """
    annual_col = "controlled_annual_nox_ef_lb_per_mmbtu"
    oz_col = "controlled_ozone_season_nox_ef_lb_per_mmbtu"
    non_oz_col = "controlled_non_ozone_season_nox_ef_lb_per_mmbtu"

    # ozone season is May - Sept (5 months).
    # To get the average emission rate for the 7 non-ozone season months, we assume
    # similar operation across all months
    # annual_avg = [(5* oz_avg) + (7 * non_oz_avg)] / 12
    weighted_nox_rates[non_oz_col] = (
        (12 * weighted_nox_rates[annual_col]) - (5 * weighted_nox_rates[oz_col])
    ) / 7

    # if there are any rates that we calculate as negative, replace with zero
    weighted_nox_rates.loc[weighted_nox_rates[non_oz_col] < 0, non_oz_col] = 0

    return weighted_nox_rates


def calculate_so2_from_fuel_consumption(
    gen_fuel_allocated: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Calculates SO2 emissions from fuel consumption data.

    Note:
        Apply emission rates in the following order based on availability:
        1. Month- and generator-specific uncontrolled SO2 emission factors based on
        boiler design parameters of all associated boilers and heat content of consumed
        fuels
        2. Adjust SO2 emissions using boiler-specific SO2 removal efficiencies if
        available. Otherwise, assume SO2 emissions are uncontrolled.

    Args:
        gen_fuel_allocated (pd.DataFrame): table enclosing generation and fuel
            consumption at generator level.
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: original data frame with SO2 emissions
    """

    # load uncontrolled emission factors
    uncontrolled_so2_factors = calculate_generator_so2_ef_per_unit_from_boiler_type(
        gen_fuel_allocated, year
    )

    # adjust for sulfur content
    uncontrolled_so2_factors = adjust_so2_efs_for_fuel_sulfur_content(
        uncontrolled_so2_factors, year
    )

    # convert to mmbtu
    uncontrolled_so2_factors = convert_ef_to_lb_per_mmbtu(
        uncontrolled_so2_factors, year, "so2"
    )

    # average the emission factors for all boilers associated with each generator
    uncontrolled_so2_factors = (
        uncontrolled_so2_factors.groupby(
            [
                "report_date",
                "plant_id_eia",
                "generator_id",
                "prime_mover_code",
                "energy_source_code",
            ],
            dropna=False,
        )["so2_ef_lb_per_mmbtu"]
        .mean()
        .reset_index()
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
    # raise a warning if we are missing emission factors for any non-zero fuel
    # consumption from non-clean fuels
    missing_ef = gen_fuel_allocated[
        gen_fuel_allocated["so2_ef_lb_per_mmbtu"].isna()
        & (gen_fuel_allocated["fuel_consumed_mmbtu"] > 0)
        & ~gen_fuel_allocated["energy_source_code"].isin(CLEAN_FUELS)
    ]

    if len(missing_ef) > 0:
        logger.warning("SO2 emission factors are missing for the above records")
        logger.warning("Missing factors for FC prime movers are currently expected")
        logger.warning(
            "\n"
            + missing_ef[
                [
                    "report_date",
                    "plant_id_eia",
                    "energy_source_code",
                    "prime_mover_code",
                    "generator_id",
                ]
            ]
            .drop_duplicates()
            .merge(
                create_plant_ba_table(year)[["plant_id_eia", "ba_code"]],
                how="left",
                on="plant_id_eia",
                validate="m:1",
            )
            .to_string()
        )
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
    gen_fuel_allocated["so2_removal_efficiency_annual"] = (
        gen_fuel_allocated["so2_removal_efficiency_annual"]
        .astype(get_dtypes()["so2_removal_efficiency_annual"])
        .fillna(0)
    )

    # calculate controlled so2 emissions
    gen_fuel_allocated["so2_mass_lb"] = gen_fuel_allocated["so2_mass_lb"] * (
        1 - gen_fuel_allocated["so2_removal_efficiency_annual"]
    )

    # remove intermediate columns in-place to avoid memory allocation issues
    gen_fuel_allocated.drop(
        columns=[
            "so2_ef_lb_per_mmbtu",
            "so2_removal_efficiency_annual",
            "so2_removal_efficiency_at_full_load",
        ],
        inplace=True,
        errors="ignore",
    )

    return gen_fuel_allocated


def calculate_generator_so2_ef_per_unit_from_boiler_type(
    gen_fuel_allocated: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Calculates a generator-specific SO2 emission factor per unit fuel based on
    boiler firing type

    Note:
        If a generator has multiple boilers, average the emission factor of all boilers.

    Args:
        gen_fuel_allocated (pd.DataFrame): table enclosing generation and fuel
            consumption at generator level.
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: original data frame with additional information regarding the
            generator's boiler type design and the associated uncotrolled SO2 emission
            factor.
    """

    # get a dataframe with all unique generator-pm-esc combinations for emitting
    # energy source types with data reported
    gen_keys_for_so2 = gen_fuel_allocated[
        (gen_fuel_allocated["fuel_consumed_mmbtu"] > 0)
    ].copy()
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
    boiler_firing_type = boiler_firing_type.drop(columns="wet_dry_bottom")
    boiler_firing_type = boiler_firing_type.drop_duplicates()

    # identify the boiler firing type for each generator
    boiler_generator_assn = load_data.load_pudl_table(
        "core_eia860__assn_boiler_generator",
        year,
        columns=["plant_id_eia", "boiler_id", "generator_id"],
    )
    # associate each boiler record with generator_id s
    boiler_firing_type = boiler_firing_type.merge(
        boiler_generator_assn,
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

    gen_so2_factors["boiler_firing_type"] = gen_so2_factors[
        "boiler_firing_type"
    ].fillna("none")

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

    # identify missing emission factors and replace with PM-fuel specific factors if
    # available
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
        logger.warning(
            "SO2 emission factors are missing for the following boiler types. "
            "A prime mover-fuel level factor will be used if available."
        )
        logger.warning("Missing factors for FC prime movers are currently expected")
        logger.warning("\n" + missing_so2_efs.to_string())
    gen_so2_factors = fill_missing_factors_based_on_pm_fuel(
        so2_emission_factors, gen_so2_factors
    )

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
        logger.warning(
            "SO2 emission factors are missing for the following boiler types. "
            "An emission factor of zero will be used for these boilers."
        )
        logger.warning("Missing factors for FC prime movers are currently expected")
        logger.warning("\n" + missing_so2_efs.to_string())
    gen_so2_factors["emission_factor"] = gen_so2_factors["emission_factor"].fillna(0)
    gen_so2_factors["multiply_by_sulfur_content"] = gen_so2_factors[
        "multiply_by_sulfur_content"
    ].fillna(0)

    return gen_so2_factors


def return_monthly_plant_fuel_sulfur_content(year: int) -> pd.DataFrame:
    """Returns the monthly plant-level sulfur content percentage.

    Note:
        Sulfur content values are on a 0-100 scale (e.g. 5.2% = 5.2)

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table giving sulfur content percentage for each plant, energy
            source code and report date in `year`.
    """
    plant_specific_fuel_sulfur_content = load_plant_specific_fuel_sulfur_content(year)
    # drop unnecessary column
    plant_specific_fuel_sulfur_content = plant_specific_fuel_sulfur_content.drop(
        columns="state"
    )
    return apply_dtypes(plant_specific_fuel_sulfur_content)


def return_multiyear_plant_fuel_sulfur_content(years: list[int]) -> pd.DataFrame:
    """Returns a multiyear averaged plant-level sulfur content percentage.

    Note:
        Sulfur content values are on a 0-100 scale (e.g. 5.2% = 5.2)

    Args:
        years (list[int]): list of four-digit years.

    Returns:
        pd.DataFrame: table giving sulfur content percentage for each plant and energy
            source code averaged over `years`.
    """
    multiyear_avg_plant_specific_fuel_sulfur_content = (
        pd.concat(
            [load_plant_specific_fuel_sulfur_content(y) for y in years], copy=False
        )
        .groupby(["plant_id_eia", "energy_source_code"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    return apply_dtypes(multiyear_avg_plant_specific_fuel_sulfur_content)


def return_multiyear_state_fuel_sulfur_content(years: list[int]) -> pd.DataFrame:
    """Returns a multiyear averaged state-level sulfur content percentage.

    Note:
        Sulfur content values are on a 0-100 scale (e.g. 5.2% = 5.2)

    Args:
        years (list[int]): list of four-digit years.

    Returns:
        pd.DataFrame: table giving sulfur content percentage for each energy source
            code and state averaged over `years`.
    """
    multiyear_avg_state_fuel_sulfur_content = (
        pd.concat(
            [load_plant_specific_fuel_sulfur_content(y) for y in years], copy=False
        )
        .drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "state"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    return apply_dtypes(multiyear_avg_state_fuel_sulfur_content)


def return_monthly_state_fuel_sulfur_content(year: int) -> pd.DataFrame:
    """Returns the monthy state-level sulfur content percentage for each fuel.

    Note:
        Sulfur content values are on a 0-100 scale (e.g. 5.2% = 5.2)

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table giving sulfur content percentage for each energy source
            code and state.
    """
    plant_specific_fuel_sulfur_content = load_plant_specific_fuel_sulfur_content(year)
    state_avg_fuel_sulfur_content = (
        plant_specific_fuel_sulfur_content.drop(columns=["plant_id_eia"])
        .groupby(["energy_source_code", "report_date", "state"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    return apply_dtypes(state_avg_fuel_sulfur_content)


def return_monthly_national_fuel_sulfur_content(year: int) -> pd.DataFrame:
    """Returns the monthly national-level sulfur content percentage for each fuel.

    Note:
        Sulfur content values are on a 0-100 scale (e.g. 5.2% = 5.2)

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table giving sulfur content percentage for each energy source
            code and report date in `year`.
    """
    plant_specific_fuel_sulfur_content = load_plant_specific_fuel_sulfur_content(year)
    national_avg_fuel_sulfur_content = (
        plant_specific_fuel_sulfur_content.drop(columns=["plant_id_eia", "state"])
        .groupby(["energy_source_code", "report_date"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )
    return apply_dtypes(national_avg_fuel_sulfur_content)


def return_annual_national_fuel_sulfur_content(year: int) -> pd.DataFrame:
    """Returns the annual national-level sulfur content percentage for each fuel.

    Note:
        Sulfur content values are on a 0-100 scale (e.g. 5.2% = 5.2)

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table giving the annual sulfur content percentage for each energy
            source code.
    """
    national_avg_fuel_sulfur_content = return_monthly_national_fuel_sulfur_content(year)

    annual_avg_fuel_sulfur_content = (
        national_avg_fuel_sulfur_content.groupby(["energy_source_code"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    # if there are any missing annual average values, attempt to fill using data from
    # previous years
    if annual_avg_fuel_sulfur_content["sulfur_content_pct"].isna().any():
        previous_year_values = load_plant_specific_fuel_sulfur_content((year - 1))
        previous_year_values = (
            previous_year_values.groupby(["energy_source_code"], dropna=False)[
                ["sulfur_content_pct"]
            ]
            .mean(numeric_only=True)
            .reset_index()
        )
        annual_avg_fuel_sulfur_content = annual_avg_fuel_sulfur_content.merge(
            previous_year_values,
            how="left",
            on="energy_source_code",
            validate="1:1",
            suffixes=(None, "_fill"),
        )
        annual_avg_fuel_sulfur_content["sulfur_content_pct"] = (
            annual_avg_fuel_sulfur_content["sulfur_content_pct"].fillna(
                annual_avg_fuel_sulfur_content["sulfur_content_pct_fill"]
            )
        )
        annual_avg_fuel_sulfur_content = annual_avg_fuel_sulfur_content.drop(
            columns=["sulfur_content_pct_fill"]
        )

    return apply_dtypes(annual_avg_fuel_sulfur_content)


def return_multiyear_national_fuel_sulfur_content(years: list[int]) -> pd.DataFrame:
    """Returns a multiyear averaged national-level sulfur content percentage.

    Note:
        Sulfur content values are on a 0-100 scale (e.g. 5.2% = 5.2)

    Args:
        years (list[int]): list of four-digit years.

    Returns:
        pd.DataFrame: table giving sulfur content percentage for each energy source
            code averaged over `years`.
    """
    multiyear_avg_national_fuel_sulfur_content = (
        pd.concat(
            [load_plant_specific_fuel_sulfur_content(y) for y in years], copy=False
        )
        .drop(columns=["plant_id_eia", "state"])
        .groupby(["energy_source_code"], dropna=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    return apply_dtypes(multiyear_avg_national_fuel_sulfur_content)


def load_plant_specific_fuel_sulfur_content(year: int) -> pd.DataFrame:
    """Calculates the weighted average sulfur content of each fuel by the fuel
    consumption.

    Args:
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: the weighted average sulfur content percentage.
    """
    plant_specific_fuel_sulfur_content = load_data.load_pudl_table(
        "out_eia923__boiler_fuel",
        year,
        columns=[
            "plant_id_eia",
            "prime_mover_code",
            "energy_source_code",
            "report_date",
            "fuel_consumed_units",
            "sulfur_content_pct",
        ],
    )

    plant_specific_fuel_sulfur_content = plant_specific_fuel_sulfur_content[
        ~plant_specific_fuel_sulfur_content["energy_source_code"].isin(CLEAN_FUELS)
    ]

    # calculate a weighted average for each generator
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

    # Add state
    plant_state = load_data.load_pudl_table(
        "core_eia__entity_plants", columns=["plant_id_eia", "state"]
    )
    plant_specific_fuel_sulfur_content = plant_specific_fuel_sulfur_content.merge(
        plant_state, on="plant_id_eia", how="left", validate="m:1"
    )

    return plant_specific_fuel_sulfur_content


def adjust_so2_efs_for_fuel_sulfur_content(
    uncontrolled_so2_factors: pd.DataFrame, year: int
) -> pd.DataFrame:
    """Adjusts SO2 emission factors for fuel sulfur content.

    Args:
        uncontrolled_so2_factors (pd.DataFrame): table enclosing the uncontrolled SO2
            emission factors.
        year (int): a four-digit year.

    Returns:
        pd.DataFrame: table enclosing the uncotrolled SO2 emission factors adjusted for
            sulfur content
    """
    plant_state = load_data.load_pudl_table(
        "core_eia__entity_plants", columns=["plant_id_eia", "state"]
    )
    uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
        plant_state, on="plant_id_eia", how="left", validate="m:1"
    )

    if year >= 2008:
        # merge in monthly/plant level
        uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
            return_monthly_plant_fuel_sulfur_content(year),
            how="left",
            on=[
                "report_date",
                "plant_id_eia",
                "energy_source_code",
                "prime_mover_code",
            ],
            validate="m:1",
        )
        # merge in monthly/state level
        uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
            return_monthly_state_fuel_sulfur_content(year),
            how="left",
            on=["report_date", "energy_source_code", "state"],
            validate="m:1",
            suffixes=(None, "_state"),
        )
        # use to fill missing values
        uncontrolled_so2_factors["sulfur_content_pct"] = uncontrolled_so2_factors[
            "sulfur_content_pct"
        ].fillna(uncontrolled_so2_factors["sulfur_content_pct_state"])

        # annual/national level
        national_avg_fuel_sulfur_content = return_annual_national_fuel_sulfur_content(
            year
        )

        # merge in annual average fuel sulfur content and use to fill missing values
        uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
            national_avg_fuel_sulfur_content,
            how="left",
            on=[
                "energy_source_code",
            ],
            validate="m:1",
            suffixes=(None, "_national"),
        )
        uncontrolled_so2_factors["sulfur_content_pct"] = uncontrolled_so2_factors[
            "sulfur_content_pct"
        ].fillna(uncontrolled_so2_factors["sulfur_content_pct_national"])

        # as a final backstop, use national-average, multi-year
        # multiyear/national level from the past 5 years
        national_multi_avg_fuel_sulfur_content = (
            return_multiyear_national_fuel_sulfur_content(
                range(current_early_release_year - 4, current_early_release_year + 1)
            )
        )

        # merge in annual average fuel sulfur content and use to fill missing values
        uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
            national_multi_avg_fuel_sulfur_content,
            how="left",
            on=[
                "energy_source_code",
            ],
            validate="m:1",
            suffixes=(None, "_national_multi"),
        )
        uncontrolled_so2_factors["sulfur_content_pct"] = uncontrolled_so2_factors[
            "sulfur_content_pct"
        ].fillna(uncontrolled_so2_factors["sulfur_content_pct_national_multi"])

    else:
        logger.info(
            f"Fuel sulfur content is not available for {year}. Using 2008-2012 "
            "data to derive a plant, state and national fuel sulfur content average"
        )
        # merge in multiyear/plant level
        uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
            return_multiyear_plant_fuel_sulfur_content(range(2008, 2013)),
            how="left",
            on=[
                "plant_id_eia",
                "energy_source_code",
            ],
            validate="m:1",
        )
        # multiyear/state level
        uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
            return_multiyear_state_fuel_sulfur_content(range(2008, 2013)),
            how="left",
            on=["energy_source_code", "state"],
            validate="m:1",
            suffixes=(None, "_state"),
        )
        # use to fill missing values
        uncontrolled_so2_factors["sulfur_content_pct"] = uncontrolled_so2_factors[
            "sulfur_content_pct"
        ].fillna(uncontrolled_so2_factors["sulfur_content_pct_state"])

        # multiyear/national level
        national_avg_fuel_sulfur_content = (
            return_multiyear_national_fuel_sulfur_content(range(2008, 2013))
        )

        # merge in annual average fuel sulfur content and use to fill missing values
        uncontrolled_so2_factors = uncontrolled_so2_factors.merge(
            national_avg_fuel_sulfur_content,
            how="left",
            on=[
                "energy_source_code",
            ],
            validate="m:1",
            suffixes=(None, "_national"),
        )
        uncontrolled_so2_factors["sulfur_content_pct"] = uncontrolled_so2_factors[
            "sulfur_content_pct"
        ].fillna(uncontrolled_so2_factors["sulfur_content_pct_national"])

    # check that we are not missing sulfur content for any generators that need it
    missing_sulfur_content = uncontrolled_so2_factors[
        uncontrolled_so2_factors["sulfur_content_pct"].isna()
        & (uncontrolled_so2_factors["multiply_by_sulfur_content"] == 1)
    ]

    if len(missing_sulfur_content) > 0:
        logger.warning("Sulfur content data is missing in EIA-923 for the below units.")
        logger.warning(
            "\n"
            + missing_sulfur_content[
                [
                    "plant_id_eia",
                    "generator_id",
                    "prime_mover_code",
                    "energy_source_code",
                ]
            ]
            .drop_duplicates()
            .merge(
                create_plant_ba_table(year)[["plant_id_eia", "ba_code"]],
                how="left",
                on="plant_id_eia",
                validate="m:1",
            )
            .to_string()
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
            "sulfur_content_pct_state",
            "sulfur_content_pct_national",
            "state",
        ]
    )

    return uncontrolled_so2_factors


def load_so2_control_efficiencies(year: int) -> pd.DataFrame:
    """Loads SO2 control efficiencies data.

    Args:
        year (int): a four-digit year.

    Raises:
        UserWarning: if one or more annual SO2 remoaval efficiencies are not in [0, 1].

    Returns:
        pd.DataFrame: table enclosing SO2 emission data for operating control equipment.
    """
    # load the emissions control data
    emissions_controls_eia923 = load_data.load_emissions_controls_eia923(year)

    # remove control efficiencies for fluidized bed combustion since this is already
    # incorporated into the uncontrolled so2 rates
    emissions_controls_eia923.loc[
        emissions_controls_eia923["equipment_tech_description"].str.contains(
            "fluidized bed"
        ),
        ["so2_removal_efficiency_annual", "so2_removal_efficiency_at_full_load"],
    ] = np.NaN

    # remove control efficiencies for fluidized bed combustion since this is already
    # incorporated into the uncontrolled so2 rates
    emissions_controls_eia923.loc[
        emissions_controls_eia923["equipment_tech_description"].str.contains(
            "fluidized bed"
        ),
        ["so2_removal_efficiency_annual", "so2_removal_efficiency_at_full_load"],
    ] = np.NaN

    # create a dataframe that contains only SO2 emission data for operating control
    # equipment
    so2_efficiency = emissions_controls_eia923.dropna(
        axis="index",
        how="all",
        subset=[
            "so2_removal_efficiency_annual",
            "so2_removal_efficiency_at_full_load",
        ],
    )
    so2_efficiency = so2_efficiency[so2_efficiency["operational_status"] == "OP"]
    so2_efficiency = so2_efficiency[
        [
            "plant_id_eia",
            "so2_control_id_eia",
            "nox_control_id_eia",
            "particulate_control_id_eia",
            "mercury_control_id_eia",
            "hours_in_service",
            "so2_removal_efficiency_annual",
            "so2_removal_efficiency_at_full_load",
        ]
    ]

    # validate that the efficiencies are between 0 and 1
    bad_efficiencies = so2_efficiency[
        (so2_efficiency["so2_removal_efficiency_annual"] < 0)
        | (so2_efficiency["so2_removal_efficiency_annual"] > 1)
    ]
    if len(bad_efficiencies) > 0:
        raise UserWarning(
            "certain loaded SO2 removal efficiencies are either negative or > 100%"
        )

    return so2_efficiency


def calculate_weighted_so2_control_efficiency(
    year: int, so2_control_efficiency: pd.DataFrame, aggregation_level: str
) -> pd.DataFrame:
    """Aggregates SO2 rates data for a given aggregation level.

    Args:
        year (int): a four-digit year.
        so2_control_efficiency (pd.DataFrame): table enclosing the SO2 control
            efficiencies.
        aggregation_level (str): aggregation level.

    Returns:
        pd.DataFrame: table enclosing SO2 control efficiency at generator level weighted
            by the number of hours the control equipment has been in service.
    """
    so2_control_efficiency = associate_control_ids_with_boiler_id(
        so2_control_efficiency,
        year,
        pollutant="so2",
    )

    if aggregation_level == "generator_id":
        boiler_generator_assn = load_data.load_pudl_table(
            "core_eia860__assn_boiler_generator",
            year,
            columns=["plant_id_eia", "boiler_id", "generator_id"],
        )
        # associate a generator_id with each record
        so2_control_efficiency = so2_control_efficiency.merge(
            boiler_generator_assn,
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


def fill_cems_missing_co2(
    cems: pd.DataFrame, year: int, subplant_emission_factors: pd.DataFrame
) -> pd.DataFrame:
    """Fills missing hourly CO2 data in CEMS based on a three-tiered approach.

    Note:
        CO2 data is considered missing if reported CO2 is zero and fuel consumption
        is positive.
        1. If a unit has non-missing emissions data for other hours in the same month,
        calculate a unit-month specific emission factor from the CEMS-reported fuel
        consumption and emissions
        2. For all remaining missing values, use the subplant and month-specific
        weighted average emission factor from `subplant_emission_factors` calculated
        from the EIA-923 data
        3. For any remaining missing values, calculate emissions based on the subplant
        primary fuel and fuel consumption

    Args:
        cems (pd.DataFrame): CEMS data frame.
        year (int): a four-digit year.
        subplant_emission_factors (pd.DataFrame): subplant emission factors data frame.

    Raises:
        UserWarning: if there are missing CO2 values

    Returns:
        pd.DataFrame: CEMS data frame with filled CO2 data.
    """

    # make a copy of the cems data so that we can validate the outputs
    cems_original = cems.copy()
    # add a new categorical option to the mass measurement code
    cems["co2_mass_measurement_code"] = cems[
        "co2_mass_measurement_code"
    ].cat.add_categories("Imputed")

    # replace all "missing" CO2 values with zero
    cems["co2_mass_lb"] = cems["co2_mass_lb"].fillna(0)

    # replace 0 reported CO2 values with missing values, if there was reported heat
    # input
    cems.loc[
        (cems["co2_mass_lb"] == 0) & (cems["fuel_consumed_mmbtu"] > 0),
        "co2_mass_lb",
    ] = np.NaN

    # create a new df with all observations with missing co2 data
    missing_co2 = cems[cems["co2_mass_lb"].isnull()]

    # First round of filling covers small gaps using non-missing emission data from
    # the same month
    ####################################################################################

    unit_months_missing_co2 = missing_co2[
        ["plant_id_eia", "emissions_unit_id_epa", "report_date"]
    ].drop_duplicates()

    # get non-missing data from cems for these unit months
    unit_months_missing_co2 = unit_months_missing_co2.merge(
        cems[
            [
                "plant_id_eia",
                "emissions_unit_id_epa",
                "report_date",
                "co2_mass_lb",
                "fuel_consumed_mmbtu",
            ]
        ],
        how="left",
        on=["plant_id_eia", "emissions_unit_id_epa", "report_date"],
        validate="1:m",
    )
    # only keep observations with non-missing, non-zero co2 emissions and fuel
    # consumption
    unit_months_missing_co2 = unit_months_missing_co2[
        (unit_months_missing_co2["co2_mass_lb"] > 0)
        & (unit_months_missing_co2["fuel_consumed_mmbtu"] > 0)
    ]

    # calculate total fuel consumption and emissions by month
    unit_month_efs = (
        unit_months_missing_co2.groupby(
            ["plant_id_eia", "emissions_unit_id_epa", "report_date"], dropna=False
        )
        .sum()
        .reset_index()
    )
    unit_month_efs["co2_lb_per_mmbtu"] = (
        unit_month_efs["co2_mass_lb"] / unit_month_efs["fuel_consumed_mmbtu"]
    )

    # merge these EFs into the missing cems data
    missing_co2 = missing_co2.merge(
        unit_month_efs[
            ["plant_id_eia", "report_date", "emissions_unit_id_epa", "co2_lb_per_mmbtu"]
        ],
        how="left",
        on=["plant_id_eia", "report_date", "emissions_unit_id_epa"],
        validate="m:1",
    ).set_index(missing_co2.index)

    # only keep observations where there is a non-missing ef
    missing_co2 = missing_co2[~missing_co2["co2_lb_per_mmbtu"].isna()]

    # calculate missing co2 data
    missing_co2["co2_mass_lb"] = (
        missing_co2["fuel_consumed_mmbtu"] * missing_co2["co2_lb_per_mmbtu"]
    )

    # update in CEMS table
    cems.update(missing_co2[["co2_mass_lb"]])

    # update the co2 mass measurement code
    cems.loc[missing_co2.index, "co2_mass_measurement_code"] = "Imputed"

    # identify all observations that are still missing co2 data
    missing_co2 = cems[cems["co2_mass_lb"].isnull()]

    # Second round of data filling using weighted average emission factors based on
    # EIA-923 heat input data
    ####################################################################################

    # merge the weighted ef into the missing data
    missing_co2 = missing_co2.merge(
        subplant_emission_factors[
            ["plant_id_eia", "report_date", "subplant_id", "co2_lb_per_mmbtu"]
        ],
        how="left",
        on=["plant_id_eia", "report_date", "subplant_id"],
        validate="m:1",
    ).set_index(missing_co2.index)

    # only keep observations where there is a non-missing ef
    missing_co2 = missing_co2[~missing_co2["co2_lb_per_mmbtu"].isna()]

    # calculate missing co2 data
    missing_co2["co2_mass_lb"] = (
        missing_co2["fuel_consumed_mmbtu"] * missing_co2["co2_lb_per_mmbtu"]
    )

    # update in CEMS table
    cems.update(missing_co2[["co2_mass_lb"]])

    # update the co2 mass measurement code
    cems.loc[missing_co2.index, "co2_mass_measurement_code"] = "Imputed"

    # identify all observations that are still missing co2 data
    missing_co2 = cems[cems["co2_mass_lb"].isnull()]

    # Third round of filling using subplant fuel codes
    ##################################################

    # for rows that have a successful fuel code match, move to a temporary dataframe
    # to hold the data
    co2_to_fill = missing_co2[~missing_co2["energy_source_code"].isna()].copy()
    fill_index = co2_to_fill.index

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

    # check that there are no missing co2 values left
    still_missing_co2 = cems[cems["co2_mass_lb"].isna()]
    if len(still_missing_co2) > 0:
        raise UserWarning(
            "There are still misssing CO2 values remaining after filling missing CO2 "
            "values in CEMS"
        )

    # check that no non-missing co2 values were modified during filling
    validation.check_non_missing_cems_co2_values_unchanged(cems_original, cems)
    del cems_original
    gc.collect()

    return cems
