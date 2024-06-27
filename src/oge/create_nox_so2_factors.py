import pandas as pd

from oge.download_data import download_helper
from oge.filepaths import downloads_folder
from oge.logging_util import get_logger
from oge.constants import ConversionFactors

logger = get_logger(__name__)


def generate_nox_emission_factor_reference_table():
    # download the nox data from EIA
    download_helper(
        "https://www.eia.gov/electricity/annual/xls/epa_a_02.xlsx",
        downloads_folder("eia/table_a2_nox_uncontrolled_emission_factors.xlsx"),
    )

    nox_factors_eia = pd.read_excel(
        downloads_folder("eia/table_a2_nox_uncontrolled_emission_factors.xlsx"),
        header=3,
        names=[
            "fuel",
            "energy_source_code",
            "data_source",
            "unit",
            "cyclone_firing",
            "fluidized_bed_firing",
            "stoker",
            "tangential_firing_dry",
            "tangential_firing_wet",
            "other_dry",
            "other_wet",
            "GT",
            "IC",
        ],
        skipfooter=1,
    )

    factor_columns = [
        "cyclone_firing",
        "fluidized_bed_firing",
        "stoker",
        "tangential_firing_dry",
        "tangential_firing_wet",
        "other_dry",
        "other_wet",
        "GT",
        "IC",
    ]
    nox_factors_eia = standardize_emission_factor_units(nox_factors_eia, factor_columns)

    nox_factors_eia = fill_out_missing_nox_so2_factor_combinations(
        nox_factors_eia, pollutant="nox"
    )

    bft_columns = [
        "cyclone_firing",
        "cyclone_firing_dry",
        "cyclone_firing_wet",
        "fluidized_bed_firing",
        "fluidized_bed_firing_dry",
        "fluidized_bed_firing_wet",
        "stoker",
        "stoker_dry",
        "stoker_wet",
        "tangential_firing_dry",
        "tangential_firing_wet",
        "tangential_firing_none",
        "other_dry",
        "other_wet",
        "other_none",
        "cell_burner_dry",
        "cell_burner_wet",
        "cell_burner_none",
        "duct_burner_dry",
        "duct_burner_wet",
        "duct_burner_none",
        "vertical_firing_dry",
        "vertical_firing_wet",
        "vertical_firing_none",
        "wall_fired_dry",
        "wall_fired_wet",
        "wall_fired_none",
        "none_dry",
        "none_wet",
    ]
    nox_factors_eia_bft = create_long_bft_specific_df(
        nox_factors_eia, bft_columns, pollutant="nox"
    )

    pm_columns = ["GT", "IC", "CA", "CC", "CS", "CT", "CE"]

    nox_factors_eia_pm = create_long_pm_specific_df(
        nox_factors_eia, pm_columns, pollutant="nox"
    )

    # concat these back toghether and organize columns
    nox_factors_eia_combined = pd.concat(
        [nox_factors_eia_bft, nox_factors_eia_pm], axis=0
    )
    nox_factors_eia_combined = nox_factors_eia_combined[
        [
            "prime_mover_code",
            "energy_source_code",
            "wet_dry_bottom",
            "boiler_firing_type",
            "emission_factor",
            "emission_factor_numerator",
            "emission_factor_denominator",
            # "data_source",
        ]
    ].sort_values(by=["energy_source_code", "prime_mover_code", "boiler_firing_type"])

    return nox_factors_eia_combined


def generate_so2_emission_factor_reference_table():
    # download the nox data from EIA
    download_helper(
        "https://www.eia.gov/electricity/annual/xls/epa_a_01.xlsx",
        downloads_folder("eia/table_a1_so2_uncontrolled_emission_factors.xlsx"),
    )

    so2_factors_eia = pd.read_excel(
        downloads_folder("eia/table_a1_so2_uncontrolled_emission_factors.xlsx"),
        header=3,
        names=[
            "fuel",
            "energy_source_code",
            "data_source",
            "unit",
            "cyclone_firing",
            "fluidized_bed_firing",
            "stoker",
            "tangential_firing",
            "other",
            "GT",
            "IC",
        ],
        skipfooter=1,
    )

    factor_columns = [
        "cyclone_firing",
        "fluidized_bed_firing",
        "stoker",
        "tangential_firing",
        "other",
        "GT",
        "IC",
    ]
    so2_factors_eia = standardize_emission_factor_units(so2_factors_eia, factor_columns)

    so2_factors_eia = fill_out_missing_nox_so2_factor_combinations(
        so2_factors_eia, pollutant="so2"
    )

    bft_columns = [
        "cyclone_firing",
        "fluidized_bed_firing",
        "stoker",
        "tangential_firing",
        "other",
        "cell_burner",
        "duct_burner",
        "vertical_firing",
        "wall_fired",
        "none",
    ]
    so2_factors_eia_bft = create_long_bft_specific_df(
        so2_factors_eia, bft_columns, pollutant="so2"
    )

    pm_columns = ["GT", "IC", "CA", "CC", "CS", "CT", "CE"]
    so2_factors_eia_pm = create_long_pm_specific_df(
        so2_factors_eia, pm_columns, pollutant="so2"
    )

    # concat these back toghether and organize columns
    so2_factors_eia_combined = pd.concat(
        [so2_factors_eia_bft, so2_factors_eia_pm], axis=0
    )
    so2_factors_eia_combined = so2_factors_eia_combined[
        [
            "prime_mover_code",
            "energy_source_code",
            "boiler_firing_type",
            "emission_factor",
            "emission_factor_numerator",
            "emission_factor_denominator",
            # "data_source",
        ]
    ].sort_values(by=["energy_source_code", "prime_mover_code", "boiler_firing_type"])

    return so2_factors_eia_combined


def standardize_emission_factor_units(
    df: pd.DataFrame, factor_columns: list[str]
) -> pd.DataFrame:
    """Splits the "unit" column (like "Lbs per Mcf") in the so2 or nox factor table into
    two separate numerator (eg "Lbs") and denominator (eg "Mcf") columns, renames "Lbs"
    to "lb", and converts the `factor_columns` to standardized units that are used by
    EIA-923:
     - MMCF (million cubic feet) converted to mcf (thousand cubic feet)
     - MG (thousand gallons) converted to (oil) barrels

    Args:
        df (pd.DataFrame): the raw nox or so2 table
        factor_columns (list[str]): the columns to convert

    Returns:
        pd.DataFrame: original dataframe with new unit columns with standardized units,
            and factor_columns converted to the proper units
    """
    # standardize the unit column
    df["emission_factor_numerator"] = (
        df["unit"].str.split(" ", expand=True)[0].replace("Lbs", "lb")
    )
    df["emission_factor_denominator"] = df["unit"].str.split(" ", expand=True)[2]

    # convert the denominator units into units used by EIA-923
    # convert million cubic feet (MMCF) to thousand cubic feet (Mcf)
    df.loc[df["emission_factor_denominator"] == "MMCF", factor_columns] /= 1000
    df.loc[
        df["emission_factor_denominator"] == "MMCF",
        "emission_factor_denominator",
    ] = "mcf"

    # convert thousand gallons (MG) to barrels
    df.loc[df["emission_factor_denominator"] == "MG", factor_columns] /= (
        ConversionFactors.kgal_to_barrel
    )
    df.loc[
        df["emission_factor_denominator"] == "MG",
        "emission_factor_denominator",
    ] = "barrel"

    return df


def fill_out_missing_nox_so2_factor_combinations(df, pollutant):
    """Fills in Nox or SO2 factors that are not explicitly specified by the EIA table,
    applying a set of assumptions about how the data should be filled.

    There are several filling steps:
    1. Missing energy_source_codes: The EIA tables do not provide NOx or SO2 factors for
    MSB, MSN, or SC fuel types. For MSB and MSN (the biogenic and nonbiogenic
    of MSW), we set the factors equal to the specified MSW factor. For SC (Coal Synfuel)
    We set the factor equal to the refined coal (RC) factor, which itself is assumed to
    be equal to the bituminous coal (BIT) factor.
    2. Factors for units missing boiler_firing_type (BFT) and wet_dry_bottom (WDB) data:
    According to the technical appendix for the EIA's Electric Power Annual (see:
    https://www.eia.gov/electricity/annual/pdf/tech_notes.pdf), whenever a unit is
    missing BFT and WDB data (in the case of NOx), they apply they assume it is a dry-
    bottom boiler with an "other" BFT. We apply the same assumption, creating a "none"
    BFT.
    3. Missing GT-specific factors: The EIA tables do not report GT-specific factors for
    certain solid fuels (like coal), which cannot be combusted in a combustion turbine.
    However, these fuel-GT combinations do occasionally show up in the data. Rather than
    leave these missing (and thus calculate no emissions), we fill missing GT factors
    using the "other" boiler firing type (BFT), or in the case of NOx, the "other" BFT
    with a dry bottom.
    4. The provided combustion-turbine specific factors apply to a range of prime movers
    that use combustion turbines, including GT, CA, CC, CS, and CT. This also applies to
    the CE prime mover code, which applies to a specific plant (7063) which is a gas
    turbine in combination with compressed air energy storage.
    5. (NOx-only) Fill missing WDB combinations for BFT-only factors: For certain BFTS,
    two factors are listed for each WDB type, for example (tangential_firing, wet),
    (tangential_firing, dry). For BFTs were WDB-specific factors are not specified, we
    treat the WDB as "none". To ensure a complete set of BFT-WDB combinations, we need
    to add "wet" and "dry" WDB options where unspecified, and "none" WBD where specified
    so that there is a wet, dry, and none option for every BFT.
    6. Besides the BFTs referenced in the EIA factor tables, there are four other
    specific BFT types enumerated in EIA data elsewhere: cell_burner, duct_burner,
    wall_fired, vertical_firing. Since these are literally "other" BFTs, we apply the
    appropriate "other" factor to each of these BFTs.

    The EIA-provided factors are listed for a combination of boiler firing types (BFT) and prime movers.
    Of the several dozen prime movers, only four types are associated with thermal generators that combust fuel:
    steam turbines (ST), combustion turbines (GT), combined cycle turbines (CS, CT, CA, CC), internal combustion engines (IC), and possibly other (OT).
    Prime-mover specific factors are provided for GT and IC prime movers, and are assumed to apply to all boiler firing types, if provided.
    Because combined cycle plants combust fuel in the combustion turbine part, GT factors are applied to all combined cycle prime mover components.
    The BFT-specific factors are used only for ST and OT prime movers
    """
    if pollutant == "nox":
        default_value = "other_dry"
        wet_dry_bottom_types = ["_wet", "_dry", "_none"]
    elif pollutant == "so2":
        default_value = "other"
        wet_dry_bottom_types = [""]

    # fill in missing energy_source_codes
    df = df.set_index("energy_source_code")
    for new_esc, esc_to_copy in {"MSB": "MSW", "MSN": "MSW", "SC": "RC"}.items():
        df.loc[new_esc, :] = df.loc[esc_to_copy, :]
    df = df.reset_index()

    # create a value for missing BFT/WDB using the default value
    df["none"] = df[default_value]

    # Fill missing GT factors using the default value
    df["GT"] = df["GT"].fillna(df[default_value])

    # apply the combustion turbine factor to all combustion-turbine-type prime movers,
    # including those at combined cycle plants, and CAES plants
    for pm in ["CA", "CC", "CS", "CT", "CE"]:
        df[pm] = df["GT"]

    # where we have a prime mover with no BFT or WDB information, we assume a dry-bottom
    # boiler in the "all other" BFT category. This is consistent with the assumption made
    # by the EIA. see: https://www.eia.gov/electricity/annual/pdf/tech_notes.pdf
    # for pm in ["ST","OT"]:
    #    df[pm] = df[default_value]

    # (for NOx) we need to fill out all possible combinations of BFT and WDB
    if pollutant == "nox":
        # for columns where BFT is specified without difference to WDB, we set all WDB options
        # equal to each other so that wet = dry = none
        for col in [
            "cyclone_firing",
            "fluidized_bed_firing",
            "stoker",
        ]:
            df[f"{col}_wet"] = df[col]
            df[f"{col}_dry"] = df[col]

        # for columns where WDB is specified, we set the "none" WDB equal to the dry bottom
        # parameter. This is consistent with the assumption EIA makes:
        # see: https://www.eia.gov/electricity/annual/pdf/tech_notes.pdf
        for col in ["tangential_firing", "other"]:
            df[f"{col}_none"] = df[f"{col}_dry"]

    # add columns for specific boiler firing types that are not mentioned specifically
    # use the "All other Boiler Types" column for this. For nox, use the appropriate
    # wet dry bottom value
    other_firing_types = [
        "cell_burner",
        "duct_burner",
        "vertical_firing",
        "wall_fired",
    ]
    for bft in other_firing_types:
        for wdb in wet_dry_bottom_types:
            df[f"{bft}{wdb}"] = df[f"other{wdb}"]
    # for nox, also create wet and dry factors for "none" bft.
    # (none, none) factors will be handled as a prime-mover specific factor
    if pollutant == "nox":
        df["none_wet"] = df["other_wet"]
        df["none_dry"] = df["other_dry"]

    return df


def create_long_bft_specific_df(df, bft_columns, pollutant):
    # separate the data by boiler-specific and PM-specific factors
    df_bft = df[
        [
            "energy_source_code",
            # "data_source",
            "emission_factor_numerator",
            "emission_factor_denominator",
        ]
        + bft_columns
    ].copy()

    # melt the data into long form
    df_bft = df_bft.melt(
        id_vars=[
            "energy_source_code",
            # "data_source",
            "emission_factor_numerator",
            "emission_factor_denominator",
        ],
        var_name="boiler_firing_type",
        value_name="emission_factor",
    ).dropna()

    # for nox only
    # create a wet_dry_bottom column and remove this suffix from the boiler firing type
    if pollutant == "nox":
        df_bft["wet_dry_bottom"] = "none"
        df_bft["boiler_firing_type"] = df_bft["boiler_firing_type"].str.replace(
            "_none", ""
        )
        df_bft.loc[
            df_bft["boiler_firing_type"].str.contains("_wet"), "wet_dry_bottom"
        ] = "wet"
        df_bft.loc[
            df_bft["boiler_firing_type"].str.contains("_dry"), "wet_dry_bottom"
        ] = "dry"
        df_bft["boiler_firing_type"] = df_bft["boiler_firing_type"].str.replace(
            "_dry", ""
        )
        df_bft["boiler_firing_type"] = df_bft["boiler_firing_type"].str.replace(
            "_wet", ""
        )

    # add prime mover columns to this data
    # these BFT-specific factors generally only apply to ST and OT prime movers.
    # However, sometimes combined cycle units can have duct burners or other boilers as part
    # of a HRSG (https://www.eia.gov/todayinenergy/detail.php?id=52778). Thus, we will also
    # add factors for these prime movers as well,
    # To add all at the same time we will create a dataframe with both and merge it in to
    # easily duplicate for each row
    df_bft = df_bft.merge(
        pd.DataFrame(
            columns=["emission_factor_numerator", "prime_mover_code"],
            data=[
                ["lb", "ST"],
                ["lb", "OT"],
                ["lb", "GT"],
                ["lb", "CA"],
                ["lb", "CS"],
                ["lb", "CT"],
            ],
        ),
        how="left",
        on="emission_factor_numerator",
    )

    return df_bft


def create_long_pm_specific_df(df, pm_columns, pollutant):
    # create a prime-mover specific table
    df_pm = df[
        [
            "energy_source_code",
            "data_source",
            "emission_factor_numerator",
            "emission_factor_denominator",
        ]
        + pm_columns
    ].copy()

    df_pm = df_pm.melt(
        id_vars=[
            "energy_source_code",
            "data_source",
            "emission_factor_numerator",
            "emission_factor_denominator",
        ],
        var_name="prime_mover_code",
        value_name="emission_factor",
    ).dropna()

    if pollutant == "nox":
        # add a column for the boiler firing type and set as "none"
        df_pm["boiler_firing_type"] = "none"
        df_pm["wet_dry_bottom"] = "none"

    return df_pm
