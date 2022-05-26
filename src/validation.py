import pandas as pd
import numpy as np


def load_egrid_plant_file(year):
    # load plant level data from egrid
    egrid_plant = pd.read_excel(
        f"../data/egrid/egrid{year}_data.xlsx",
        sheet_name=f"PLNT{str(year)[-2:]}",
        header=1,
        usecols=[
            "BACODE",
            "PSTATABB",
            "PLPRMFL",
            "ORISPL",
            "PNAME",
            "PLGENATN",
            "PLGENATR",
            "PLHTIANT",
            "UNCO2",
            "UNHTIT",
            "PLCO2AN",
        ],
    )
    # calculate total net generation from reported renewable and nonrenewable generation
    egrid_plant["net_generation_mwh"] = (
        egrid_plant["PLGENATN"] + egrid_plant["PLGENATR"]
    )
    egrid_plant = egrid_plant.drop(columns=["PLGENATN", "PLGENATR"])
    # rename the columns
    egrid_plant = egrid_plant.rename(
        columns={
            "BACODE": "ba_code",
            "PSTATABB": "state",
            "PLPRMFL": "energy_source_code",
            "ORISPL": "plant_id_egrid",
            "PNAME": "plant_name",
            "UNHTIT": "fuel_consumed_mmbtu",
            "PLHTIANT": "fuel_consumed_for_electricity_mmbtu",
            "UNCO2": "co2_mass_tons",
            "PLCO2AN": "co2_mass_tons_adjusted",
        }
    )

    # if egrid has a missing value for co2 for a clean plant, replace with zero
    clean_fuels = ["SUN", "MWH", "WND", "WAT", "WH", "PUR", "NUC"]
    egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_tons_adjusted"
    ] = egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_tons_adjusted"
    ].fillna(
        0
    )
    egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_tons"
    ] = egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_tons"
    ].fillna(
        0
    )

    # reorder the columns
    egrid_plant = egrid_plant[
        [
            "ba_code",
            "state",
            "plant_id_egrid",
            "plant_name",
            "energy_source_code",
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
            "co2_mass_tons",
            "co2_mass_tons_adjusted",
        ]
    ]

    # remove any plants that have no reported data
    # NOTE: it seems that egrid includes a lot of proposed projects that are not yet operating, but just has missing data for them
    plants_with_no_data_in_egrid = list(
        egrid_plant[
            egrid_plant[
                [
                    "net_generation_mwh",
                    "fuel_consumed_mmbtu",
                    "fuel_consumed_for_electricity_mmbtu",
                    "co2_mass_tons",
                    "co2_mass_tons_adjusted",
                ]
            ].sum(axis=1)
            == 0
        ]["plant_id_egrid"]
    )
    egrid_plant = egrid_plant[
        ~egrid_plant["plant_id_egrid"].isin(plants_with_no_data_in_egrid)
    ]

    # We also want to remove any plants that are located in Puerto Rico
    egrid_plant = egrid_plant[(egrid_plant["state"] != "PR")]

    # create a column for eia id
    egrid_plant = add_egrid_plant_id(egrid_plant, from_id="egrid", to_id="eia")

    return egrid_plant


def add_egrid_plant_id(df, from_id, to_id):
    # For plants that have different EPA and EIA plant IDs, the plant ID in eGRID is usually the EPA ID, but sometimes the EIA ID
    # however, there are sometime 2 EIA IDs for a single eGRID ID, so we need to group the data in the EIA table by the egrid id
    # We need to update all of the egrid plant IDs to the EIA plant IDs
    egrid_crosswalk = pd.read_csv(
        "../data/egrid/egrid_static_tables/2020/table_C5_crosswalk_of_EIA_ID_to_EPA_ID.csv"
    )
    id_map = dict(
        zip(
            list(egrid_crosswalk[f"plant_id_{from_id}"]),
            list(egrid_crosswalk[f"plant_id_{to_id}"]),
        )
    )

    df[f"plant_id_{to_id}"] = df[f"plant_id_{from_id}"]
    df[f"plant_id_{to_id}"].update(df[f"plant_id_{to_id}"].map(id_map))

    return df


def test_for_negative_values(df, columns_to_test):
    for column in columns_to_test:
        negative_test = df[df[column] < 0]
        if not negative_test.empty:
            print(
                f"Warning: There are {len(negative_test)} records where {column} is negative. Check `negative_test` for complete list"
            )
    return negative_test


def test_for_missing_fuel(df, generation_column):
    missing_fuel_test = df[
        (df[generation_column] > 0)
        & (
            (df["fuel_consumed_for_electricity_mmbtu"].isnull())
            | (df["fuel_consumed_for_electricity_mmbtu"] == 0)
        )
    ]
    if not missing_fuel_test.empty:
        print(
            f"Warning: There are {len(missing_fuel_test)} records where {generation_column} is positive but no fuel consumption is reported. Check `missing_fuel_test` for complete list"
        )

    return missing_fuel_test


def test_chp_allocation(df):
    chp_allocation_test = df[
        df["fuel_consumed_for_electricity_mmbtu"] > df["fuel_consumed_mmbtu"]
    ]
    if not chp_allocation_test.empty:
        print(
            f"Warning: There are {len(chp_allocation_test)} records where fuel consumed for electricity is greater than total fuel consumption. Check `chp_allocation_test` for complete list"
        )

    return chp_allocation_test


def test_for_missing_co2(df):
    missing_co2_test = df[
        df["co2_mass_tons"].isna() & ~df["fuel_consumed_mmbtu"].isna()
    ]
    if not missing_co2_test.empty:
        print(
            f"Warning: There are {len(missing_co2_test)} records where co2 data is missing. Check `missing_co2_test` for complete list"
        )
    return missing_co2_test


def test_for_missing_data(df, columns_to_test):
    missing_data_test = df[df[columns_to_test].isnull().all(axis=1)]
    if not missing_data_test.empty:
        print(
            f"Warning: There are {len(missing_data_test)} records for which no data was reported. Check `missing_data_test` for complete list"
        )
    return missing_data_test


def test_for_outlier_heat_rates(df):
    # check heat rates
    print("Heat Rate Test")
    # remove non-fossil fuel types
    thermal_generators = df[
        ~df["energy_source_code"].isin(["SUN", "MWH", "WND", "WAT", "WH", "PUR"])
    ]
    for fuel_type in sorted(
        list(thermal_generators.energy_source_code.dropna().unique())
    ):
        generators = thermal_generators[
            thermal_generators["energy_source_code"] == fuel_type
        ]
        heat_rate = (
            generators["fuel_consumed_for_electricity_mmbtu"]
            / generators["net_generation_mwh"]
        )
        heat_rate_stats = heat_rate[(heat_rate >= 0) & (heat_rate != np.inf)].describe()
        outlier_threshold = (
            (heat_rate_stats["75%"] - heat_rate_stats["25%"]) * 1.5
        ) + heat_rate_stats["75%"]
        heat_rate_test = generators[
            (heat_rate > outlier_threshold)
            & (heat_rate != np.inf)
            & (generators["net_generation_mwh"] >= 1)
        ]
        if not heat_rate_test.empty:
            print(
                f'    Warning: {len(heat_rate_test)} of {len(generators)} {fuel_type} records have heat rate > {outlier_threshold.round(2)} mmbtu/MWh (median = {heat_rate_stats["50%"].round(2)})'
            )


def test_for_missing_energy_source_code(df):
    missing_esc_test = df[df["energy_source_code"].isna()]
    if not missing_esc_test.empty:
        print(
            f"Warning: There are {len(missing_esc_test)} records where there is a missing energy source code. Check `missing_esc_test` for complete list"
        )

    return missing_esc_test


def test_for_zero_data(df, columns_to_test):
    zero_data_test = df[
        (~df[columns_to_test].isnull().all(axis=1))
        & (df[columns_to_test].sum(axis=1) == 0)
    ]
    if not zero_data_test.empty:
        print(
            f"Warning: There are {len(zero_data_test)} records where all operating data are zero. Check `zero_data_test` for complete list"
        )
    return zero_data_test


def test_for_missing_subplant_id(df):
    missing_subplant_test = df[df["subplant_id"].isna()]
    if not missing_subplant_test.empty:
        print(
            f"Warning: There are {len(missing_subplant_test)} records without a subplant ID. See `missing_subplant_test` for details"
        )
    return missing_subplant_test


def test_gtn_results(df):
    gtn_test = df[df["net_generation_mwh"] > df["gross_generation_mwh"]]
    if not gtn_test.empty:
        print(
            f"Warning: There are {len(gtn_test)} records where net generation > gross generation. See `gtn_test` for details"
        )
    return gtn_test
