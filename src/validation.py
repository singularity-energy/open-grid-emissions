import pandas as pd
import numpy as np

import src.load_data as load_data


def load_egrid_plant_file(year):
    # load plant level data from egrid
    egrid_plant = pd.read_excel(
        f"../data/downloads/egrid/egrid{year}_data.xlsx",
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
            "UNCO2": "co2_mass_lb",  # this is actually in tons, but we are converting in the next step
            "PLCO2AN": "co2_mass_lb_adjusted",  # this is actually in tons, but we are converting in the next step
        }
    )

    # convert co2 mass tons to lb
    egrid_plant["co2_mass_lb"] = egrid_plant["co2_mass_lb"] * 2000
    egrid_plant["co2_mass_lb_adjusted"] = egrid_plant["co2_mass_lb_adjusted"] * 2000

    # if egrid has a missing value for co2 for a clean plant, replace with zero
    clean_fuels = ["SUN", "MWH", "WND", "WAT", "WH", "PUR", "NUC"]
    egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_lb_adjusted"
    ] = egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_lb_adjusted"
    ].fillna(
        0
    )
    egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_lb"
    ] = egrid_plant.loc[
        egrid_plant["energy_source_code"].isin(clean_fuels), "co2_mass_lb"
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
            "co2_mass_lb",
            "co2_mass_lb_adjusted",
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
                    "co2_mass_lb",
                    "co2_mass_lb_adjusted",
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
        "../data/manual/egrid_static_tables/table_C5_crosswalk_of_EIA_ID_to_EPA_ID.csv"
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
    missing_co2_test = df[df["co2_mass_lb"].isna() & ~df["fuel_consumed_mmbtu"].isna()]
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


def test_for_missing_incorrect_prime_movers(df, year):

    # cehck for incorrect PM by comparing to EIA-860 data
    pudl_out = load_data.initialize_pudl_out(year)
    pms_in_eia860 = pudl_out.gens_eia860()[
        ["plant_id_eia", "generator_id", "prime_mover_code"]
    ]
    incorrect_pm_test = df.copy()[["plant_id_eia", "generator_id", "prime_mover_code"]]
    incorrect_pm_test = incorrect_pm_test.merge(
        pms_in_eia860,
        how="left",
        on=["plant_id_eia", "generator_id"],
        suffixes=("_allocated", "_eia860"),
    )
    incorrect_pm_test[
        incorrect_pm_test["prime_mover_code_allocated"]
        != incorrect_pm_test["prime_mover_code_eia860"]
    ]
    if not incorrect_pm_test.empty:
        print(
            f"Warning: There are {len(incorrect_pm_test)} records for which the allocated prime mover does not match the reported prime mover. Check `incorrect_pm_test` for complete list"
        )

    # check for missing PM code
    missing_pm_test = df[df["prime_mover_code"].isna()]
    if not missing_pm_test.empty:
        print(
            f"Warning: There are {len(missing_pm_test)} records for which no prime mover was assigned. Check `missing_pm_test` for complete list"
        )

    return incorrect_pm_test, missing_pm_test


def test_for_outlier_heat_rates(df):
    # check heat rates
    print("Heat Rate Test")
    # remove non-fossil fuel types
    thermal_generators = df[
        ~df["energy_source_code"].isin(["SUN", "MWH", "WND", "WAT", "WH", "PUR"])
    ]
    heat_rate_test_all = []
    for fuel_type in sorted(
        list(thermal_generators.energy_source_code.dropna().unique())
    ):
        # identify all generators with a given fuel type
        generators = thermal_generators[
            thermal_generators["energy_source_code"] == fuel_type
        ]
        # identify all unique prime mover codes for generators of that fuel type
        for pm in sorted(list(generators.prime_mover_code.dropna().unique())):
            generators_with_pm = generators[generators["prime_mover_code"] == pm]
            # calculate a heat rate for each generator of the given fuel type
            heat_rate = (
                generators_with_pm["fuel_consumed_for_electricity_mmbtu"]
                / generators_with_pm["net_generation_mwh"]
            )
            # calculate descriptive statistics for all nonnegative heat rates
            heat_rate_stats = heat_rate[
                (heat_rate >= 0) & (heat_rate != np.inf)
            ].describe()
            # set the outlier threhshold to 1.5 x IQR
            outlier_threshold = (
                (heat_rate_stats["75%"] - heat_rate_stats["25%"]) * 1.5
            ) + heat_rate_stats["75%"]
            # identify all generators whose heatrate is outside this threshold
            heat_rate_test = generators_with_pm[
                (
                    (heat_rate > outlier_threshold)
                    & (heat_rate != np.inf)
                    & (generators_with_pm["net_generation_mwh"] >= 1)
                    | (heat_rate.round(2) == 0)
                )
            ]
            if not heat_rate_test.empty:
                print(
                    f"    Warning: {len(heat_rate_test)} of {len(generators_with_pm)} records for {fuel_type} generators with {pm} prime mover have heat rate of zero or > {outlier_threshold.round(2)} mmbtu/MWh"
                )
                print(
                    f'             median = {heat_rate_stats["50%"].round(2)}, max = {heat_rate_stats["max"].round(2)}, min = {heat_rate_stats["min"].round(2)}'
                )
                heat_rate_test_all.append(heat_rate_test)

    heat_rate_test_all = pd.concat(heat_rate_test_all, axis=0)[
        [
            "report_date",
            "plant_id_eia",
            "generator_id",
            "energy_source_code",
            "prime_mover_code",
            "net_generation_mwh",
            "fuel_consumed_mmbtu",
            "fuel_consumed_for_electricity_mmbtu",
        ]
    ]
    heat_rate_test_all["heat_rate"] = (
        heat_rate_test_all["fuel_consumed_for_electricity_mmbtu"]
        / heat_rate_test_all["net_generation_mwh"]
    )
    return heat_rate_test_all


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
            f"Warning: There are {round(len(gtn_test)/len(df)*100, 1)}% of records where net generation > gross generation. See `gtn_test` for details"
        )
    return gtn_test


def co2_source_metric(cems, partial_cems, shaped_eia_data):
    """Calculates what percent of CO2 emissions mass came from each source."""
    # determine the source of the co2 data
    co2_from_eia = (
        partial_cems["co2_mass_lb"].sum() + shaped_eia_data["co2_mass_lb"].sum()
    )
    co2_from_eia = pd.DataFrame(
        [{"co2_mass_measurement_code": "EIA Calculated", "co2_mass_lb": co2_from_eia}]
    )

    co2_from_cems = (
        cems.groupby("co2_mass_measurement_code", dropna=False)["co2_mass_lb"].sum().reset_index()
    )
    co2_from_cems["co2_mass_measurement_code"] = "CEMS " + co2_from_cems[
        "co2_mass_measurement_code"
    ].astype(str)

    co2_source = pd.concat([co2_from_cems, co2_from_eia])
    co2_source["percent"] = (
        co2_source["co2_mass_lb"] / co2_source["co2_mass_lb"].sum() * 100
    )
    co2_source = co2_source.round(2)
    return co2_source


def net_generation_method_metric(cems, partial_cems, shaped_eia_data):
    """Calculates what percent of net generation mwh was calculated using each method."""
    # determine the method for the net generation data
    data_metric = "net_generation_mwh"

    eia_ng_method = (
        shaped_eia_data.groupby("profile_method", dropna=False)[data_metric]
        .sum()
        .reset_index()
        .rename(columns={"profile_method": "method"})
    )
    cems_ng_method = (
        cems.groupby("gtn_method", dropna=False)[data_metric]
        .sum()
        .reset_index()
        .rename(columns={"gtn_method": "method"})
    )
    partial_cems_ng_method = partial_cems[data_metric].sum()
    partial_cems_ng_method = pd.DataFrame(
        [{"method": "partial_cems", data_metric: partial_cems_ng_method}]
    )

    ng_method = pd.concat([cems_ng_method, partial_cems_ng_method, eia_ng_method])
    ng_method["percent"] = ng_method[data_metric] / ng_method[data_metric].sum() * 100
    ng_method = ng_method.round(2)
    return ng_method


def hourly_profile_source_metric(cems, partial_cems, shaped_eia_data):
    """Calculates the percentage of data whose hourly profile was determined by method"""
    data_metric = "co2_mass_lb"

    # determine the source of the hourly profile
    profile_from_cems = cems[data_metric].sum()
    profile_from_partial_cems = partial_cems[data_metric].sum()
    profile_from_eia = (
        shaped_eia_data.groupby("profile_method", dropna=False)[data_metric].sum().reset_index()
    )

    profile_from_cems = pd.DataFrame(
        [
            {"profile_method": "cems", data_metric: profile_from_cems},
            {"profile_method": "partial_cems", data_metric: profile_from_partial_cems},
        ]
    )

    profile_source = pd.concat([profile_from_cems, profile_from_eia])
    profile_source["percent"] = (
        profile_source[data_metric] / profile_source[data_metric].sum() * 100
    )
    profile_source = profile_source.round(2)
    return profile_source
