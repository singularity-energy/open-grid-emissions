import pandas as pd
import numpy as np

import src.load_data as load_data
import src.column_checks as column_checks

GENERATED_EMISSION_RATE_COLS = [
    "generated_co2_rate_lb_per_mwh_for_electricity",
    "generated_ch4_rate_lb_per_mwh_for_electricity",
    "generated_n2o_rate_lb_per_mwh_for_electricity",
    "generated_nox_rate_lb_per_mwh_for_electricity",
    "generated_so2_rate_lb_per_mwh_for_electricity",
    "generated_co2_rate_lb_per_mwh_adjusted",
    "generated_ch4_rate_lb_per_mwh_adjusted",
    "generated_n2o_rate_lb_per_mwh_adjusted",
    "generated_nox_rate_lb_per_mwh_adjusted",
    "generated_so2_rate_lb_per_mwh_adjusted",
]


def output_intermediate_data(df, file_name, path_prefix, year):

    df.to_csv(f"../data/outputs/{path_prefix}{file_name}_{year}.csv", index=False)
    column_checks.check_columns(f"../data/outputs/{path_prefix}{file_name}_{year}.csv")


def output_to_results(df, file_name, subfolder, path_prefix):

    df.to_csv(f"../data/results/{path_prefix}{subfolder}{file_name}.csv", index=False)


def write_generated_averages(ba_fuel_data, path_prefix, year):
    avg_fuel_type_production = (
        ba_fuel_data.groupby(["fuel_category"]).sum().reset_index()
    )
    # Add row for total before taking rates
    total = avg_fuel_type_production.mean(numeric_only=True).to_frame().T
    total.loc[0, "fuel_category"] = "total"
    avg_fuel_type_production = pd.concat([avg_fuel_type_production, total], axis=0)

    # Find rates
    for emission_type in ["_for_electricity", "_adjusted"]:
        for emission in ["co2", "ch4", "n2o", "nox", "so2"]:
            avg_fuel_type_production[
                f"generated_{emission}_rate_lb_per_mwh{emission_type}"
            ] = (
                (
                    avg_fuel_type_production[f"{emission}_mass_lb{emission_type}"]
                    / avg_fuel_type_production["net_generation_mwh"]
                )
                .fillna(0)
                .replace(np.inf, np.NaN)
                .replace(-np.inf, np.NaN)
                .replace(
                    np.NaN, 0
                )  # TODO: temporary placeholder while solar is broken. Eventually there should be no NaNs.
            )
    output_intermediate_data(
        avg_fuel_type_production,
        "annual_generation_averages_by_fuel",
        path_prefix,
        year,
    )


def write_power_sector_results(ba_fuel_data, path_prefix):
    """
    Helper function to write combined data by BA
    """

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

    for ba in list(ba_fuel_data.ba_code.unique()):

        # filter the data for a single BA
        ba_table = ba_fuel_data[ba_fuel_data["ba_code"] == ba].drop(columns="ba_code")

        # convert the datetime_utc column back to a datetime
        ba_table["datetime_utc"] = pd.to_datetime(ba_table["datetime_utc"], utc=True)

        # calculate a total for the BA
        ba_total = (
            ba_table.groupby(["datetime_utc"], dropna=False)
            .sum()[data_columns]
            .reset_index()
        )
        ba_total["fuel_category"] = "total"

        # concat the totals to the fuel-specific totals
        ba_table = pd.concat([ba_table, ba_total], axis=0, ignore_index=True)

        # round all values to one decimal place
        ba_table = ba_table.round(2)

        for emission_type in ["_for_electricity", "_adjusted"]:
            for emission in ["co2", "ch4", "n2o", "nox", "so2"]:
                ba_table[f"generated_{emission}_rate_lb_per_mwh{emission_type}"] = (
                    (
                        ba_table[f"{emission}_mass_lb{emission_type}"]
                        / ba_table["net_generation_mwh"]
                    )
                    .fillna(0)
                    .replace(np.inf, np.NaN)
                    .replace(-np.inf, np.NaN)
                )

        # create a local datetime column
        try:
            local_tz = load_data.ba_timezone(ba, "local")
            ba_table["datetime_local"] = ba_table["datetime_utc"].dt.tz_convert(
                local_tz
            )
        # TODO: figure out what to do for missing ba
        except ValueError:
            ba_table["datetime_local"] = pd.NaT

        # re-order columns
        ba_table = ba_table[
            ["fuel_category", "datetime_local", "datetime_utc"]
            + data_columns
            + GENERATED_EMISSION_RATE_COLS
        ]

        # export to a csv
        ba_table.to_csv(
            f"../data/results/{path_prefix}power_sector_data/{ba}.csv", index=False
        )
