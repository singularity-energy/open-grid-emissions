"""Helper functions for visualization."""
import pandas as pd
import plotly.express as px

from filepaths import downloads_folder, results_folder

"""
From a timeseries, make a day-hour array with hours in EST. Also return day labels for plotting.

Assumes timeseries is hourly, because we don't know how to resample to hourly if it's not.
"""


def day_hour_heatmap(timeseries: pd.Series, year: int = 2022):
    timeseries = timeseries.rename("a")  # needs to be named so we can merge

    # Make a new dataframe that has all full days in year in EST
    timeseries.index = timeseries.index.tz_convert("EST")
    hours_index = pd.DataFrame(
        index=pd.date_range(
            f"{year}-01-01 T00:00", f"{year}-12-31 T23:00", freq="H"
        ).tz_localize("EST")
    )
    hours_index = hours_index.merge(
        timeseries, how="left", left_index=True, right_index=True
    )

    # Get times for x labels
    times = hours_index.resample("D").first().index.to_numpy()
    hours = hours_index.to_numpy()
    hours = hours.reshape((24, len(hours) // 24), order="F")

    return (hours, times)


def graph_hourly_data_by_fuel_category(
    hourly_data, ba, column_name, fuel_category_name, plot_type
):
    fuel_color = {
        "natural_gas": "sienna",
        "coal": "black",
        "nuclear": "red",
        "biomass": "green",
        "geothermal": "orange",
        "wind": "blue",
        "solar": "gold",
        "petroleum": "purple",
        "hydro": "skyblue",
        "other": "lightgrey",
        "waste": "pink",
    }

    fuel_order = [
        "nuclear",
        "geothermal",
        "hydro",
        "other",
        "coal",
        "biomass",
        "petroleum",
        "waste",
        "solar",
        "wind",
        "natural_gas",
    ]

    hourly_data = hourly_data[hourly_data[fuel_category_name] != "total"]

    if plot_type == "area":
        plot = px.area(
            hourly_data,
            x="datetime_utc",
            y=column_name,
            color=fuel_category_name,
            color_discrete_map=fuel_color,
            template="plotly_white",
            title=f"Hourly {column_name} data for {ba} by fuel type",
            category_orders={fuel_category_name: fuel_order},
        ).update_traces(line={"width": 0})
    elif plot_type == "line":
        plot = px.line(
            hourly_data,
            x="datetime_utc",
            y=column_name,
            color=fuel_category_name,
            color_discrete_map=fuel_color,
            template="plotly_white",
            title=f"Hourly {column_name} data for {ba} by fuel type",
            category_orders={fuel_category_name: fuel_order},
        ).update_traces(line={"width": 0})

    return plot


def load_ba_ef_data_to_graph(ba, year, pollutant, rate_type, show_egrid):
    if show_egrid:
        egrid_ba = pd.read_excel(
            f"{downloads_folder()}egrid/egrid{year}_data.xlsx",
            sheet_name=f"BA{str(year)[-2:]}",
            header=1,
            usecols=[
                "BANAME",
                "BACODE",
                "BANOXRTO",
                "BASO2RTA",
                "BACO2RTA",
                "BAC2ERTA",
            ],
        ).rename(
            columns={
                "BANAME": "ba_name",
                "BACODE": "ba_code",
                "BANOXRTO": "generated_nox_rate_lb_per_mwh_for_electricity_adjusted",
                "BASO2RTA": "generated_so2_rate_lb_per_mwh_for_electricity_adjusted",
                "BACO2RTA": "generated_co2_rate_lb_per_mwh_for_electricity_adjusted",
                "BAC2ERTA": "generated_co2e_rate_lb_per_mwh_for_electricity_adjusted",
            }
        )

    hourly_consumed = pd.read_csv(
        f"{results_folder()}{year}/carbon_accounting/hourly/us_units/{ba}.csv",
        parse_dates=["datetime_local"],
    ).set_index("datetime_local")
    hourly_produced = pd.read_csv(
        f"{results_folder()}{year}/power_sector_data/hourly/us_units/{ba}.csv",
        parse_dates=["datetime_local"],
    ).set_index("datetime_local")

    data_to_graph = hourly_consumed[
        [f"consumed_{pollutant}_rate_lb_per_mwh_{rate_type}"]
    ].rename(
        columns={f"consumed_{pollutant}_rate_lb_per_mwh_{rate_type}": "hourly_consumed"}
    )
    data_to_graph["month"] = data_to_graph.index.astype(str).str[:7]
    data_to_graph["annual_consumed"] = data_to_graph[["hourly_consumed"]].mean().item()
    monthly_consumed = (
        data_to_graph.groupby("month")["hourly_consumed"]
        .mean()
        .reset_index()
        .rename(columns={"hourly_consumed": "monthly_consumed"})
    )
    orig_index = data_to_graph.index
    data_to_graph = data_to_graph.merge(
        monthly_consumed, how="left", on="month", validate="m:1"
    ).set_index(orig_index)
    data_to_graph = data_to_graph.merge(
        hourly_produced.loc[
            hourly_produced["fuel_category"] == "total",
            f"generated_{pollutant}_rate_lb_per_mwh_{rate_type}",
        ],
        how="inner",
        left_index=True,
        right_index=True,
    ).rename(
        columns={
            f"generated_{pollutant}_rate_lb_per_mwh_{rate_type}": "hourly_produced"
        }
    )
    data_to_graph["annual_produced"] = data_to_graph[["hourly_produced"]].mean().item()
    if show_egrid:
        egrid_value = egrid_ba.loc[
            egrid_ba["ba_code"] == ba,
            f"generated_{pollutant}_rate_lb_per_mwh_{rate_type}",
        ].item()
        data_to_graph["egrid_value"] = egrid_value

    return data_to_graph
