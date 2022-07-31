import pandas as pd


def calculate_ba_gross_loss(year):
    """Calculates the grid gross loss value for each Balancing Area.

    Grid Gross Loss = Total Energy Losses / Total Disposition

    Losses and disposition are reported for each utility-state. Thus, we need to map this
    disposition data to each BA. However, each utility may sell energy into multiple BAs,
    so we need a way to allocate utility disposition to each BA.

    In a separate table in EIA-861, utilities report their total retail sales to customers in each BA.
    We can use this data to calculate what fraction of a utility's total retail sales are in each BA,
    and use this fraction of sales to allocate total disposition and losses from each utility to each BA.

    Once we allocate losses and disposition to each BA, we can sum these values by BA and then
    calculate GGL using the formula above.
    """

    # calculate total disposition by utility
    utility_disposition = calculate_utility_disposition(year)

    # calculate the allocation fraction based on retail sales
    retail_sales_by_ba = calculate_retail_sales_fraction(year)

    # allocate utility disposition to each ba
    ba_disposition = allocate_disposition_to_ba(
        retail_sales_by_ba, utility_disposition, year
    )

    # calculate grid gross loss
    ba_disposition["grid_gross_loss"] = ba_disposition["total_energy_losses_mwh"] / (
        ba_disposition["total_disposition_mwh"] - ba_disposition["sales_for_resale_mwh"]
    )
    # fill missing ggl due to divide by zero with zero
    ba_disposition["grid_gross_loss"] = ba_disposition["grid_gross_loss"].fillna(0)

    # check that the grid gross loss values are reasonable
    # we flag negative values and values greater than 20%
    bad_ggl = ba_disposition[
        (ba_disposition["grid_gross_loss"] < 0)
        | (ba_disposition["grid_gross_loss"] > 0.2)
    ]
    if len(bad_ggl) > 0:
        print("Warning: Certain BAs have abnormal grid gross loss values:")
        print(bad_ggl)

    return ba_disposition


def calculate_utility_disposition(year):
    utility_op_data = load_utility_operational_data(year)

    # sum relevant disposition data by utility
    utility_disposition = utility_op_data.copy()[
        [
            "utility_id_eia",
            "state",
            "nerc_region_code",
            "sales_for_resale_mwh",
            "total_energy_losses_mwh",
            "total_disposition_mwh",
        ]
    ]
    utility_disposition.loc[
        :, ["sales_for_resale_mwh", "total_energy_losses_mwh", "total_disposition_mwh"]
    ] = utility_disposition.loc[
        :, ["sales_for_resale_mwh", "total_energy_losses_mwh", "total_disposition_mwh"]
    ].fillna(
        0
    )
    utility_disposition = (
        utility_disposition.groupby(["utility_id_eia"], dropna=False)[
            ["sales_for_resale_mwh", "total_energy_losses_mwh", "total_disposition_mwh"]
        ]
        .sum()
        .reset_index()
    )
    return utility_disposition


def calculate_retail_sales_fraction(year):
    # figure out what fraction of each utility's retail sales occur in each BA
    retail_sales = load_utility_retail_sales_data(year)

    # replace all zeros with a small number to prevent divide by zero errors
    retail_sales["retail_sales_mwh"] = retail_sales["retail_sales_mwh"].replace(
        0, 0.0001
    )
    # calculate what fraction of a utility's total retail sales are in each BA
    retail_sales_by_ba = (
        retail_sales.groupby(["utility_id_eia", "ba_code"], dropna=False)[
            "retail_sales_mwh"
        ].sum()
        / retail_sales.groupby("utility_id_eia", dropna=False)["retail_sales_mwh"].sum()
    ).round(6)
    # clean up the dataframe
    retail_sales_by_ba = retail_sales_by_ba.reset_index().dropna(thresh=2)
    retail_sales_by_ba = retail_sales_by_ba.rename(
        columns={"retail_sales_mwh": "fraction_of_sales"}
    )
    return retail_sales_by_ba


def allocate_disposition_to_ba(retail_sales_by_ba, utility_disposition, year):

    # merge disposition into the allocation fraction dataframe
    ba_disposition = retail_sales_by_ba.merge(
        utility_disposition,
        how="outer",
        on="utility_id_eia",
        validate="m:1",
        indicator="source",
    )
    # if there is no retail sales data reported by a ba, set allocation fraction to 1
    ba_disposition.loc[
        ba_disposition["source"] == "right_only", "fraction_of_sales"
    ] = ba_disposition.loc[
        ba_disposition["source"] == "right_only", "fraction_of_sales"
    ].fillna(
        1
    )

    # for these bas with missing fractions, they will also be missing ba codes
    # let's try and fill any missing ba codes that we can using data reported in a different table
    utility_rto = load_utility_rto(year)
    ba_disposition = ba_disposition.merge(
        utility_rto, how="left", on="utility_id_eia", validate="m:1"
    )
    ba_disposition["ba_code"] = ba_disposition["ba_code"].fillna(
        ba_disposition["rto_code"]
    )

    # allocate disposition to each BA
    for col in [
        "sales_for_resale_mwh",
        "total_energy_losses_mwh",
        "total_disposition_mwh",
    ]:
        ba_disposition[col] = ba_disposition[col] * ba_disposition["fraction_of_sales"]

    # total the values by BA
    ba_disposition = (
        ba_disposition.groupby("ba_code", dropna=False)[
            ["sales_for_resale_mwh", "total_energy_losses_mwh", "total_disposition_mwh"]
        ]
        .sum()
        .reset_index()
    )

    return ba_disposition


def load_utility_operational_data(year):
    """Loads data from the EIA-861 Operational Data spreadsheet.

    Note: this is temporary until this form is integrated into PUDL
    """
    column_names = [
        "report_date",
        "utility_id_eia",
        "utility_name_eia",
        "state",
        "utility_ownership_type",
        "nerc_region_code",
        "peak_demand_summer",
        "peak_demand_winter",
        "net_generation_mwh",
        "wholesale_power_purchases_mwh",
        "exchange_energy_received_mwh",
        "exchange_energy_delivered_mwh",
        "net_energy_exchanged_mwh",
        "wheeled_power_received_mwh",
        "wheeled_power_delivered_mwh",
        "net_wheeled_power_mwh",
        "transmission_losses_by_others_mwh",
        "total_sources_mwh",
        "retail_sales_mwh",
        "sales_for_resale_mwh",
        "furnished_without_charge_mwh",
        "self_consumed_mwh",
        "total_energy_losses_mwh",
        "total_disposition_mwh",
    ]

    utility_op_data = pd.read_excel(
        f"../data/downloads/eia861/f861{year}/Operational_Data_{year}.xlsx",
        sheet_name="States",
        header=2,
        na_values=".",
        names=column_names,
        usecols="A:X",
    )

    return utility_op_data


def load_utility_retail_sales_data(year):
    """Loads the two retail sales tables from EIA-861 and concats them together."""

    column_names = ["utility_id_eia", "state", "ba_code", "retail_sales_mwh"]

    retail_sales = pd.read_excel(
        f"../data/downloads/eia861/f861{year}/Sales_Ult_Cust_{year}.xlsx",
        sheet_name="States",
        header=2,
        na_values=".",
        names=column_names,
        usecols="B,G,I,W",
    )
    retail_sales_cs = pd.read_excel(
        f"../data/downloads/eia861/f861{year}/Sales_Ult_Cust_CS_{year}.xlsx",
        sheet_name="Customer Sited",
        header=2,
        na_values=".",
        names=column_names,
        usecols="D,F,G,U",
    )

    # concat the two dataframes toegether
    retail_sales = pd.concat([retail_sales, retail_sales_cs], axis=0)

    # create a missing ba code for any records where a ba code is missins
    retail_sales["ba_code"] = retail_sales["ba_code"].fillna(
        retail_sales["state"] + "MS"
    )

    retail_sales["retail_sales_mwh"] = retail_sales["retail_sales_mwh"].fillna(0)

    return retail_sales


def load_utility_rto(year):
    """Loads utilities that report belonging to a single RTO to help fill missing ba codes."""

    utility_rto = (
        pd.read_excel(
            f"../data/downloads/eia861/f861{year}/Utility_Data_{year}.xlsx",
            sheet_name="States",
            header=1,
            na_values=".",
            names=[
                "utility_id_eia",
                "state",
                "CISO",
                "ERCO",
                "PJM",
                "NYIS",
                "SWPP",
                "MISO",
                "ISNE",
            ],
            usecols="B,D,O:U",
        )
        .replace("Y", 1)  # replace Y/N with 1/0
        .replace("N", 0)
        .fillna(0)
    )

    # only keep data where there is a single RTO identified
    utility_rto = utility_rto[
        utility_rto[["CISO", "ERCO", "PJM", "NYIS", "SWPP", "MISO", "ISNE"]].sum(axis=1)
        == 1
    ]

    # melt the data to get the ba_code for each utility
    utility_rto = utility_rto.melt(
        id_vars=["utility_id_eia", "state"], var_name="rto_code"
    )
    utility_rto = utility_rto[utility_rto["value"] == 1].drop(
        columns=["value", "state"]
    )
    utility_rto = utility_rto.drop_duplicates()

    return utility_rto
