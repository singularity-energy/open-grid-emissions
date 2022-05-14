import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from pathlib import Path
import src.data_cleaning as data_cleaning
import sqlalchemy as sa

import pudl.analysis.epa_crosswalk as epa_crosswalk
import pudl.output.pudltabl


def load_cems_gross_generation(start_year, end_year):
    """Loads hourly CEMS gross generation data for multiple years."""
    cems_all = []

    for year in range(start_year, end_year + 1):
        print(f'loading {year}')
        # specify the path to the CEMS data
        cems_path = f"../data/pudl/pudl_data/parquet/epacems/year={year}"

        # specify the columns to use from the CEMS database
        cems_columns = [
            "plant_id_eia",
            "unitid",
            "operating_datetime_utc",
            "operating_time_hours",
            "gross_load_mw"]

        # load the CEMS data
        cems = pd.read_parquet(cems_path, columns=cems_columns)

        # rename cems plant_id_eia to plant_id_epa (PUDL simply renames the ORISPL_CODE column from the raw CEMS data as 'plant_id_eia' without actually crosswalking to the EIA id)
        # rename the heat content column to use the convention used in the EIA data
        cems = cems.rename(columns={"plant_id_eia": "plant_id_epa",})

        # if the unitid has any leading zeros, remove them
        cems["unitid"] = cems["unitid"].str.lstrip("0")

        # crosswalk the plant IDs and add a plant_id_eia column
        cems = data_cleaning.crosswalk_epa_eia_plant_ids(cems, year)

        # fill any missing values for operating time or steam load with zero
        cems["operating_time_hours"] = cems["operating_time_hours"].fillna(0)

        # calculate gross generation by multiplying gross_load_mw by operating_time_hours
        cems["gross_generation_mwh"] = cems["gross_load_mw"] * cems["operating_time_hours"]

        cems_all.append(cems)

    cems = pd.concat(cems_all, axis=0)

    # separate out data that only has zeros reported
    # this will help speed up calculations and allow us to add this data back later
    #cems_zero_hours = cems[cems['gross_load_mw'] == 0]
    #cems = cems[cems['gross_load_mw'] > 0]

    # add a report date
    print('adding report date')
    cems = data_cleaning.add_report_date(cems)

    return cems  #, cems_zero_hours

def identify_subplants(start_year, end_year, gen_fuel_allocated):
    """
    Groups units and generators into unique subplant groups.

    This function consists of three primary parts:
    1. Identify a list of all unique plant-units that exist in the CEMS data
        for the years in question. This will be used to filter the crosswalk.
    2. Load the EPA-EIA crosswalk and filter it based on the units that exist
        in the CEMS data for the years in question
    3. Use graph analysis to identify distinct groupings of EPA units and EIA
        generators based on 1:1, 1:m, m:1, or m:m relationships.
    
    """

    # get a list of unique keys from epacems
    ids_all = []
    for year in range(start_year, end_year + 1):
        # specify the path to the CEMS data
        cems_path = f"../data/pudl/pudl_data/parquet/epacems/year={year}"

        # specify the columns to use from the CEMS database
        cems_columns = [
            "plant_id_eia",
            "unitid",
            'unit_id_epa']

        # load the CEMS data
        ids = pd.read_parquet(cems_path, columns=cems_columns)

        ids = ids[["plant_id_eia", "unitid", "unit_id_epa"]].drop_duplicates()

        ids_all.append(ids)

    ids = pd.concat(ids_all, axis=0)

    ids = ids[["plant_id_eia", "unitid", "unit_id_epa"]].drop_duplicates()

    # load the crosswalk and filter it by the data that actually exists in cems
    crosswalk = pudl.output.epacems.epa_crosswalk()
    filtered_crosswalk = epa_crosswalk.filter_crosswalk(crosswalk, ids)[['plant_id_eia','unitid','unit_id_epa','CAMD_PLANT_ID','CAMD_UNIT_ID','CAMD_GENERATOR_ID','EIA_PLANT_ID','EIA_GENERATOR_ID']]

    # change the plant id to an int
    filtered_crosswalk['EIA_PLANT_ID'] = filtered_crosswalk['EIA_PLANT_ID'].astype(int)

    # filter out generators that are retired
    """
    # load the generator retirement date
    pudl_db = "sqlite:///../data/pudl/pudl_data/sqlite/pudl.sqlite"
    pudl_engine = sa.create_engine(pudl_db)
    pudl_out = pudl.output.pudltabl.PudlTabl(pudl_engine, freq="MS", start_date=f"{start_year}-01-01", end_date=f"{end_year}-12-31")
    generator_retirement_date =  pudl_out.gens_eia860().loc[:,['plant_id_eia','generator_id','retirement_date']]
    # add the retirement date to the crosswalk
    filtered_crosswalk = filtered_crosswalk.merge(generator_retirement_date, how='left', left_on=['EIA_PLANT_ID','EIA_GENERATOR_ID'], right_on=['plant_id_eia','generator_id'])
    filtered_crosswalk = filtered_crosswalk.drop(columns=['plant_id_eia_y','generator_id']).rename(columns={'plant_id_eia_x':'plant_id_eia'})
    # only keep data for plants that have not already retired before the start year
    filtered_crosswalk = filtered_crosswalk[~(filtered_crosswalk['retirement_date'].dt.year < start_year)]
    """
    # filter to generators that exist in the EIA data
    # get a list of unique generators in the EIA-923 data
    unique_eia_ids = gen_fuel_allocated[['plant_id_eia','generator_id']].drop_duplicates()
    filtered_crosswalk = unique_eia_ids.merge(
            filtered_crosswalk,
            left_on=["plant_id_eia", "generator_id"],
            right_on=["EIA_PLANT_ID", "EIA_GENERATOR_ID"],
            how="inner",
            suffixes=("_actual",None)
        ).drop(columns=['plant_id_eia_actual','generator_id'])

    crosswalk_with_subplant_ids = epa_crosswalk.make_subplant_ids(filtered_crosswalk)
    # fix the column names
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.drop(columns=['plant_id_eia','unitid','unit_id_epa','CAMD_GENERATOR_ID'])
    crosswalk_with_subplant_ids = crosswalk_with_subplant_ids.rename(
            columns={
                "CAMD_PLANT_ID": "plant_id_epa",
                "EIA_PLANT_ID": "plant_id_eia",
                "CAMD_UNIT_ID": "unitid",
                "EIA_GENERATOR_ID": "generator_id",
            }
        )
    # change the eia plant id to int
    crosswalk_with_subplant_ids['plant_id_eia'] = crosswalk_with_subplant_ids['plant_id_eia'].astype(int)

    return crosswalk_with_subplant_ids

def gross_to_net_ratios(cems_df, generators, plant_entity):
    """
    Calculates a monthly gross to net generation ratio (net generation / gross generation) that can be used to convert gross generation to net generation

    Args:
        arg
    Returns:
        output
    """

    # add a monthly report data to the cems data and get a list of unit_ids that reported to CEMS in each month
    reporting_units = data_cleaning.add_report_date(cems_df, plant_entity)[["report_date", "plant_id_eia", "unitid"]].drop_duplicates()

    #load the EPA-EIA crosswalk data
    crosswalk = load_data.load_epa_eia_crosswalk()[['plant_id_eia','generator_id', 'unitid']]

    # merge the generator_id into the reporting units list, keeping duplicate matches if multiple generators are associated with a single unit
    reporting_units = reporting_units.merge(crosswalk, how='left', on=['plant_id_eia', 'unitid'])

    # only keep the columns we need for grouping the data
    cems_df = cems_df[['plant_id_eia',
                       'operating_datetime_utc', 'gross_generation_mwh']]

    # set the datetimeindex
    cems_df = cems_df.set_index(pd.DatetimeIndex(
        cems_df['operating_datetime_utc']))

    # create a report date column to match the eia_923 data
    cems_df = data_cleaning.add_report_date(cems_df, plant_entity)

    # sum the cems data by plant and month
    cems_df = cems_df.groupby(
        ['plant_id_eia', 'report_date']).sum().reset_index() #['plant_id_eia', pd.Grouper(freq='M')]).sum().reset_index()

    # re-format the report date column in EIA-923 from YYYY-MM-DD to YYYY-MM
    generators['report_date'] = (pd.to_datetime(
        generators['report_date'], format='%Y-%m-%d')).dt.strftime('%Y-%m')

    #drop any observations from EIA-923 that are not in ["report_date", "plant_id_eia", "generator_id"] from reporting_units
    generators = generators.merge(reporting_units, how='inner', on=["report_date", "plant_id_eia", "generator_id"])
    generators= generators.drop_duplicates(subset=["report_date", "plant_id_eia", "generator_id"])

    # aggregate EIA 293 data to plant level
    generators = generators[['plant_id_eia', 'report_date', 'net_generation_mwh']].groupby(
        ['plant_id_eia', 'report_date']).sum().reset_index()

    # merge the net generation data into the cems data
    cems_df = cems_df.merge(generators, how='left', on=[
                            'plant_id_eia', 'report_date'])
    cems_df['plant_id_eia'] = cems_df['plant_id_eia'].astype(int)

    # create a new dataframe to hold the gross-to-net ratios, and delete any observations with missing data
    gtn_regression = cems_df.copy().dropna(how='any', axis=0)
    gtn_regression = gtn_regression.groupby(
        'plant_id_eia').apply(model_gross_to_net)
    gtn_regression = pd.DataFrame(gtn_regression.tolist(
    ), index=gtn_regression.index, columns=['gtn_linear', 'rsquared']).reset_index()

    # only keep the results with rsquared values greater than 0.90
    gtn_regression = gtn_regression[gtn_regression['rsquared'] >= 0.9]
    gtn_regression = gtn_regression[['plant_id_eia', 'gtn_linear']]

    # get month-by-month gross-to-net ratio
    cems_df['gtn_ratio'] = cems_df['net_generation_mwh'] / \
        cems_df['gross_generation_mwh']

    # trim values that are outliers
    cems_df.loc[cems_df['gtn_ratio'] > 1.5, 'gtn_ratio'] = np.NaN
    cems_df.loc[cems_df['gtn_ratio'] < 0.5, 'gtn_ratio'] = np.NaN

    # merge in regression results
    cems_df = cems_df.merge(gtn_regression, how='left', on='plant_id_eia')

    # fill missing values using the ratio from the regression results
    cems_df['gtn_ratio'] = cems_df['gtn_ratio'].fillna(cems_df['gtn_linear'])

    gtn_ratios = cems_df[['plant_id_eia', 'report_date', 'gtn_ratio']]

    gtn_fill_values = gtn_ratios[['plant_id_eia', 'gtn_ratio']].groupby(['plant_id_eia']).mean().rename(columns={'gtn_ratio':'gtn_fill'})

    return gtn_ratios, gtn_fill_values

def model_gross_to_net(df):
    """
    Description

    Args:
        arg
    Returns:
        output
    """
    # get a linear model for the data points
    model = smf.ols('net_generation_mwh ~ gross_generation_mwh', data=df).fit()

    # find and remove any outliers
    try:
        outliers = model.outlier_test()
        corrected = df[~df.index.isin(
            outliers[outliers['bonf(p)'] < 0.5].index)]

        # get a linear model of the corrected data
        model = smf.ols(
            'net_generation_mwh ~ gross_generation_mwh', data=corrected).fit()
    except ValueError:
        pass
    slope = model.params[1]
    rsquared = model.rsquared

    return slope, rsquared
