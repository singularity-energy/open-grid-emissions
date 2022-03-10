import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from pathlib import Path
import src.data_cleaning as data_cleaning

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
    crosswalk = data_cleaning.get_epa_eia_crosswalk()[['plant_id_eia','generator_id', 'unitid']]

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
