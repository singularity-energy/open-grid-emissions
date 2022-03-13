import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from pathlib import Path
from pandas import DataFrame

def clean_eia_930(df:DataFrame):
    """
    Args:
       df (pd.DataFrame): dataframe containing rows of EIA-930 in the format provided by balance
       sheets.
    Returns:
       cleaned df with same format as input
    """
    ## Remove bad data (negative and zero fossil fuel generation)
    fossil_cols = ["Net Generation (MW) from Coal",\
                   "Net Generation (MW) from Natural Gas",\
                   "Net Generation (MW) from All Petroleum Products"]
    for col in fossil_cols:
        df[df[col]<0] = np.nan

    # TODO other forms of cleaning as needed

    return df


def add_report_date(df, plant_entity):
    """
    Add a report date column to the cems data based on the plant's local timezone

    Args:
        df (pd.Dataframe): dataframe containing 'plant_id_eia' and 'operating_datetime_utc' columns
        plant_entity (pd.Dataframe): dataframe containing 'plant_id_eia' and 'timezone'columns
    Returns:
        Original dataframe with 'report_date' column added
    """
    # get timezone
    df = df.merge(
        plant_entity[['plant_id_eia', 'timezone']], how='left', on='plant_id_eia')

    # create a datetimeindex from the operating_datetime_utc column
    datetime_utc = pd.DatetimeIndex(df['operating_datetime_utc'])

    # create blank column to hold local datetimes
    df['report_date'] = np.NaN

    # get list of unique timezones
    timezones = list(df['timezone'].unique())

    # convert timezones to GMT equivalent so that we don't have to deal with daylight saving
    # according to https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
    # however need to do GMT+ rather than GMT- according to https://github.com/pandas-dev/pandas/issues/21509
    tz_to_gmt = {
        'America/New_York': 'Etc/GMT+5',
        'America/Denver': 'Etc/GMT+7',
        'America/North_Dakota/New_Salem': 'Etc/GMT+6',
        'America/Chicago': 'Etc/GMT+6',
        'America/Los_Angeles': 'Etc/GMT+8',
        'America/Indiana/Indianapolis': 'Etc/GMT+5',
        'America/Phoenix': 'Etc/GMT+7',
        'America/Kentucky/Louisville': 'Etc/GMT+5',
        'America/Detroit': 'Etc/GMT+5',
        'America/Indiana/Tell_City': 'Etc/GMT+6',
        'US/Central': 'Etc/GMT+6',
        'America/Anchorage': 'Etc/GMT+9',
        'Pacific/Honolulu': 'Etc/GMT+10',
        'US/Eastern': 'Etc/GMT+5',
        'America/Boise': 'Etc/GMT+7',
        'America/Sitka': 'Etc/GMT+9',
        'US/Pacific': 'Etc/GMT+8',
        'US/Mountain': 'Etc/GMT+7',
        'US/Arizona': 'Etc/GMT+7',
        'America/Juneau': 'Etc/GMT+9',
        'America/Indiana/Vincennes': 'Etc/GMT+5',
        'America/Nome': 'Etc/GMT+9',
        'America/North_Dakota/Beulah': 'Etc/GMT+6',
        'America/Indiana/Petersburg': 'Etc/GMT+5',
        'US/Alaska': 'Etc/GMT+9',
        'US/Hawaii': 'Etc/GMT+10',
        'America/North_Dakota/Center': 'Etc/GMT+6',
        'America/Menominee': 'Etc/GMT+6'
    }

    # convert UTC to the local timezone and change to YYYY-MM format
    for tz in timezones:
        tz_mask = df['timezone'] == tz
        df.loc[tz_mask, 'report_date'] = datetime_utc[tz_mask].tz_convert(
            tz_to_gmt.get(tz, tz)).to_series(index=df[tz_mask].index).astype(str).str[:7]

    # drop the operating_datetime_local column
    df = df.drop(columns=['timezone'])

    return df


def monthly_fuel_types(cems_df, boiler_fuel_eia923, plant_entity):
    """
    Assigns a fuel type to each EPA unitid. This takes a sequential filtering approach to
    associate the fuel type at the boiler level from EIA-923 to each EPA unit, as there is not
    a direct match between EPA units and EIA boilers or generators.

    Args:
        cems_df (DataFrame): a DataFrame selected from the epacems parquet database
            for a single year, with columns including plant_id_eia,	unitid,	operating_datetime_utc,
            operating_time_hours, gross_load_mw, co2_mass_tons, co2_mass_measurement_code, heat_content_mmbtu.
        boiler_fuel_eia923

    TO DO:
    - [ ] check if I can speed up by changing apply to merge
    - [ ] match fuel type from eGRID?
    - [ ] deal with multi-fuels
    """
    # Step 1: Identify which plants have boilers with a single fuel type
    # filter to units that are missing co2 data
    epa_units = epa_unit_list(cems_df, plant_entity)
    primary_fuel = primary_fuel_bf_eia923(
        boiler_fuel_eia923, fuel_thresh=0.9, level='boiler')
    single_fuel_plant = identify_single_fuel_plants(primary_fuel)

    # merge single_fuel_plant data into epa_units
    epa_units = epa_units.merge(single_fuel_plant, how='left', on=[
                                'plant_id_eia', 'report_date'])

    # replace missing values with "unknown"
    epa_units['primary_fuel'] = epa_units['primary_fuel'].fillna('unknown')

    # filter rows where the primary fuel is still missing
    missing = epa_units.query('primary_fuel == "unknown"')

    # Step 2: For the remaining units, check if there is a 1:1 association between boilers and units at a plant
    # For each remaining plant, get a list of boiler unitids and compare that list to boiler ids

    # create a list of all unitids associated with each plant
    unit_list = epa_units[['plant_id_eia', 'unitid']].drop_duplicates()
    unit_list = unit_list.groupby(['plant_id_eia'])['unitid'].apply(list)
    unit_list = unit_list.reset_index()

    # create a list of all boiler_ids associated with each plant

    boiler_list = primary_fuel[['plant_id_eia', 'boiler_id']].drop_duplicates()
    boiler_list = boiler_list.groupby(['plant_id_eia'])[
        'boiler_id'].apply(list)
    boiler_list = boiler_list.reset_index()

    # for each plant in unit_list, check if the list of unitids matches the list of boiler_ids from the corresponding plant in boiler_list
    # create column for match
    unit_list['boiler_match'] = unit_list.apply(
        unit_to_boiler_match, axis=1, args=(boiler_list,))

    # for the units with a 1:1 unitid to boiler_id match, match the boiler primary fuel to the unit primary fuel
    # could I create a mask, then merge the primary fuel?
    # GP: this throws warnings because we assign to a copy, but it's merged later so still works.
    missing['primary_fuel'] = missing.apply(
        assign_boiler_fuel_type, axis=1, args=(unit_list, primary_fuel))

    # update the main dataframe with the filled-in results
    epa_units.update(missing)

    # filter rows where the primary fuel is still missing
    missing = epa_units.query('primary_fuel == "unknown"')

    # load epa-eia crosswalk data
    crosswalk = get_epa_eia_crosswalk()

    # Step 3: Check for a boiler-unit map in the epa-eia crosswalk document
    missing['primary_fuel'] = missing.apply(
        crosswalk_match, axis=1, args=(crosswalk, primary_fuel))

    # update the main dataframe with the filled-in results
    epa_units.update(missing)

    # filter rows where the primary fuel is still missing
    missing = epa_units.query('primary_fuel == "unknown"')

    # Step 4: Pull in fuel types from Crosswalk
    missing['primary_fuel'] = missing.apply(
        fuel_code_lookup, axis=1, level='unit', crosswalk=crosswalk)

    # update the main dataframe with the filled-in results
    epa_units.update(missing)

    # filter rows where the primary fuel is still missing
    missing = epa_units.query('primary_fuel == "unknown"')

    # change plant_id_eia column to integer from float
    epa_units['plant_id_eia'] = epa_units['plant_id_eia'].astype(int)

    return epa_units


def epa_unit_list(cems_df, plant_entity):
    """
    Description.

    Args:
        arg
    Returns:
        output
    """
    # create a dataframe of unique plant_id and unitid combinations from cems
    epa_units = cems_df[['plant_id_eia', 'unitid',
                         'operating_datetime_utc']].drop_duplicates()

    '''
    # Set operating_datetime_utc as a Datetime
    epa_units['operating_datetime_utc'] = pd.to_datetime(
        epa_units['operating_datetime_utc'])
    '''

    # extract month and year data from datetime
    epa_units = add_report_date(epa_units, plant_entity)

    # drop duplicate rows and reset the index
    epa_units = epa_units[['plant_id_eia',
                           'unitid', 'report_date']].drop_duplicates()
    epa_units = epa_units.reset_index(drop=True)

    return epa_units


def primary_fuel_bf_eia923(boiler_fuel_eia923, fuel_thresh, level='boiler'):
    """
    Determines the primary fuel of each boiler or each plant for each month and adds a new primary_fuel column

    Args:
        boiler_fuel_eia923 (pd.DataFrame): a DataFrame of boiler_fuel_eia923 containing the following columns:
            - 'plant_id_eia'
            - 'boiler_id'
            - 'energy_source_code'
            - 'report_date'
            - 'fuel_consumed_units'
            - 'fuel_mmbtu_per_unit'

        fuel_threshold = the threshold percentage above which a fuel is assigned as primary. For example, if set at 0.9, then if a single fuel makes up more than 90% of the heat input, it will be set as primary.

        level (string): specify whether you want boiler fuel proportion at plant level or boiler level. Default is boiler.

    Returns:
        pd.Dataframe: with columns ['plant_id_eia', ['boiler_id',] 'report_date', 'primary_fuel']

    """

    # Figure out the heat content proportions of each fuel received:
    pf_by_heat = fuel_proportions_bf_eia923(boiler_fuel_eia923, level=level)

    if level == 'boiler':
        # On a per boiler, per month basis, identify the fuel that made the largest
        # contribution to the plant's overall heat content consumed. If that
        # proportion is greater than fuel_thresh, set the primary_fuel to be
        # that fuel.  Otherwise, set it to "MULTI".
        pf_by_heat = pf_by_heat.set_index(
            ['plant_id_eia', 'boiler_id', 'report_date'])
    elif level == 'plant':
        # On a per plant, per month basis, identify the fuel that made the largest
        # contribution to the plant's overall heat content consumed. If that
        # proportion is greater than fuel_thresh, set the primary_fuel to be
        # that fuel.  Otherwise, set it to "MULTI".
        pf_by_heat = pf_by_heat.set_index(
            ['plant_id_eia', 'report_date'])

    # identify where percentage is greater than the threshold
    mask = pf_by_heat >= fuel_thresh
    pf_by_heat = pf_by_heat.where(mask)
    # create a new primary fuel column based on the name of the column where the
    pf_by_heat['primary_fuel'] = pf_by_heat.idxmax(axis=1)
    pf_by_heat['primary_fuel'] = pf_by_heat['primary_fuel'].fillna(
        value='unknown')

    return pf_by_heat[['primary_fuel', ]].reset_index()


def fuel_proportions_bf_eia923(boiler_fuel_eia923, level='boiler'):
    """
    Calculates the percentage of heat input from each fuel type for each boiler for each month.

    Args:
        boiler_fuel_eia923 (pd.DataFrame): a DataFrame of boiler_fuel_eia923 containing the following columns:
            - 'plant_id_eia'
            - 'boiler_id'
            - 'energy_source_code'
            - 'report_date'
            - 'fuel_consumed_units'
            - 'fuel_mmbtu_per_unit'
        level (string): specify whether you want boiler fuel proportion at plant level or boiler level. Default is boiler.
    """

    # calculate fuel consumption in mmbtu
    boiler_fuel_eia923['fuel_consumed_mmbtu'] = boiler_fuel_eia923['fuel_consumed_units'] * \
        boiler_fuel_eia923['fuel_mmbtu_per_unit']

    if level == 'boiler':
        # drop fuel_consumed_units and fuel_mmbtu_per_unit columns
        boiler_fuel_eia923 = boiler_fuel_eia923[['report_date',
                                                 'plant_id_eia',
                                                 'boiler_id',
                                                 'energy_source_code',
                                                 'fuel_consumed_mmbtu']]

        # Set report_date as a DatetimeIndex
        boiler_fuel_eia923 = boiler_fuel_eia923.set_index(
            pd.DatetimeIndex(boiler_fuel_eia923['report_date']))

        # Group by report_date (monthly), plant_id_eia, boiler_id, and fuel_type
        bf_gb = boiler_fuel_eia923.groupby(
            ['plant_id_eia', 'boiler_id', pd.Grouper(freq='M'), 'energy_source_code'])

        # Add up all the MMBTU for each boiler & month. At this point each record
        # in the dataframe contains only information about a single fuel.
        heat_df = bf_gb.agg(np.sum)

        # Simplfy the DF a little before we turn it into a pivot table.
        heat_df = heat_df.reset_index()

        # Take the individual rows organized by energy_source_code, and turn them
        # into columns, each with the total MMBTU for that fuel, month, and boiler.
        fuel_pivot = heat_df.pivot_table(
            index=['report_date', 'plant_id_eia', 'boiler_id'],
            columns='energy_source_code',
            values='fuel_consumed_mmbtu')

    elif level == 'plant':
        # drop fuel_consumed_units and fuel_mmbtu_per_unit columns
        boiler_fuel_eia923 = boiler_fuel_eia923[['report_date',
                                                 'plant_id_eia',
                                                 'energy_source_code',
                                                 'fuel_consumed_mmbtu']]

        # Set report_date as a DatetimeIndex
        boiler_fuel_eia923 = boiler_fuel_eia923.set_index(
            pd.DatetimeIndex(boiler_fuel_eia923['report_date']))

        # Group by report_date (monthly), plant_id_eia, and fuel_type
        bf_gb = boiler_fuel_eia923.groupby(
            ['plant_id_eia', pd.Grouper(freq='M'), 'energy_source_code'])

        # Add up all the MMBTU for each boiler & month. At this point each record
        # in the dataframe contains only information about a single fuel.
        heat_df = bf_gb.agg(np.sum)

        # Simplfy the DF a little before we turn it into a pivot table.
        heat_df = heat_df.reset_index()

        # Take the individual rows organized by energy_source_code, and turn them
        # into columns, each with the total MMBTU for that fuel, month, and boiler.
        fuel_pivot = heat_df.pivot_table(
            index=['report_date', 'plant_id_eia'],
            columns='energy_source_code',
            values='fuel_consumed_mmbtu')

    # Add a column that has the *total* heat content of all fuels:
    fuel_pivot['total'] = fuel_pivot.sum(axis=1, numeric_only=True)

    # Replace any NaN values we got from pivoting with zeros.
    fuel_pivot = fuel_pivot.fillna(value=0)

    # drop any months where zero fuel input was recorded
    fuel_pivot = fuel_pivot.drop(fuel_pivot[fuel_pivot['total'] == 0].index)

    # Divide all columns by the total heat content, giving us the proportions
    # for each fuel instead of the heat content.
    fuel_pivot = fuel_pivot.divide(fuel_pivot.total, axis='index')

    # Replace any NaN values we got from dividing by zero with zeros.
    fuel_pivot = fuel_pivot.fillna(value=0)

    # Drop the total column (it's nothing but 1.0 values) and clean up the
    # index and columns a bit before returning the DF.
    fuel_pivot = fuel_pivot.drop('total', axis=1)
    fuel_pivot = fuel_pivot.reset_index()
    fuel_pivot.columns.name = None

    # change the report_date format from YYYY-MM-DD to YYYY-MM
    fuel_pivot['report_date'] = fuel_pivot['report_date'].dt.strftime('%Y-%m')

    return fuel_pivot


def identify_single_fuel_plants(primary_fuel):
    """
    Description.

    Args:
        arg
    Returns:
        output
    """
    # Step 1: figure out which plants use a single fuel type for all boilers
    # Group by plant_id and fuel_type, and check which lists have one unique value
    single_fuel_plant = primary_fuel.groupby(
        ['plant_id_eia', 'report_date', 'primary_fuel']).primary_fuel.nunique()

    # determine where there is a single primary fuel per period
    single_fuel_plant = single_fuel_plant.groupby(
        ['plant_id_eia', 'report_date']).agg(np.sum).eq(1).reset_index()

    # drop rows where primary_fuel is False (there are different fuels being burned in each boiler at a plant)
    single_fuel_plant = single_fuel_plant[single_fuel_plant.primary_fuel == True]

    # replace the primary_fuel column with the actual fuel type
    single_fuel_plant = single_fuel_plant[['plant_id_eia', 'report_date']]
    single_fuel_plant = single_fuel_plant.merge(primary_fuel[[
                                                'plant_id_eia', 'report_date', 'primary_fuel']].drop_duplicates(), how='left', on=['plant_id_eia', 'report_date'])

    # replace missing values with "unknown"
    single_fuel_plant['primary_fuel'] = single_fuel_plant['primary_fuel'].fillna(
        'unknown')

    return single_fuel_plant


def unit_to_boiler_match(row, boiler_list):
    """
    For each row in the unit_list dataframe:
    For each plant_id, checks if the list of unitid's from epacems matches the list of boiler_ids from bf_eia923
    If these two lists match exactly, it is likely that the epa unitid corresponds directly with the eia boiler_id

    ToDo:
        [ ] Pass boiler_list into the function
    """
    # get the plant_id of the current row
    plant_id_eia = row['plant_id_eia']
    # check if the list of unitids and boiler_ids match
    try:
        match = (sorted(row['unitid']) == sorted(
            boiler_list[boiler_list['plant_id_eia'] == plant_id_eia]['boiler_id'].iloc[0]))

    # if there is no match found in the boiler list, return False
    except IndexError:
        match = False

    return match


def assign_boiler_fuel_type(row, unit_list, primary_fuel):
    """
    Assigns a primary fuel type to each CEMS unit.
    Todo:
        [x] Figure out how to pass unit_list and boiler_list into this function
    """
    # get the plant_id_eia
    plant_id_eia = row['plant_id_eia']

    # check if there is a 1:1 unitid to boiler_id match for a given plant
    try:
        if unit_list[unit_list['plant_id_eia'] == plant_id_eia]['boiler_match'].iloc[0] == True:
            # get the unit_id and report_date
            unitid = row['unitid']
            report_date = row['report_date']

            # look up the primary fuel for this boiler_id
            fuel_type = primary_fuel.query(
                'plant_id_eia == @plant_id_eia and boiler_id == @unitid and report_date == @report_date')['primary_fuel'].iloc[0]

        else:
            fuel_type = 'unknown'
    except IndexError:
        fuel_type = 'unknown'
    return fuel_type


def get_epa_eia_crosswalk():
    """
    Read in the manual EPA-EIA Crosswalk table.
    """
    '''
    map_eia_epa_file = importlib.resources.open_binary(
        'pudl.package_data.glue', 'eia_epa_id_crosswalk.csv')

    return pd.read_csv(
        map_eia_epa_file,
        usecols=['plant_id_epa', 'plant_id_eia', 'unitid',
                 'generator_id', 'boiler_id', 'energy_source_code'],
        dtype={'plant_id_epa': 'int32', 'plant_id_eia': 'int32'})'''

    crosswalk = pd.read_csv('../data/epa/epa_eia_crosswalk.csv', usecols=[
                            'CAMD_PLANT_ID', 'EIA_PLANT_ID', 'CAMD_UNIT_ID', 'EIA_GENERATOR_ID', 'EIA_BOILER_ID', 'EIA_FUEL_TYPE', 'CAMD_FUEL_TYPE'])

    # rename the columns
    crosswalk = crosswalk.rename(columns={'CAMD_PLANT_ID': 'plant_id_epa',
                                          'EIA_PLANT_ID': 'plant_id_eia',
                                          'CAMD_UNIT_ID': 'unitid',
                                          'EIA_GENERATOR_ID': 'generator_id',
                                          'EIA_BOILER_ID': 'boiler_id',
                                          'EIA_FUEL_TYPE': 'energy_source_code',
                                          'CAMD_FUEL_TYPE': 'epa_fuel_name'})

    # replace the epa fuel name with the corresponding eia fuel code
    fuel_type_dict = {'Pipeline Natural Gas': 'NG',  # natural gas
                      'Coal': 'CBL',  # Coal, blended
                      'Natural Gas': 'NG',  # natural gas
                      'Other Oil': 'WO',  # Waste/Other Oil.
                      # Residual Fuel Oil. Including No. 5 & 6 fuel oils and bunker C fuel oil.
                      'Residual Oil': 'RFO',
                      # Distillate Fuel Oil. Including diesel, No. 1, No. 2, and No. 4 fuel oils.
                      'Diesel Oil': 'DFO',
                      # Wood/Wood Waste Solids. Including paper pellets, railroad ties, utility polies, wood chips, bark, and other wood waste solids.
                      'Wood': 'WDS',
                      'Other Gas': 'OG',  # Other Gas
                      'Process Gas': 'OG',  # Other Gas
                      'Petroleum Coke': 'PC',  # Petroleum Coke
                      # Waste/Other Coal. Including anthracite culm, bituminous gob, fine coal, lignite waste, waste coal.
                      'Coal Refuse': 'WC',
                      'Other Solid Fuel': 'OBS',  # Other Biomass Solids
                      'Tire Derived Fuel': 'TDF'  # Tire-derived Fuels
                      }
    crosswalk['epa_fuel_code'] = crosswalk['epa_fuel_name'].map(fuel_type_dict)

    # fill missing values in the energy_source_code column with values from the energy_source_code column
    crosswalk['energy_source_code'] = crosswalk['energy_source_code'].fillna(
        crosswalk['energy_source_code'])

    return crosswalk


def crosswalk_match(row, crosswalk, primary_fuel):
    """
    For each row in the epa_units dataframe:
    Checks the EPA-EIA crosswalk table to see if the unitid has been manually matched to a boiler_id or generator_id
    If so, this assigns the fuel type from the primary_fuel table
    """
    # get the plant_id_eia and unit_id
    plant_id_eia = row['plant_id_eia']
    unitid = row['unitid']

    # check if the crosswalk contains a boiler_id for the associated plant_id and unitid
    try:
        boiler_id = crosswalk.query(
            'plant_id_epa == @plant_id_eia and unitid == @unitid')['boiler_id'].iloc[0]

        # if there is no boiler_id matched, do nothing
        if pd.isnull(boiler_id):
            fuel_type = 'unknown'
        else:
            # get the report_date
            report_date = row['report_date']

            # look up the primary_fuel using the boiler_id and report_date
            fuel_type = primary_fuel.query(
                'plant_id_eia == @plant_id_eia and boiler_id == @boiler_id and report_date == @report_date')['primary_fuel'].iloc[0]

    # if the query fails, try to match based on the generator id
    except IndexError:
        try:
            generator_id = crosswalk.query(
                'plant_id_epa == @plant_id_eia and unitid == @unitid')['generator_id'].iloc[0]

            # if there is no generator_id matched, do nothing
            if pd.isnull(generator_id):
                fuel_type = 'unknown'

            # check if the generator_id from the crosswalk matches a boiler_id from primary_fuel
            else:
                # get the report_date
                report_date = row['report_date']

                # look up the primary_fuel using the boiler_id and report_date
                fuel_type = primary_fuel.query(
                    'plant_id_eia == @plant_id_eia and boiler_id == @generator_id and report_date == @report_date')['primary_fuel'].iloc[0]
        # if the generator_id does not match a boiler_id, do nothing
        except IndexError:
            fuel_type = 'unknown'
    return fuel_type


def fuel_code_lookup(row, level, crosswalk):
    """
    Description.

    Args:
        arg
    Returns:
        output
    """
    plant_id_eia = row['plant_id_eia']

    if level == 'unit':
        unitid = row['unitid']

        try:
            fuel_type = crosswalk.query('plant_id_eia == @plant_id_eia and unitid == @unitid')[
                'energy_source_code'].iloc[0]
        except IndexError:
            fuel_type = 'unknown'

    elif level == 'plant':
        try:
            fuel_type = crosswalk.query('plant_id_eia == @plant_id_eia')[
                'energy_source_code'].iloc[0]
        except IndexError:
            fuel_type = ''

    return fuel_type

def fill_missing_co2(cems_df):
    """
    This function uses fuel type data filled into cems_df by the monthly_fuel_types()
    function to fill in co2 data where missing.
    To do this, it multiplies the heat_input_mmbtu by the emission factor for that fuel type
    """

    missing = cems_df[cems_df['co2_mass_tons'].isnull()]
    # get emission factors
    fuel_ef_per_mmbtu = get_emissions_factors(
    )[['energy_source_code', 'co2_tons_per_mmbtu']]
    # add emission factor to missing df
    missing = missing.merge(fuel_ef_per_mmbtu, how='left',
                            left_on='primary_fuel', right_on='energy_source_code')
    # calculate missing co2 data
    missing['co2_mass_tons'] = missing['heat_content_mmbtu'] * \
        missing['co2_tons_per_mmbtu']

    return missing

def get_emissions_factors():
    """
    Read in the table of emissions factors
    TODO: add year as an input to get the emissions factors from the correct egrid year
    """
    year = 2019
    return pd.read_csv(f'../data/egrid/egrid{year}_static_tables/table_C1_emission_factors_for_CO2_CH4_N2O.csv')

def calculate_heat_input_weighted_ef(boiler_fuel_eia923, level):
    # Figure out the heat content proportions of each fuel received:
    fp_by_heat = fuel_proportions_bf_eia923(boiler_fuel_eia923, level=level)
    if level == 'boiler':
        fp_by_heat = fp_by_heat.set_index(
            ['plant_id_eia', 'boiler_id', 'report_date'])
    elif level == 'plant':
        fp_by_heat = fp_by_heat.set_index(
            ['plant_id_eia', 'report_date'])

    #create a weighted emission factor
    fuel_ef_per_mmbtu = get_emissions_factors(
    )[['fuel_code', 'co2_tons_per_mmbtu']]

    #get list of unique column names
    fuel_list = list(fp_by_heat.columns.unique())

    #multiply the percentage in each column by the associated fuel emission factor
    for fuel in fuel_list:
        try:
            fp_by_heat[fuel] = fp_by_heat[fuel] * fuel_ef_per_mmbtu.loc[fuel_ef_per_mmbtu['fuel_code'] == fuel, 'co2_tons_per_mmbtu'].to_numpy()[0]
        except IndexError:
            fp_by_heat[fuel] = fp_by_heat[fuel] * 0

    #calculate the weighted EF by summing across columns
    fp_by_heat['fuel_weighted_ef_tons_per_mmbtu'] = fp_by_heat.sum(axis=1)

    weighted_ef = fp_by_heat[['fuel_weighted_ef_tons_per_mmbtu']]

    return weighted_ef
