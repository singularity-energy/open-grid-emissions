import pandas as pd
import src.data_cleaning as data_cleaning

# Map from 923 fuel types (added to cems data in data_pipeline) to 930 fuel types
# Based on EIA 923 instructions: https://www.eia.gov/survey/form/eia_923/instructions.pdf
fuel_code_map = {'NG': "NG",
                 'BIT': "COL",  # bituminous
                 'SUB': "COL",  # subbituminous
                 'BLQ': "OTH",  # black liquor
                 'RC': "COL",  # refined coal
                 'WDS': "OTH",  # wood / wood waste solids
                 'SUN': "SUN",  # solar. TODO: why is this in cems???
                 'KER': "OIL",  # kerosene
                 'RFO': "OIL",  # residual fuel oil
                 'JF': "OIL",  # jet fuel
                 'DFO': "OIL",  # distillate fuel oil
                 'SGC': "NG",  # synthetic gas derived from coal ???
                 'LIG': "COL",  # lignite coal
                 'OG': "NG",  # other gas
                 'PC': "COL",  # petroleum coke
                 'BFG': "NG",  # blast furnace gas
                 'WC': "COL"}  # waste coal

"""
    reformat_chalendar
Reformat wide-format data (one row per time stamp) from Chalendar to long (one row per data point)
Drop all columns that are not fuel-specific generation 
"""


def reformat_chalendar(raw):
    # where we have variable (NG = net generation) and fuel type
    target_cols = [c for c in raw.columns if len(c.split(".")) == 5]
    print("Filtering")
    cleaned = raw.loc[:, target_cols].melt(
        ignore_index=False, value_name="generation").reset_index()
    print("Expanding cols")
    cleaned[["dtype", "BA", "other BA", "var", "fuel", "interval"]
            ] = cleaned["variable"].str.split(r"[.-]", expand=True, regex=True)
    print("Dropping and renaming")
    cleaned = cleaned.drop(columns=["dtype", "var", "interval", "other BA"])
    cleaned = cleaned.rename(columns={"index": "datetime_utc"})
    # Correct to start-of-hour
    #cleaned.loc[~cleaned.BA.isin(START_OF_HOUR_REGIONS), "datetime_utc"] += dt.timedelta(hours=-1)
    return cleaned


def load_chalendar(fname: str, year: int = 2020):
    raw = pd.read_csv(fname, index_col=0, parse_dates=True)
    raw = raw[raw.index.year == year]
    return reformat_chalendar(raw)


def load_chalendar_for_pipeline(year):
    """
    Loads and formats cleaned hourly net generation data for use in the data pipeline
    """
    # read the data, only keeping net generation columns
    data = pd.read_csv('../data/eia930/chalendar/EBA_elec.csv',
                       index_col=0, parse_dates=True).filter(like='-ALL.NG.')

    # name the index
    data.index = data.index.rename('datetime_utc')

    # only keep data for the single year we are interested, plus or minus one day, to account for conversion to local time later
    data = data.loc[(data.index >= f'{year-1}-12-31')
                    & (data.index < f'{year+1}-01-02'), :]

    # remove columns for total net generation
    data = data.loc[:, [col for col in data.columns if col not in data.filter(
        like='-ALL.NG.H').columns]]

    # melt the data into long format
    data = data.reset_index().melt(id_vars='datetime_utc',
                                   value_name='net_generation_mwh_930')

    # create columns for ba_code and fuel category
    data[['ba_code', 'fuel_category']] = data["variable"].str.split(
        r"[.-]", expand=True, regex=True)[[1, 4]]

    # drop BAs not located in the United States
    foreign_bas = ['AESO', 'BCHA', 'CEN',
                   'CFE', 'HQT', 'IESO', 'MHEB', 'SPC', ]
    data = data[~data['ba_code'].isin(foreign_bas)]

    # create a local datetime column
    # TODO: Convert to local prevailing time, not local standard time
    data['datetime_local'] = data['datetime_utc']
    for ba in list(data['ba_code'].unique()):
        data.loc[data.ba_code == ba, 'datetime_local'] = data.loc[data.ba_code ==
                                                                  ba, 'datetime_utc'].dt.tz_convert(data_cleaning.ba_timezone(ba, 'GMT'))

    # create a report date column
    data['report_date'] = data['datetime_local'].astype(str).str[:7]
    data['report_date'] = pd.to_datetime(data['report_date'])

    # rename the fuel categories using format in data/manual/energy_source_groups
    fuel_categories = {'COL':'coal', 
                    'NG':'natural_gas', 
                    'OTH':'other', 
                    'WAT':'hydro', 
                    'WND':'wind', 
                    'SUN':'solar', 
                    'NUC':'nuclear', 
                    'OIL':'petroleum', 
                    'BIO':'biomass',
                    'GEO':'geothermal'}

    data['fuel_category'] = data['fuel_category'].map(fuel_categories)

    # reorder the columns and remove the variable column
    data = data[['ba_code', 'fuel_category', 'datetime_utc',
                 'datetime_local', 'report_date', 'net_generation_mwh_930']]

    return data
