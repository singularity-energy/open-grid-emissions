import pandas as pd

# Map from 923 fuel types (added to cems data in data_pipeline) to 930 fuel types
# Based on EIA 923 instructions: https://www.eia.gov/survey/form/eia_923/instructions.pdf
fuel_code_map = {'NG':"NG",\
     'BIT':"COL", # bituminous
     'SUB':"COL", # subbituminous
     'BLQ':"OTH", # black liquor
     'RC':"COL", # refined coal
     'WDS':"OTH", # wood / wood waste solids
     'SUN':"SUN", # solar. TODO: why is this in cems??? 
     'KER':"OIL", # kerosene
     'RFO':"OIL", # residual fuel oil
     'JF':"OIL", # jet fuel
     'DFO':"OIL", # distillate fuel oil
     'SGC':"NG", # synthetic gas derived from coal ???
     'LIG':"COL", #lignite coal
     'OG':"NG", # other gas 
     'PC':"COL", # petroleum coke
     'BFG':"NG", # blast furnace gas 
     'WC':"COL"} # waste coal

"""
    reformat_chalendar
Reformat wide-format data (one row per time stamp) from Chalendar to long (one row per data point)
Drop all columns that are not fuel-specific generation 
"""
def reformat_chalendar(raw):
    target_cols = [c for c in raw.columns if len(c.split("."))==5] # where we have variable (NG = net generation) and fuel type
    print("Filtering")
    cleaned = raw.loc[:,target_cols].melt(ignore_index=False, value_name="generation").reset_index()
    print("Expanding cols")
    cleaned[["dtype","BA", "other BA", "var","fuel","interval"]] = cleaned["variable"].str.split(r"[.-]", expand=True, regex=True)
    print("Dropping and renaming")
    cleaned = cleaned.drop(columns=["dtype","var","interval", "other BA"])
    cleaned = cleaned.rename(columns={"index":"datetime_utc"})
    # Correct to start-of-hour
    #cleaned.loc[~cleaned.BA.isin(START_OF_HOUR_REGIONS), "datetime_utc"] += dt.timedelta(hours=-1)
    return cleaned 

def load_chalendar(fname:str, year:int=2020):
    raw = pd.read_csv(fname,index_col=0, parse_dates=True)
    raw = raw[raw.index.year==year]
    return reformat_chalendar(raw)