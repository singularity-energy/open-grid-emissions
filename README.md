# Hourly eGRID

## Modules
- `data_pipeline`: main script for running the data pipeline from start to finish
- `download_data`: functions that download data from the internet
- `load_data`: functions for loading data from downloaded files
- `data_cleaning`: functions that clean loaded data
- `gross_to_net_generation`: Functions for identifying subplants and gross to net generation conversion factors
- `eia930`: functions for cleaning and formatting EIA-930 data
- `impute_hourly_profiles`: functions related to assigning an hourly profile to monthly data
- `output_data`: functions for writing intermediate and final data to csvs
- `column_checks`: functions that check that all data outputs have the correct column names
- `validation`: functions for testing and validating data outputs
- `visualization`: functions for visualizing data in notebooks

## Notebooks
- `data_validation` is for running validatation tests on data

## Data Structure
- `data/downloads` contains all files that are downloaded by functions in `load_data`
- `data/manual` contains all manually-created files, including the egrid static tables
- `data/outputs` contains intermediate outputs from the data pipeline... any files created by our code that are not final results
- `data/results` contains all final output files that will be published

## Installation

Python packages may be managed using `conda`

### Conda

Open anaconda prompt, use `cd` to navigate to the directory where your local files are stored (e.g. "GitHub/hourly-egrid"), and then run:

```
conda update conda
conda env create -f environment.yml
```

Installation requires that the conda channel-priority be set to "flexible". This is the default behavior, 
so if you've never manually changed this, you shouldn't have to worry about this. However, 
if you receive an error message like "Found conflicts!" when trying to install the environment,
try setting your channel priority to flexible by running the following command:
`conda config --set channel_priority flexible` and then re-running the above commands.
