# Repo Structure
hourly-egrid  
|- data: directory for local versions of downloaded files  
|   |- egrid  
|   |   |- egrid2019_static_tabes  
|- notebooks: jupyter notebooks for testing and running code  
|   |- data_pipeline.ipynb  
|- src: python code called by notebooks  
|   |- data_cleaning.py  
|   |- gross_to_net_generation.py  
|- test  


# Installation

## Pipenv
This package uses pipenv and virtualenvwrapper to manage python packages.
Once you have installed both, you can set up the packages required for this code
by running `pipenv install` from the top-level project directory.
To also include packages used in testing and development, run `pipenv install --dev`

## Conda

Open anaconda prompt, use cd to navigate to the directory where your local files are stored (e.g. "GitHub/hourly-egrid"), and then run:

```
conda env create -f environment.yml
```

# Running

EIA data access looks for an API key in environment variable EIA_API_KEY.
You can sign up for an API key at https://www.eia.gov/opendata/register.php


# Testing

To test, run `pytest` from the top-level project directory. 
