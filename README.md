# Installation

This package uses pipenv and virtualenvwrapper to manage python packages.
Once you have installed both, you can set up the packages required for this code
by running `pipenv install` from the top-level project directory.
To also include packages used in testing and development, run `pipenv install --dev`

# Running

EIA data access looks for an API key in environment variable EIA_API_KEY.
You can sign up for an API key at https://www.eia.gov/opendata/register.php

# Testing

To test, run `pytest` from the top-level project directory. 
