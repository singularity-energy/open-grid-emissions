# Open Grid Emissions Initiative
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7062459.svg)](https://doi.org/10.5281/zenodo.7062459)

The Open Grid Emissions Initiative seeks to fill a critical need for high-quality, publicly-accessible, hourly grid emissions data that can be used for GHG accounting, policymaking, academic research, and energy attribute certificate markets. The initiative includes this repository of open-source grid emissions data processing tools that use peer-reviewed, well-documented, and validated methodologies to create the accompanying public dataset of hourly, monthly, and annual U.S. electric grid generation, GHG, and air pollution data.

Please check out [our documentation](https://docs.singularity.energy/docs/open-grid-emissions) for more details about the Open Grid Emissions methodology.

The Open Grid Emissions Dataset can be [downloaded here](https://singularity.energy/open-grid-emissions). An archive of previous versions of the dataset and intermediate data outputs (for research and validation purposes) can be found on [Zenodo](https://zenodo.org/communities/singularity-energy?page=1&size=20).

## Installing and running the data pipeline
To manage the code environment necessary to run the OGE data pipeline, either `pipenv` or `conda` may be used. Currently, we utilize `pipenv` as our preferred environment manager for running the pipeline that is used for data releases, but `conda` will also work if you are more familiar with `conda`. 

First, navigate to the folder where you want to save the repository and run the following commands:

### If you are using pipenv
Note that this option requires to have Python and git installed on your machine.
```bash
pip install pipenv
git clone https://github.com/singularity-energy/open-grid-emissions.git
cd open-grid-emissions
pipenv sync
pipenv shell
pip install build
python -m build
pip install .
```

### If you are using conda
```bash
conda install git
git clone https://github.com/singularity-energy/open-grid-emissions.git
conda update conda
cd open-grid-emissions
conda env create -f environment.yml
conda activate open_grid_emissions
pip install build
python -m build
pip install .
```

The pipeline can be run as follows:
```bash
cd src/oge
python data_pipeline.py --year 2022
```
independently of the installation method you chose.

A more detailed walkthrough of these steps can be found below in the "Development Setup" section.

## Data Availability and Release Schedule
The latest release includes data for year 2005-2022 covering the contiguous United States, Alaska, and Hawaii. In future releases, we plan to expand the geographic coverage to additional U.S. territories (dependent on data availability).

Parts of the input data used for the Open Grid Emissions dataset is released by the U.S. Energy Information Administration in the Autumn following the end of each year (2022 data was published in September 2023). Each release will include the most recent year of available data as well as updates of all previous available years based on any updates to the OGE methodology. All previous versions of the data will be archived on Zenodo.

Updated datasets will also be published whenever a new version of the open-grid-emissions repository is released.

### Running the pipeline with early release data
The OGE pipeline can be used to generate data using Early Release EIA data as soon as it is integrated into the PUDL nightly builds. In order to do that, `constants.latest_validated_year` must be changed to match `constants.current_early_release_year` before running the pipeline.

In addition, you will need to download and use the pudl nightly build data until the data becomes available through a stable release. To do so, you need to set your `PUDL_BUILD` environment variable to "nightly". You can do this through the command line using `set PUDL_BUILD=nightly` (for Windows), or by adding the following to the `__init__.py` file in `src/oge`:
```python
import os

os.environ["PUDL_BUILD"] = "nightly"
```

## Contribute
There are many ways that you can contribute!
 - Tell us how you are using the dataset or python tools
 - Request new features or data outputs by submitting a feature request or emailing us at <>
 - Tell us how we can make the datasets even easier to use
 - Ask a question about the data or methods in our [discussion forum](https://github.com/singularity-energy/open-grid-emissions/discussions)
 - [Submit an issue](https://github.com/singularity-energy/open-grid-emissions/issues) if you've identified a way the methods or assumptions could be improved
 - Contribute your subject matter expertise to the discussion about [open issues and questions](https://github.com/singularity-energy/open-grid-emissions/issues?q=is%3Aissue+is%3Aopen+label%3Aquestion)
 - Submit a pull request to help us fix open issues

## Repository Structure
### Modules
- `anomaly_screening`: classes use to flag timeseries for anomalies as proposed in Tyler H. Ruggles et al. Developing reliable hourly electricity demand data through screening and imputation (2020)
- `column_checks`: functions that check that all data outputs have the correct column names
- `constants`: specifies conversion factors and constants used across all modules
- `data_pipeline`: main script for running the data pipeline from start to finish
- `download_data`: functions that download data from the internet
- `data_cleaning`: functions that clean loaded data
- `eia930`: functions for cleaning and formatting EIA-930 data
- `emissions`: functions used for imputing emissions data
- `filepaths`: used to identify where repository files are located on the user's computer
- `gross_to_net_generation`: functions for identifying subplants and gross to net generation conversion factors
- `helpers`: functions that are used across modules
- `impute_hourly_profiles`: functions related to assigning an hourly profile to monthly data
- `load_data`: functions for loading data from downloaded files
- `output_data`: functions for writing intermediate and final data to csvs
- `subplant_identification`: functions for identifying subplant IDs
- `validation`: functions for testing and validating data outputs
- `visualization`: functions for visualizing data in notebooks

### Notebooks
Notebooks are organized into five directories based on their purpose
- `explore_data`: notebooks used for exploring data outputs and results
- `explore_methods`: notebooks that can be used to explore specific methods step-by-step
- `manual_data`: notebooks that are used to create/update certain files in `data/manual`
- `validation`: notebooks related to validating results
- `visualization`: notebooks used to visualize data
- `work_in_progress`: temporary notebooks being used for development purposes on specific branches

### Data Structure
All manual reference tables are stored in `src/oge/reference_tables`.

All files downloaded/created as part of the pipeline are stored in your HOME directory (e.g. users/user.name/):
- `$HOME/open_grid_emissions_data/downloads` contains all files that are downloaded by functions in `load_data`
- `$HOME/open_grid_emissions_data/outputs` contains intermediate outputs from the data pipeline... any files created by our code that are not final results
- `$HOME/open_grid_emissions_data/results` contains all final output files that will be published

## Importing OGE as a Package in your Project
OGE is not yet available on PyPi but can be installed from GitHub. For example, this can be done by adding `oge = {git="https://github.com/singularity-energy/open-grid-emissions.git"}` to your Pipfile if you are using `pipenv` for your project.

Note that you don't need to run the pipeline to generate the output data as these are available on Amazon Simple Storage Service (S3). Simply, set the `OGE_DATA_STORE` environment variable to `s3` in the **\_\_init\_\_.py** file of your project to fetch OGE data from Amazon S3.
To summarize, your **\_\_init\_\_.py** file would then look like this:
```python
import os

os.environ["OGE_DATA_STORE"] = "s3"
```

## Development Setup
If you would like to run the code on your own computer and/or contribute updates to the code, the following steps can help get you started.

### Setup with conda
This installation is recommended if you are unfamiliar with git and Python.
#### Install conda and python
We suggest using miniconda or Anaconda to manage the packages needed to run the Open Grid Emissions code. Anaconda and Miniconda install a similar environment, but Anaconda installs more packages by default and Miniconda installs them as needed. These can be downloaded from [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution)

#### Install and setup git software manager
In order to download the repository, you will need to use git. You can either install Git Bash from https://git-scm.com/downloads, or you can install it using conda. To do so, after installing Anaconda or Miniconda, open an Anaconda Command Prompt (Windows) or Terminal.app (Mac) and type the following command:
```bash
conda install git
```
Then you will need set up git following these instructions: https://docs.github.com/en/get-started/quickstart/set-up-git

#### Download the codebase to a local repository
Using Anaconda command prompt or Git Bash, use the `cd` and `mkdir` commands to create and/or enter the directory where you would like to download the code (e.g. "Users/myusername/GitHub"). Then run:
```bash
git clone https://github.com/singularity-energy/open-grid-emissions.git
```

#### Setup the conda environment
Open anaconda prompt, use `cd` to navigate to the directory where your local files are stored (e.g. "GitHub/open-grid-emissions"), and then run:
```bash
conda update conda
conda env create -f environment.yml
```
Installation requires that the conda channel-priority be set to "flexible". This is the default behavior, so if you've never manually changed this, you shouldn't have to worry about this. However, if you receive an error message like "Found conflicts!" when trying to install the environment, try setting your channel priority to flexible by running the following command:`conda config --set channel_priority flexible` and then re-running the above commands.

The final step is to install the `oge` package itself in the conda environment. To do so, run:
```bash
conda activate open_grid_emissions
pip install build
python -m build
pip install --editable .
```

The open_grid_emissions conda environment should now be set up and ready to run.

### Setup with pipenv
#### Install python and git
We recommend that you use Python 3.11. If you don't have Python installed, we recommend that you use [**pyenv**](https://github.com/pyenv/pyenv). It lets you easily switch between multiple versions of Python. You will also need to use git to clone the repository. It can be installed from https://git-scm.com/downloads,

#### Install pipenv
This can be done via:
```bash
pip install pipenv
```

#### Download the codebase
As mentioned previously, clone the repository with:
```bash
git clone https://github.com/singularity-energy/open-grid-emissions.git
```
and navigate to the root of the directory:
```bash
cd open-grid-emissions
```

#### Setup the environment
In the root of the directory, create and activate the environment with:
```bash
# set up virtual environment (use whichever version of python 3.11 you have installed)
pipenv --python 3.11.4

# if you have updated the pipfile and need to update pipfile.lock, run
pipenv install
# Otherwise, if you just want to install packages from the pipfile.lock, run
pipenv sync

# activate virtual environment
pipenv shell

# install an editable version of the oge package
pip install build
python -m build
pip install –-editable .
```

If you ever need to remove and reinstall the environment, run `pipenv --rm` from the root directory then follow the directions above.

### Running the complete data pipeline
If you would like to run the full data pipeline to generate all intermediate outputs and results files, navigate to `open-grid-emissions/src/oge`, and run the following (replacing 2022 with whichever year you want to run):
```bash
python data_pipeline.py --year 2022
```

### Keeping the code updated
From time to time, the code will be updated on GitHub. To ensure that you are keeping your local version of the code up to date, open git bash and follow these steps:
```bash
# change the directory to where ever your local git repository is saved
# after hitting enter, it should show the name of the git branch (e.g. "(main)")
cd GitHub/open-grid-emissions  

# save any changes that you might have made locally to your copy of the code
git add .

# fetch and merge the updated code from github
git pull origin main
```

### Install a code editor
If you want to edit the code and do not already have an integrated development environment (IDE) installed, one good option is Visual Studio Code (download: https://code.visualstudio.com/).

## Contribution Guidelines
If you plan on contributing edits to the codebase that will be merged into the main branch, please follow these best practices:

1. Please do not make edits directly to the main branch. Any new features or edits should be completed in a new branch. To do so, open git bash, navigate to your local repo (e.g. `cd GitHub/open-grid-emissions`), and create a new branch, giving it a descriptive name related to the edit you will be doing:

	`git checkout -b branch_name`

2. As you code, it is a good practice to 'save' your work frequently by opening git bash, navigating to your local repo (`cd GitHub/open-grid-emissions`), making sure that your current feature branch is active (you should see the feature name in parentheses next to the command line), and running

	`git add .`

3. You should commit your work to the branch whenever you have working code or whenever you stop working on it using:

	`git add .`  
	`git commit -m "short message about updates"`

4. Once you are done with your edits, save and commit your code using step #3 and then push your changes:

	`git push`

5. Now open the GitHub repo web page. You should see the branch you pushed up in a yellow bar at the top of the page with a button to "Compare & pull request".
	- Click "Compare & pull request". This will take you to the "Open a pull request" page.
	- From here, you should write a brief description of what you actually changed.
	- Click "Create pull request"
	- The changes will be reviewed and discussed. Once any edits have been made, the code will be merged into the main branch.

## Conventions and standards
- We generally follow the naming conventions used by the Public Utility Data Liberation Project: https://catalystcoop-pudl.readthedocs.io/en/latest/dev/naming_conventions.html
- Functions should include descriptive docstrings (using the Google style guide https://google.github.io/styleguide/pyguide.html#383-functions-and-methods), inline comments should be used to describe individual steps, and variable names should be made descriptive (e.g. `cems_plants_with_missing_co2_data` not `cems_missing` or `cpmco2`)
- All pandas merge operations should include the `validate` parameter to ensure that unintentional duplicate entries are not created (https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html)
- All pandas groupby operations should include the `dropna=False` parameter so that data with missing groupby keys are not unintentionally dropped from the data.
- All code should be formatted using `ruff`, running `ruff format` in the root of the repository will format all files according to the set of configurations enclosed in the `pyproject.toml` file.
- Clear all outputs from notebooks before committing your work.
- Any manual changes to reported categorical data, conversion factors, or manual data mappings should be loaded from a .csv file `src/oge/reference_tables` rather than stored in a dictionary or variable in the code.
