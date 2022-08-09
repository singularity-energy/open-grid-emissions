# Open Grid Emissions Initiative
The Open Grid Emissions Initiative seeks to fill a critical need for high-quality, publicly-accessible, hourly grid emissions data that can be used for GHG accounting, policymaking, academic research, and energy attribute certificate markets. The initiative includes this repository of open-source grid emissions data processing tools that use peer-reviewed, well-documented, and validated methodologies to create the accompanying public dataset of hourly, monthly, and annual U.S. electric grid generation, GHG, and air pollution data.

## Contribute
There are many ways that you can contribute!
 - Tell us how you are using the dataset or python tools
 - Request new features or data outputs by submitting a feature request or emailing us at ...
 - Tell us how we can make the datasets even easier to use
 - Ask a question about the data or methods in our discussion forum
 - Submit an issue if you've identified a way the methods or assumptions could be improved
 - Contribute your subject matter expertise to the discussion about open issues and questions
 - Submit a pull request to help us fix open issues

# Repository Structure
### Modules
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

### Notebooks
Notebooks are organized into five directories based on their purpose
- `explore_data`: notebooks used for exploring data outputs and results
- `explore_methods`: notebooks that can be used to explore specific methods step-by-step
- `manual_data`: notebooks that are used to create/update certain files in `data/manual`
- `validation`: notebooks related to validating results
- `work_in_progress`: temporary notebooks being used for development purposes on specific branches

### Data Structure
- `data/downloads` contains all files that are downloaded by functions in `load_data`
- `data/manual` contains all manually-created files, including the egrid static tables
- `data/outputs` contains intermediate outputs from the data pipeline... any files created by our code that are not final results
- `data/results` contains all final output files that will be published

# Development Setup and contribution guidelines

If you would like to run the code on your own computer and/or contribute updates to the code, the following steps can help get you started.

## Users unfamiliar with git / python

### Install conda and python

We suggest using miniconda or Anaconda to manage the packages needed to run the Open Grid Emissions code. Anaconda and Miniconda install a similar environment, but Anaconda installs more packages by default and Miniconda installs them as needed. These can be downloaded from [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution)

### Install a code editor

If you want to edit the code and do not already have an integrated development environment (IDE) installed, one good option is Visual Studio Code (download: https://code.visualstudio.com/). 

### Install and setup git software manager

In order to download the repository, you will need to use git. You can either install Git Bash from https://git-scm.com/downloads, or you can install it using conda. To do so, fter installing Anaconda or Miniconda, open an Anaconda Command Prompt (Windows) or Terminal.app (Mac) and type the following command:

```
conda install git
```

Then you will need set up git following these instructions: https://docs.github.com/en/get-started/quickstart/set-up-git

## Once you have git and conda installed

### Download the codebase to a local repository

Using Anaconda command prompt or Git Bash, use the `cd` and `mkdir` commands to create and/or enter the directory where you would like to download the code (e.g. "Users/myusername/GitHub"). Then run:

```
git clone https://github.com/singularity-energy/open-grid-emissions.git
```

### Setup the conda environment

Open anaconda prompt, use `cd` to navigate to the directory where your local files are stored (e.g. "GitHub/open-grid-emissions"), and then run:

```
conda update conda
conda env create -f environment.yml
```

Installation requires that the conda channel-priority be set to "flexible". This is the default behavior, 
so if you've never manually changed this, you shouldn't have to worry about this. However, 
if you receive an error message like "Found conflicts!" when trying to install the environment,
try setting your channel priority to flexible by running the following command:
`conda config --set channel_priority flexible` and then re-running the above commands.

## Running the complete data pipeline

If you would like to run the full data pipeline to generate all intermediate outputs and results files, open anaconda prompt, navigate to the directory where the repository is saved, and run the following (replacing 2020 with whichever year you want to run):

```
conda activate open_grid_emissions
python src/data_pipeline.py --year 2020
```

## Keeping the code updated

From time to time, the code will be updated on GitHub. To ensure that you are keeping your local version of the code up to date, open git bash and follow these steps:
```
# change the directory to where ever your local git repository is saved
# after hitting enter, it should show the name of the git branch (e.g. "(master)")
cd GitHub/open-grid-emissions  

# save any changes that you might have made locally to your copy of the code
git add .

# fetch and merge the updated code from github
git pull origin master
```

## Contributing edits to the code

If you plan on contributing edits to the codebase that will be merged into the master branch, please follow these best practices:

1. Please do not make edits directly to the master branch. Any new features or edits should be completed in a new branch. To do so, open git bash, navigate to your local repo (e.g. `cd GitHub/open-grid-emissions`), and create a new branch, giving it a descriptive name related to the edit you will be doing:

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
	- The changes will be reviewed and discussed. Once any edits have been made, the code will be merged into the master branch.
