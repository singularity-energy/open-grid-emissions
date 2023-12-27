---
stoplight-id: contributing
---

## Contribute
There are many ways that you can contribute!
 - Tell us how you are using the dataset or python tools
 - Request new features or data outputs by submitting a feature request or emailing us at <>
 - Tell us how we can make the datasets even easier to use
 - Ask a question about the data or methods in our [discussion forum](https://github.com/singularity-energy/open-grid-emissions/discussions)
 - [Submit an issue](https://github.com/singularity-energy/open-grid-emissions/issues) if you've identified a way the methods or assumptions could be improved
 - Contribute your subject matter expertise to the discussion about [open issues and questions](https://github.com/singularity-energy/open-grid-emissions/issues?q=is%3Aissue+is%3Aopen+label%3Aquestion)
 - Submit a pull request to help us fix open issues

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
- We generally follow the [naming conventions](https://catalystcoop-pudl.readthedocs.io/en/latest/dev/naming_conventions.html) used by the Public Utility Data Liberation Project:
- Functions should include descriptive docstrings (using the [Google style guide](https://google.github.io/styleguide/pyguide.html#383-functions-and-methods)), inline comments should be used to describe individual steps, and variable names should be made descriptive (e.g. `cems_plants_with_missing_co2_data` not `cems_missing` or `cpmco2`)
- All pandas [merge operations](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html) should include the `validate` parameter to ensure that unintentional duplicate entries are not created
- All pandas groupby operations should include the `dropna=False` parameter so that data with missing groupby keys are not unintentionally dropped from the data.
- All code should be formatted using `ruff`
- Clear all outputs from notebooks before committing your work.
- Any manual changes to reported categorical data, conversion factors, or manual data mappings should be loaded from a .csv file `data/manual` rather than stored in a dictionary or variable in the code.