# --------------------------------------------------------------------------------------
# Options:
#  - Use `pytest -rP` to show print statements from PASSED tests after they finish
#  - Use `pytest -s` to direct print statements to the console
# 
# Run a specific test case with:
# pytest name_of_this_file.py -rP -k 'name_of_test_function'
#
# NOTE: The required input data must be downloaded first. See data_pipeline.py.

import sys
sys.path.append('../')

import src.data_cleaning as data_cleaning


def test_export_eia923_2019_and_2020():
    """These are the nominal years where data cleaning is well-tested."""
    data_cleaning.clean_eia923(2019, False, add_subplant_id=False)
    data_cleaning.clean_eia923(2020, False, add_subplant_id=False)


def test_export_eia923_2018():
    data_cleaning.clean_eia923(2018, False, add_subplant_id=False)
