# Running Tests

We use `pytest` for unit testing.

You can run a specific test file with: `pytest path_to_test_file.py`.

Note that input data for the data pipeline may need to be downloaded before you can run
certain units tests. See `data_pipeline.py`.

**Options:**
 - Normally, pytest will only show print statements if a test case FAILS
 - Use `pytest -rP` to show print statements from PASSED tests after they finish
 - Use `pytest -s` to direct print statements to the console as they happen
 - Run a specific test function with `pytest -k name_of_function`
