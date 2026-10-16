import filecmp
import gzip
import os
import shutil

import pytest

from oge.physics_reconciliation import make_dataset

FIXTURE_FOLDER = os.path.join(os.path.dirname(__file__), "data")

# files written by `make_dataset` in every run
OUTPUT_FILES = [
    "eia930_basic.csv",
    "eia930_elec.csv",
    "eia930_objective.csv",
    "eia930_opt.csv",
    "eia930_raw.csv",
    "eia930_rolling.csv",
    "eia930_weights.csv",
]
# files written only when the data includes hours before 2018-07-01
PRE_2018_07_OUTPUT_FILES = [
    "eia930_objective_no_src.csv",
    "eia930_opt_no_src.csv",
]


@pytest.fixture(scope="module")
def gridemissions_make_dataset(tmp_path_factory):
    """`make_dataset` from the original `gridemissions` package.

    `gridemissions` is a dev-only dependency, kept so that `oge.physics_reconciliation`
    can be checked against the code it was copied from. Importing it writes config and
    logging files, so its config folder is pointed at a temporary folder first. The
    test is skipped if the package isn't installed.
    """
    config_folder = tmp_path_factory.mktemp("gridemissions_config")
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GRIDEMISSIONS_CONFIG_DIR_PATH", str(config_folder))
        mp.setenv("GRIDEMISSIONS_CONFIG_FILE_PATH", str(config_folder / "config.json"))
        workflows = pytest.importorskip("gridemissions.workflows")
    return workflows.make_dataset


def copy_raw_fixture(fixture_name: str, folder: str):
    """Decompresses a raw EIA-930 fixture to `eia930_raw.csv` in `folder`.

    Args:
        fixture_name (str): file name of the gzipped fixture in `test/data`.
        folder (str): folder to write `eia930_raw.csv` to. Created if it doesn't exist.
    """
    os.makedirs(folder, exist_ok=True)
    with (
        gzip.open(os.path.join(FIXTURE_FOLDER, fixture_name), "rb") as f_in,
        open(os.path.join(folder, "eia930_raw.csv"), "wb") as f_out,
    ):
        shutil.copyfileobj(f_in, f_out)


@pytest.mark.parametrize(
    "fixture_name, expected_files",
    [
        (
            "eia930_raw_2024-11-01_2024-11-14.csv.gz",
            OUTPUT_FILES,
        ),
        (
            "eia930_raw_2018-06-24_2018-07-07.csv.gz",
            OUTPUT_FILES + PRE_2018_07_OUTPUT_FILES,
        ),
    ],
)
def test_make_dataset_matches_gridemissions(
    fixture_name, expected_files, gridemissions_make_dataset, tmp_path
):
    """Every file written by `make_dataset` should be byte-for-byte identical to the
    output of the original `gridemissions` package.

    `oge.physics_reconciliation` is a copy of the `gridemissions` cleaning code, and
    later refactors must not change its results. Each fixture is two weeks of real
    EIA-930 data for every BA, in the format `clean_930` writes to `eia930_raw.csv`.
    Every BA is included because the optimization step couples BAs through
    interchange, so dropping a BA would change the results for the others.

    - The 2024 fixture covers the first weeks of storage reporting in several BAs
      (e.g. ISNE, ERCO, NEVP), when negative storage values started to appear.
    - The 2018 fixture crosses 2018-07-01, so it exercises both the optimization
      without generation by source (before) and with it (after).
    """
    gridemissions_folder = str(tmp_path / "gridemissions")
    oge_folder = str(tmp_path / "oge")
    copy_raw_fixture(fixture_name, gridemissions_folder)
    copy_raw_fixture(fixture_name, oge_folder)

    gridemissions_make_dataset(
        start="",
        end="",
        file_name="eia930",
        tmp_folder=gridemissions_folder,
        folder_hist=gridemissions_folder,
        scrape=False,
        add_ca_fuels=False,
        calc_consumed=False,
    )
    make_dataset(
        file_name="eia930",
        tmp_folder=oge_folder,
        folder_hist=oge_folder,
        add_ca_fuels=False,
    )

    # both versions write the same set of files, including the pre-2018-07 files
    # only when the fixture has hours before 2018-07-01
    assert sorted(os.listdir(gridemissions_folder)) == sorted(expected_files)
    assert sorted(os.listdir(oge_folder)) == sorted(expected_files)
    # each file is identical
    for file_name in expected_files:
        assert filecmp.cmp(
            os.path.join(gridemissions_folder, file_name),
            os.path.join(oge_folder, file_name),
            shallow=False,
        ), f"{file_name} differs from gridemissions output"
