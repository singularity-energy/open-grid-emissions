"""Convenience functions for paths."""

import os
from importlib.metadata import version


def get_data_store():
    """Set data location"""
    store = os.getenv("OGE_DATA_STORE")
    if store is None:
        return os.path.join(os.path.expanduser("~"), "open_grid_emissions_data")
    elif store == "1" or store.lower() == "local":
        return os.path.join(os.path.expanduser("~"), "open_grid_emissions_data")
    elif store == "2" or store.lower() == "s3":
        # get the data version based on the current OGE code version.
        # data versions are only published as X.X.0 versions
        oge_data_version = (
            f"{version('oge').split('.')[0]}.{version('oge').split('.')[1]}.0"
        )
        return f"s3://open-grid-emissions/{oge_data_version}/open_grid_emissions_data"


def top_folder(rel=""):
    """Returns a path relative to the top-level repo folder. This will work regardless
    of where the function is imported or called from.
    """
    return os.path.join(
        os.path.abspath(os.path.join(os.path.realpath(__file__), "../")), rel
    )


def reference_table_folder(rel=""):
    return os.path.join(top_folder("reference_tables"), rel).replace("\\", "/")


def data_folder(rel=""):
    """Returns a path relative to the `data` folder."""
    return os.path.join(get_data_store(), rel).replace("\\", "/")


def downloads_folder(rel=""):
    return os.path.join(data_folder("downloads"), rel).replace("\\", "/")


def outputs_folder(rel=""):
    return os.path.join(data_folder("outputs"), rel).replace("\\", "/")


def results_folder(rel=""):
    return os.path.join(data_folder("results"), rel).replace("\\", "/")


def containing_folder(filepath: str) -> str:
    """Returns the folder containing `filepath`."""
    return os.path.dirname(os.path.realpath(filepath))


def make_containing_folder(filepath: str):
    """Make sure the the folder where `filepath` goes exists."""
    os.makedirs(containing_folder(filepath), exist_ok=True)
