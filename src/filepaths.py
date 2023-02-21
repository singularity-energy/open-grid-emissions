"""Convenience functions for paths."""
import os


def top_folder(rel=""):
    """
    Returns a path relative to the top-level repo folder.

    This will work regardless of where the function is imported or called from.
    """
    return os.path.join(
        os.path.abspath(os.path.join(os.path.realpath(__file__), "../../")), rel
    )


def data_folder(rel=""):
    """Returns a path relative to the `data` folder."""
    return os.path.join(top_folder("data"), rel)


def downloads_folder(rel=""):
    return os.path.join(data_folder("downloads"), rel)


def manual_folder(rel=""):
    return os.path.join(data_folder("manual"), rel)


def results_folder(rel=""):
    return os.path.join(data_folder("results"), rel)


def outputs_folder(rel=""):
    return os.path.join(data_folder("outputs"), rel)
