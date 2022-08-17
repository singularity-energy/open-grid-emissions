import os

# Convenience functions for paths.


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


class InputDataFilenames:
    # EIA-860 filenames.
    EIA_860_ENVIRO_ASSOC_FILE_FMT = '6_1_EnviroAssoc_Y{}.xlsx'  # .format(year)
    EIA_860_ENVIRO_EQUIP_FILE_FMT = '6_2_EnviroEquip_Y{}.xlsx'  # .format(year)
