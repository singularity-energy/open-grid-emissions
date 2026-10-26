"""Configure logging for the OGE codebase."""

import logging
import os

import coloredlogs

from oge.filepaths import make_containing_folder


def get_logger(name: str) -> logging.Logger:
    """Helper function to append `oge` to the logger name and return a logger.

    As a result, all returned loggers a children of the top-level `oge` logger.
    """
    return logging.getLogger(f"oge.{name}")


def configure_root_logger(logfile: str | None = None, level: str = "INFO"):
    """Configure the OGE logger to print to the console, and optionally to a file.

    This function is safe to call multiple times, since it will check if logging
    handlers have already been installed and skip them if so.

    Logging is printed with the same format as PUDL:
    ```
    2023-02-21 16:10:44 [INFO] oge.test:21 This is an example
    ```
    """
    root_logger = logging.getLogger()

    oge_logger = logging.getLogger("oge")
    log_format = "%(asctime)s [%(levelname)4s] %(name)s:%(lineno)s %(message)s"

    # Direct the output of the OGE logger to the terminal (and color it). Make
    # sure this hasn't been done already to avoid adding duplicate handlers.
    if len(oge_logger.handlers) == 0:
        coloredlogs.install(fmt=log_format, level=level, logger=oge_logger)
        oge_logger.addHandler(logging.NullHandler())

    # Send everything to the log file by adding a file handler to the root logger.
    if logfile is not None:
        make_containing_folder(logfile)
        # if a handler already writes to this file (e.g. when this function is called
        # twice), close and replace it so that each line is only written once. Other
        # handlers on the root logger are left alone.
        logfile_path = os.path.abspath(logfile)
        for handler in list(root_logger.handlers):
            if (
                isinstance(handler, logging.FileHandler)
                and handler.baseFilename == logfile_path
            ):
                root_logger.removeHandler(handler)
                handler.close()
        file_logger = logging.FileHandler(logfile, mode="w")
        file_logger.setFormatter(logging.Formatter(log_format))
        root_logger.addHandler(file_logger)
