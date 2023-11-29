# Set up the OGE logging configuration once.
import logging
from .logging_util import configure_root_logger
from .filepaths import outputs_folder

configure_root_logger(outputs_folder("logfile.txt"), logging.INFO)
