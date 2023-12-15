# Set up the OGE logging configuration once.
import logging
from oge.logging_util import configure_root_logger
from oge.filepaths import outputs_folder

configure_root_logger(outputs_folder("logfile.txt"), logging.INFO)
