import sys
import logging

sys.path.append('../src')
sys.path.append('..')

import src.eia930 as eia930
from src.filepaths import top_folder

from src.logging_util import get_logger, configure_logger

pudl_logger = logging.getLogger(name="catalystcoop.pudl")

configure_logger(logfile=top_folder('test/test_logfile.txt'), level=logging.INFO)
# If you call this again, nothing bad should happen. Logging statements should
# still only show up once.
configure_logger(logfile=top_folder('test/test_logfile.txt'), level=logging.INFO)
logger = get_logger('test')


def main():
  """These statements should each be printed once in a nice format."""
  logger.info('This is the OGE logger')
  pudl_logger.info('This is the PUDL logger')


if __name__ == '__main__':
  main()
