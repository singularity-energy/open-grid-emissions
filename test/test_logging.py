import logging
import os

import pandas as pd

from oge.logging_util import get_logger, configure_root_logger

pudl_logger = logging.getLogger(name="catalystcoop.pudl")
LOGFILE = os.path.join(os.path.dirname(__file__), "test_logfile.txt")

configure_root_logger(logfile=LOGFILE, level=logging.INFO)
# If you call this again, nothing bad should happen. Logging statements should
# still only show up once.
configure_root_logger(logfile=LOGFILE, level=logging.INFO)
logger = get_logger("test")


def main():
    """These statements should each be printed once in a nice format."""
    logger.info("This is the OGE logger")
    pudl_logger.info("This is the PUDL logger")

    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    logger.info("\n" + df.to_string())


if __name__ == "__main__":
    main()


def test_configure_root_logger_twice_keeps_one_file_handler(tmp_path):
    """Configuring the same log file twice should leave a single handler writing to it.

    Otherwise every line would be written to the file twice, and the first handler's
    file would stay open. Handlers that write elsewhere should be left alone.
    """
    root_logger = logging.getLogger()
    other_handler = logging.NullHandler()
    root_logger.addHandler(other_handler)
    logfile = str(tmp_path / "test.log")
    try:
        configure_root_logger(logfile=logfile, level=logging.INFO)
        configure_root_logger(logfile=logfile, level=logging.INFO)
        file_handlers = [
            h
            for h in root_logger.handlers
            if isinstance(h, logging.FileHandler)
            and h.baseFilename == os.path.abspath(logfile)
        ]
        # one handler for this file, and the unrelated handler is still there
        assert len(file_handlers) == 1
        assert other_handler in root_logger.handlers

        get_logger("test").info("written once")
        file_handlers[0].flush()
        with open(logfile) as f:
            assert f.read().count("written once") == 1
    finally:
        for h in list(root_logger.handlers):
            if h is other_handler or (
                isinstance(h, logging.FileHandler)
                and h.baseFilename == os.path.abspath(logfile)
            ):
                root_logger.removeHandler(h)
                h.close()
