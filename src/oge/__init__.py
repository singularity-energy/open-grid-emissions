import os
import sqlalchemy as sa

from oge.filepaths import pudl_folder
from oge.logging_util import configure_root_logger

# Set PUDL engine
PUDL_ENGINE = (
    sa.create_engine("sqlite:///" + pudl_folder("pudl.sqlite"))
    if os.getenv("PUDL_DATA_STORE") not in ["2", "s3"]
    else None
)


# Set up the OGE logging configuration once.
configure_root_logger(logfile=None)
