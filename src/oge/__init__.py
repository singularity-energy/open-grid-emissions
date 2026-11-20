import sqlalchemy as sa
import warnings

from oge.filepaths import pudl_folder, use_pudl_s3
from oge.logging_util import configure_root_logger

# Suppress specific PUDL UserWarning about schema field shadowing
warnings.filterwarnings(
    "ignore",
    message='Field name "schema" in "Resource" shadows an attribute in parent "PudlMeta"',
    category=UserWarning,
    module="pudl.metadata.classes",
)

# Set PUDL engine
PUDL_ENGINE = (
    None
    if use_pudl_s3()
    else sa.create_engine("sqlite:///" + pudl_folder("pudl.sqlite"))
)


# Set up the OGE logging configuration once.
configure_root_logger(logfile=None)
