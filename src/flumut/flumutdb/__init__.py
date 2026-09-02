from flumut.core.globals import DATABASE_PROXY
from flumut.core.options import DatabaseOptions
from flumut.flumutdb import loader
from flumut.flumutdb.initializer import initialize
from flumut.flumutdb.models import (
    Annotation,
    BaseModel,
    DbVersion,
    Effect,
    Evidence,
    Host,
    Mapping,
    Marker,
    Mutation,
    MutationType,
    Paper,
    Protein,
    Reference,
    Segment,
    Subtype,
)


def open_database(options: DatabaseOptions) -> None:
    """Connect to the database described by ``options`` and read it into memory.

    Args:
        options: Which database to open, and whether it may be written to.
    """
    initialize(options.path, options.read_only)
    loader.load()


def close_database() -> None:
    """Close this thread's connection to the database, if it has one."""
    if DATABASE_PROXY.obj is None:
        return
    DATABASE_PROXY.close()
