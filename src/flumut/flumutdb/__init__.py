from peewee import DatabaseError

from flumut.core.globals import DATABASE_PROXY, DB_FILE
from flumut.core.options import DatabaseOptions
from flumut.flumutdb import loader
from flumut.flumutdb.exceptions import InvalidDatabaseError
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

#: The models and helpers the rest of FluMut, and the database editor, import
#: from here rather than reaching into ``flumut.flumutdb.models``.
__all__ = [
    'Annotation',
    'BaseModel',
    'DbVersion',
    'Effect',
    'Evidence',
    'Host',
    'Mapping',
    'Marker',
    'Mutation',
    'MutationType',
    'Paper',
    'Protein',
    'Reference',
    'Segment',
    'Subtype',
    'close_database',
    'initialize',
    'loader',
    'open_database',
    'read_version',
]


def read_version(options: DatabaseOptions) -> DbVersion | None:
    """Connect to the database described by ``options`` and read its version.

    Args:
        options: Which database to open, and whether it may be written to.

    Returns:
        The version of the database, or None if it does not declare one.

    Raises:
        DatabaseNotFoundError: If a read-only database file does not exist.
        InvalidDatabaseError: If the file is not a SQLite database, or has no version table.
    """
    initialize(options.path, options.read_only)
    try:
        return DbVersion.get_or_none()
    except DatabaseError as e:
        raise InvalidDatabaseError(options.path or DB_FILE, e) from e


def open_database(options: DatabaseOptions) -> None:
    """Connect to the database described by ``options`` and read it into memory.

    Args:
        options: Which database to open, and whether it may be written to.
    """
    read_version(options)
    DbVersion.is_compatible()
    loader.load()


def close_database() -> None:
    """Close this thread's connection to the database, if it has one."""
    if DATABASE_PROXY.obj is None:
        return
    DATABASE_PROXY.close()
