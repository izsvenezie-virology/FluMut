from peewee import SqliteDatabase

from flumut.core.globals import DB_FILE
from flumut.flumutdb import DATABASE_PROXY
from flumut.flumutdb.models import DbVersion, Marker, Segment


def initialize(path: str | None = None, read_only: bool = True) -> None:
    if not path:
        path = DB_FILE
    if read_only:
        path = f'file:{path}?mode=ro'
    new_db = SqliteDatabase(path, pragmas={'foreign_keys': 1})

    DATABASE_PROXY.initialize(new_db)

    DbVersion.is_compatible()
    Segment.clear_cache()
    Marker.clear_cache()
