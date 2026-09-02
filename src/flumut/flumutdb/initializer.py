from peewee import SqliteDatabase

from flumut.core.globals import DATABASE_PROXY, DB_FILE
from flumut.flumutdb.models import DbVersion


def initialize(path: str | None = None, read_only: bool = True) -> None:
    if not path:
        path = DB_FILE

    if read_only:
        new_db = SqliteDatabase(f'file:{path}?mode=ro', uri=True, pragmas={'foreign_keys': 1})
    else:
        new_db = SqliteDatabase(path, pragmas={'foreign_keys': 1})

    DATABASE_PROXY.initialize(new_db)

    DbVersion.is_compatible()
