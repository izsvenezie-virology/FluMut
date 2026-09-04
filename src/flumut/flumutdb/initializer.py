import errno
import os

from peewee import SqliteDatabase

from flumut.core.globals import DATABASE_PROXY, DB_FILE


def initialize(path: str | None = None, read_only: bool = True) -> None:
    if not path:
        path = DB_FILE

    if read_only:
        if not os.path.isfile(path):
            raise FileNotFoundError(errno.ENOENT, 'Database file not found', path)
        new_db = SqliteDatabase(f'file:{path}?mode=ro', uri=True, pragmas={'foreign_keys': 1})
    else:
        new_db = SqliteDatabase(path, pragmas={'foreign_keys': 1})

    DATABASE_PROXY.initialize(new_db)
