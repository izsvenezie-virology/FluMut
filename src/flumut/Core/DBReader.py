import sqlite3
from typing import Callable

from importlib_resources import files


class DBReaderClass:
    def __init__(self) -> None:
        self._db_path: str = files('flumutdb').joinpath('flumut_db.sqlite')
        self._connection: sqlite3.Connection = None

    @property
    def db_path(self) -> str:
        return self._db_path

    @db_path.setter
    def db_path(self, db_path: str) -> None:
        self._db_path = db_path
        self.close_db()

    @property
    def connection_string(self) -> str:
        return f'file:{self.db_path}?mode=ro'

    @property
    def connection(self) -> sqlite3.Connection:
        if self._connection is None:
            self._connection = sqlite3.connect(self.connection_string, uri=True)
        return self._connection

    def execute_query(self, query: str, row_factory: Callable = None) -> sqlite3.Cursor:
        cursor = self.connection.cursor()
        cursor.row_factory = row_factory
        return cursor.execute(query)

    def __del__(self):
        self.close_db()

    def close_db(self) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None


DBReader = DBReaderClass()
