import sqlite3
from typing import Callable, Self

from importlib_resources import files


class DBReader:
    def __init__(self) -> None:
        self._db_path: str = files('flumutdb').joinpath('flumut_db.sqlite')
        self._connection: sqlite3.Connection = None

    def __new__(cls) -> Self:
        if not hasattr(cls, 'instance'):
            cls.instance = super(DBReader, cls).__new__(cls)
        return cls.instance

    @property
    def db_path(self) -> str:
        return self._db_path

    @db_path.setter
    def db_path(self, db_path: str) -> None:
        self.close_db()
        self._db_path = db_path

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
