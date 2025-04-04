import logging
import sqlite3
import sys
from typing import Any, Callable, Dict, Tuple

from importlib_resources import files

REQUIRED_DB_VERSION = 7
'''The major version of FluMutDB required to use this version of FluMut.'''


class DBConnection:
    _default_db: str = files('flumutdb').joinpath('flumut_db.sqlite')
    _db_file: str = _default_db
    '''Default database is the one stored in flumutdb package.'''
    _connection: sqlite3.Connection = None

    def __new__(cls):
        '''Singleton pattern.'''
        if not hasattr(cls, 'instance'):
            cls.instance = super(DBConnection, cls).__new__(cls)
        return cls.instance

    @property
    def db_file(self) -> str:
        '''The database file.'''
        return self._db_file

    @db_file.setter
    def db_file(self, db_path: str) -> None:
        self.close_db()
        self._db_file = db_path

    @property
    def connection_string(self) -> str:
        '''Connection string to open the database file in read only mode'''
        return f'file:{self.db_file}?mode=ro'

    @property
    def connection(self) -> sqlite3.Connection:
        '''The connection to the database.'''
        if self._connection is None:
            self._connection = sqlite3.connect(self.connection_string, uri=True)
            self.check_db_version()
        return self._connection

    @property
    def version(self) -> Tuple[str, str, str]:
        '''
        Major, minor versions and the date of the last release.

        :return `str` major: The major version of the database. Is bumped when the structure of the database changes.
        :return `str` minor: Minor version of the database. Is bumped when data change.
        :return `str` date: Date of the latest release.
        '''
        major, minor, date = self.execute_query('SELECT * FROM db_version').fetchone()
        return int(major), int(minor), date

    @property
    def version_string(self) -> str:
        '''Database version as a string.'''
        major, minor, date = self.version
        return f'v.{major}.{minor}, released on {date}'

    @property
    def is_pip_package(self) -> bool:
        '''Checks if the database is the flumutdb package one.'''
        return self._default_db == self.db_file

    def check_db_version(self):
        '''
        Checks if the FluMutDB version is compatible with this version of FluMut.
        Exits if versions are not compatible.
        '''
        db_major_version, _, _ = self.version
        if db_major_version < REQUIRED_DB_VERSION:
            logging.critical(f'FluMutDB {self.version_string} is too old.')
            logging.critical(f'This version of FluMut works only with FluMutDB v.{REQUIRED_DB_VERSION}.x.')
            logging.critical(f'Please, update FluMutDB with `flumut --update`.')
            sys.exit(1)
        if db_major_version > REQUIRED_DB_VERSION:
            logging.critical(f'FluMutDB {self.version_string} is too recent.')
            logging.critical(f'This version of FluMut works only with FluMutDB v.{REQUIRED_DB_VERSION}.x.')
            logging.critical(f'Please, search for FluMut updates.')
            sys.exit(1)

    def execute_query(self, query: str, row_factory: Callable = None) -> sqlite3.Cursor:
        '''
        Execute a query in the database.

        :param `str` query: The SQL query.
        :param `Callable` row_factory: Changes the format of the result.
        :return `sqlite3.Cursor`: The formatted result of the query.
        '''
        cursor = self.connection.cursor()
        cursor.row_factory = row_factory
        return cursor.execute(query)

    def __del__(self):
        self.close_db()

    def close_db(self) -> None:
        '''Closes the connection with the database.'''
        if self._connection is not None:
            self._connection.close()
            self._connection = None


def execute_query(query: str, row_factory: Callable = None) -> sqlite3.Cursor:
    return DBConnection().execute_query(query, row_factory)


def to_dict(cursor: sqlite3.Cursor, row: sqlite3.Row) -> Dict[str, Any]:
    '''Row factory to have each entry formatted as a dictionary with database column name as key and the value as value.'''
    result = {}
    for idx, col in enumerate(cursor.description):
        result[col[0]] = row[idx]
    return result
