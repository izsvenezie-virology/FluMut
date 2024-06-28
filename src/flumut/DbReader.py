import sqlite3
from importlib_resources import files
from typing import Any, Callable, Dict

_db_file: str = files('flumutdb').joinpath('flumut_db.sqlite')
_connection: sqlite3.Connection
_cursor: sqlite3.Cursor

def set_db_file(db_path: str) -> None:
    global _db_file
    _db_file = db_path

def get_db_version():
    open_connection()
    major, minor, date = _cursor.execute('SELECT * FROM db_version').fetchone()
    close_connection()
    return major, minor, date

def open_connection() -> None:
    global _connection
    global _cursor
    _connection = sqlite3.connect(_db_file)
    _cursor = _connection.cursor()

def close_connection() -> None:
    global _connection
    global _cursor
    _connection.close()
    _connection = None
    _cursor = None

def execute_query(query: str, row_factory: Callable = None):
    _cursor.row_factory = row_factory
    return _cursor.execute(query)

def to_dict(cursor: sqlite3.Cursor, row: sqlite3.Row) -> Dict[str, Any]:
    result = {}
    for idx, col in enumerate(cursor.description):
        result[col[0]] = row[idx]
    return result
