import sqlite3
import urllib.request
from typing import Any, Callable, Dict


_connection: sqlite3.Connection
_cursor: sqlite3.Cursor

def update_db(db_path: str) -> None:
    url = 'https://github.com/izsvenezie-virology/FluMutDB/releases/latest/download/FluMutDB.sqlite3'
    _, _ = urllib.request.urlretrieve(url, db_path)

def open_connection(db_path: str) -> None:
    global _connection
    global _cursor
    _connection = sqlite3.connect(db_path)
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
