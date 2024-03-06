import sqlite3


_connection: sqlite3.Connection
_cursor: sqlite3.Cursor

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

def execute_query(query: str):
    return _cursor.execute(query)
