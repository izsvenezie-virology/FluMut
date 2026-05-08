import sqlite3
from pathlib import Path
from typing import Optional, Tuple

from peewee import SqliteDatabase

import flumut.flumutdb.core as core
from flumut.flumutdb.core import BUNDLED_DB, DATABASE_PROXY, REQUIRED_MAJOR_VERSION, USER_DB
from flumut.flumutdb.models import DbVersion, Marker, Segment


def _read_db_version(path: Path) -> Optional[Tuple[int, int, str]]:
    if not path.exists():
        return None
    try:
        conn = sqlite3.connect(f'file:{path}?mode=ro', uri=True)
        try:
            cursor = conn.execute('SELECT major, minor, date FROM dbversion LIMIT 1')
            return cursor.fetchone()
        finally:
            conn.close()
    except Exception:
        return None


def _pick_best_db() -> Path:
    user_ver = _read_db_version(USER_DB)
    if user_ver is None or user_ver[0] != REQUIRED_MAJOR_VERSION:
        return BUNDLED_DB

    bundled_ver = _read_db_version(BUNDLED_DB)
    if bundled_ver is None:
        return USER_DB

    user_key = (user_ver[1], user_ver[2])
    bundled_key = (bundled_ver[1], bundled_ver[2])
    return USER_DB if user_key >= bundled_key else BUNDLED_DB


def initialize(db_path: str | None = None, read_only: bool = True) -> None:
    """Connect to the FluMutDB database.

    Args:
        db_path: Path to a SQLite database file. If not provided, automatically
                 selects the most up-to-date compatible database between the
                 bundled DB and any user-local DB in the platform data directory.
        read_only: Whether to open the database in read-only mode.

    Raises:
        FileNotFoundError: If an explicit db_path is given but does not exist.
        MissingVersionError: If the selected DB has no version record.
        IncompatibleVersionError: If the selected DB's major version is incompatible.
    """
    if db_path is not None:
        path = Path(db_path)
        if not path.exists():
            raise FileNotFoundError(f'Database not found: {path}')
    else:
        path = _pick_best_db()

    old_db = DATABASE_PROXY.obj
    if old_db is not None and not old_db.is_closed():
        old_db.close()

    if read_only:
        new_db = SqliteDatabase(f'file:{path}?mode=ro', uri=True, pragmas={'foreign_keys': 1})
    else:
        new_db = SqliteDatabase(path, pragmas={'foreign_keys': 1})

    DATABASE_PROXY.initialize(new_db)
    core.CURRENT_DB_PATH = path

    DbVersion.is_compatible()
    Segment.clear_cache()
    Marker.clear_cache()
