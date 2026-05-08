import os
import sys
from pathlib import Path

from peewee import DatabaseProxy
from platformdirs import user_data_dir

REQUIRED_MAJOR_VERSION = 7


DATABASE_PROXY = DatabaseProxy()
ENV_DB: Path = Path(sys.prefix) / "share" / "flumutdb" / "flumut_db.sqlite"
BUNDLED_DB: Path = Path(__file__).parent / "data" / "flumut_db.sqlite"
USER_DB: Path = Path(user_data_dir("flumutdb")) / "flumut_db.sqlite"
CURRENT_DB_PATH: Path | None = None


def _resolve_user_db() -> Path:
    if getattr(sys, "frozen", False):
        return Path(user_data_dir("flumutdb")) / "flumut_db.sqlite"
    env_db = 
    check = env_db.parent
    while not check.exists():
        check = check.parent
    if os.access(check, os.W_OK):
        return env_db
    return Path(user_data_dir("flumutdb")) / "flumut_db.sqlite"
