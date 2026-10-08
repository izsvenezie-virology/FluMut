"""Checks that the bundled database and the code that ships with it agree.

These are the gates the release cycle leans on. They run in the ordinary test
suite rather than only in a workflow, so a database change that forgets a step
fails on the contributor's machine, before it reaches a pull request.
"""

import importlib.util
from contextlib import closing
from pathlib import Path

import pytest

from flumut import __version__
from flumut.core.globals import DB_FILE, DB_MAJOR_VERSION
from flumut.core.updates import parse_version
from flumut.flumutdb import DbVersion

REPO_ROOT = Path(__file__).resolve().parents[1]
DUMP_SCRIPT = REPO_ROOT / 'scripts' / 'dump_db.py'
DUMP_FILE = REPO_ROOT / 'db' / 'flumut_db.sql'


@pytest.fixture(scope='module')
def db_version() -> DbVersion:
    """The version row of the bundled database."""
    from flumut.flumutdb import initialize

    initialize(None, read_only=True)
    version = DbVersion.get_or_none()
    assert version is not None, 'The bundled database has no dbversion row.'
    return version


def test_bundled_database_exists() -> None:
    """The package ships the database it analyses with."""
    assert Path(DB_FILE).is_file(), f'No bundled database at {DB_FILE}'


def test_database_major_matches_the_code(db_version: DbVersion) -> None:
    """``DB_MAJOR_VERSION`` declares the schema the bundled database actually has.

    These drift apart when a migration bumps one and not the other, and the
    result is a FluMut that refuses to open its own database.
    """
    assert db_version.major == DB_MAJOR_VERSION, (
        f'The bundled database is schema v{db_version.major} but DB_MAJOR_VERSION is {DB_MAJOR_VERSION}. '
        f'Bump DB_MAJOR_VERSION in src/flumut/core/globals.py, or the dbversion row.'
    )


def test_database_has_exactly_one_version_row() -> None:
    """A second version row would make ``DbVersion.get_or_none()`` a coin toss."""
    assert DbVersion.select().count() == 1


def test_flumut_version_is_readable() -> None:
    """``__version__`` parses, so the update check can compare against it."""
    assert parse_version(__version__) is not None, f'__version__ = {__version__!r} is not a MAJOR.MINOR.PATCH version'


def test_sql_dump_is_committed() -> None:
    """The dump is tracked in git, not just generated locally."""
    assert DUMP_FILE.is_file(), f'{DUMP_FILE} is missing. Run: python scripts/dump_db.py'


def test_sql_dump_matches_the_database() -> None:
    """The dump describes the bundled database as it is now.

    Pull requests review a database change through the dump, and the release
    attaches it, so a change that skips ``scripts/dump_db.py`` would show and
    ship a stale picture of the data.
    """
    spec = importlib.util.spec_from_file_location('dump_db', DUMP_SCRIPT)
    dump_db = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dump_db)

    with closing(dump_db.connect(Path(DB_FILE))) as connection:
        rendered = dump_db.render(connection)

    assert DUMP_FILE.read_text(encoding='utf-8') == rendered, f'{DUMP_FILE} is out of date. Run: python scripts/dump_db.py'
