"""End-to-end tests for user errors.

A mistake in the input, the options or the database must end in a one-line
message and a non-zero exit status, never in a raw traceback.
"""

import shutil
import sqlite3
from contextlib import closing

import pytest

from flumut.core.exceptions import FluMutError
from flumut.core.globals import DATABASE_PROXY, DB_FILE, DB_MAJOR_VERSION
from flumut.flumutdb.exceptions import FluMutDBError


@pytest.fixture(autouse=True)
def restore_database():
    """Put back the database the session was using, as these tests point FluMut at broken ones."""
    previous = DATABASE_PROXY.obj
    yield
    DATABASE_PROXY.initialize(previous)


def _invalid_character(tmp_path) -> dict:
    fasta = tmp_path / 'invalid.fa'
    fasta.write_text('>sample_PB2\nACGTZZACGT\n', encoding='utf-8')
    return {'fasta': fasta}


def _invalid_name_regex(tmp_path) -> dict:
    return {'extra_args': ('-n', '(')}


def _missing_database(tmp_path) -> dict:
    return {'extra_args': ('-D', str(tmp_path / 'missing.sqlite'))}


def _not_a_sqlite_file(tmp_path) -> dict:
    database = tmp_path / 'text.sqlite'
    database.write_text('this is not a database, and is long enough not to look like an empty one', encoding='utf-8')
    return {'extra_args': ('-D', str(database))}


def _sqlite_without_version_table(tmp_path) -> dict:
    database = tmp_path / 'other.sqlite'
    with closing(sqlite3.connect(database)) as connection:
        connection.execute('CREATE TABLE unrelated (id INTEGER)')
        connection.commit()
    return {'extra_args': ('-D', str(database))}


def _incompatible_major_version(tmp_path) -> dict:
    database = tmp_path / 'incompatible.sqlite'
    shutil.copy(DB_FILE, database)
    with closing(sqlite3.connect(database)) as connection:
        connection.execute('UPDATE dbversion SET major = ?', (DB_MAJOR_VERSION + 1,))
        connection.commit()
    return {'extra_args': ('-D', str(database))}


@pytest.mark.parametrize(
    'make_arguments',
    [
        _invalid_character,
        _invalid_name_regex,
        _missing_database,
        _not_a_sqlite_file,
        _sqlite_without_version_table,
        _incompatible_major_version,
    ],
    ids=lambda make_arguments: make_arguments.__name__.lstrip('_'),
)
def test_user_error_is_reported_without_a_traceback(flumut, tmp_path, make_arguments) -> None:
    """The CLI turns the error into ``SystemExit(1)``; anything else escaping would print a traceback."""
    result = flumut(outputs=('markers',), **make_arguments(tmp_path))

    assert result.exit_code == 1
    assert isinstance(result.exception, SystemExit)
    assert isinstance(result.exception.__context__, (FluMutError, FluMutDBError))
