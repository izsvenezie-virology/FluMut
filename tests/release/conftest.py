"""A throwaway FluMut checkout to run the release flow against.

The release code drives real ``git`` and a real SQLite file, and the parts most
worth covering are exactly the ones a mock would paper over, so these tests
build a small but genuine repository instead.
"""

import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from flumut.core.globals import DATABASE_PROXY
from flumut.flumutdb.initializer import initialize
from flumut_db_editor.release import DATABASE_FILE, DUMP_SCRIPT, GLOBALS_FILE, VERSION_FILE

REPO_ROOT = Path(__file__).resolve().parents[2]

SCHEMA = """
CREATE TABLE "dbversion" ("id" INTEGER NOT NULL PRIMARY KEY, "notes" TEXT, "major" INTEGER, "minor" INTEGER, "patch" INTEGER);
CREATE TABLE "evidence" ("id" INTEGER NOT NULL PRIMARY KEY, "notes" TEXT, "paper_id" INTEGER);
CREATE TABLE "reference" ("id" INTEGER NOT NULL PRIMARY KEY, "notes" TEXT, "name" TEXT);
INSERT INTO "dbversion" ("id", "major", "minor", "patch") VALUES (1, 7, 0, 0);
INSERT INTO "evidence" ("id", "paper_id") VALUES (1, 1);
INSERT INTO "reference" ("id", "name") VALUES (1, 'A/reference/1');
"""

GLOBALS_SOURCE = """\
# Trimmed copy of src/flumut/core/globals.py for the release tests.
DB_MAJOR_VERSION = 7
"""


def run(*args: str, cwd: Path) -> str:
    return subprocess.run(args, cwd=cwd, check=True, capture_output=True, text=True).stdout


@pytest.fixture
def sandbox(tmp_path: Path):
    """A committed FluMut-shaped checkout on ``main``, whose database the editor has open.

    Its ``origin`` is a local bare repository, so a release can push without the network.
    """
    root = tmp_path / 'FluMut'
    for path in (DATABASE_FILE, VERSION_FILE, GLOBALS_FILE, DUMP_SCRIPT):
        (root / path).parent.mkdir(parents=True, exist_ok=True)

    (root / '.gitignore').write_text('__pycache__/\n', encoding='utf-8')
    (root / VERSION_FILE).write_text("__version__ = '1.0.0'\n__author__ = 'Test'\n", encoding='utf-8')
    (root / GLOBALS_FILE).write_text(GLOBALS_SOURCE, encoding='utf-8')
    shutil.copy(REPO_ROOT / DUMP_SCRIPT, root / DUMP_SCRIPT)

    with sqlite3.connect(root / DATABASE_FILE) as connection:
        connection.executescript(SCHEMA)
    run(sys.executable, DUMP_SCRIPT.as_posix(), cwd=root)

    run('git', 'init', '--bare', str(tmp_path / 'origin.git'), cwd=tmp_path)
    run('git', 'init', '--initial-branch', 'main', cwd=root)
    run('git', 'config', 'user.email', 'test@example.com', cwd=root)
    run('git', 'config', 'user.name', 'Test', cwd=root)
    run('git', 'remote', 'add', 'origin', str(tmp_path / 'origin.git'), cwd=root)
    run('git', 'add', '--all', cwd=root)
    run('git', 'commit', '--message', 'initial', cwd=root)

    previous = DATABASE_PROXY.obj
    initialize(str(root / DATABASE_FILE), read_only=False)
    yield root.resolve()
    DATABASE_PROXY.close()
    DATABASE_PROXY.initialize(previous)


@pytest.fixture
def edit(sandbox: Path):
    """Run SQL against the sandbox database, standing in for the editor."""

    def apply_sql(*statements: str) -> None:
        with sqlite3.connect(sandbox / DATABASE_FILE) as connection:
            for statement in statements:
                connection.execute(statement)

    return apply_sql
