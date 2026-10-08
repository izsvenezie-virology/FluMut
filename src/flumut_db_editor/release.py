"""Propose a FluMutDB release from the database edited in this checkout.

A release puts the new versions, the SQL dump and the re-recorded snapshot on a
``db/vX.Y.Z`` branch, pushes it, and opens GitHub's pull request page. It never
merges or tags: merging the pull request is the review gate. The bump rules are
the ones written down in ``VERSIONING.md``.
"""

import re
import runpy
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlencode

from flumut.core.globals import DATABASE_PROXY
from flumut.flumutdb import DbVersion

BASE_BRANCH = 'main'

#: Paths a database release may change, relative to the checkout root. The
#: globals hold DB_MAJOR_VERSION, which a major release changes by hand.
DATABASE_FILE = Path('src/flumut/data/flumut_db.sqlite')
DUMP_FILE = Path('db/flumut_db.sql')
VERSION_FILE = Path('src/flumut/__init__.py')
GLOBALS_FILE = Path('src/flumut/core/globals.py')
SNAPSHOT_DIR = Path('tests/e2e/data/single_sample')
RELEASE_PATHS = (DATABASE_FILE, DUMP_FILE, VERSION_FILE, GLOBALS_FILE, SNAPSHOT_DIR)

DUMP_SCRIPT = Path('scripts/dump_db.py')

PREPROCESS_TABLES = {'segment', 'protein', 'reference', 'annotation'}

DUMP_ROW = re.compile(r'^INSERT INTO "(?P<table>[^"]+)" \([^)]*\) VALUES \((?P<key>[^,)]+)')
VERSION_ASSIGNMENT = re.compile(r"""^(__version__\s*=\s*['"])([^'"]+)""", re.MULTILINE)
GITHUB_REMOTE = re.compile(r'github\.com[:/](?P<repo>[^/]+/[^/]+?)(?:\.git)?/?$')

REVIEW_CHECKLIST = """
## Review

- [ ] The change to `db/flumut_db.sql` is the intended one.
- [ ] The change to the recorded outputs in `tests/e2e/data/` reflects the new evidence.
- [ ] The bump level matches `VERSIONING.md`.
"""


class ReleaseError(RuntimeError):
    """Raised when a release cannot be proposed or carried out."""


@dataclass(frozen=True)
class Proposal:
    """The release the edited database calls for."""

    root: Path
    level: str
    changes: dict[str, int]
    db_version: str
    flumut_version: str

    @property
    def new_db_version(self) -> str:
        return bump(self.db_version, self.level)

    @property
    def new_flumut_version(self) -> str:
        # A database update ships as a FluMut patch: the data got better, the code did not change.
        return bump(self.flumut_version, 'patch')

    @property
    def branch(self) -> str:
        return f'db/v{self.new_db_version}'

    @property
    def title(self) -> str:
        return f'FluMutDB v.{self.new_db_version}'

    @property
    def default_summary(self) -> str:
        changed = ', '.join(f'{table} ({rows})' for table, rows in self.changes.items()) or 'schema only'
        return f'{self.level.capitalize()} database update: {changed}.'


def bump(version: str, level: str) -> str:
    """``version`` moved up by ``level``, with the lower parts reset."""
    major, minor, patch = map(int, version.split('.'))
    if level == 'major':
        return f'{major + 1}.0.0'
    if level == 'minor':
        return f'{major}.{minor + 1}.0'
    return f'{major}.{minor}.{patch + 1}'


def bump_db_version(level: str) -> str:
    """Move the version recorded in the open database up by ``level``, and return it."""
    version = DbVersion.get()
    version.major, version.minor, version.patch = map(int, bump(str(version), level).split('.'))
    version.save()
    return str(version)


def changed_tables(before: str, after: str) -> tuple[dict[str, int], bool]:
    """Compare two SQL dumps: the rows changed in each table, and whether the schema changed.

    Every row is dumped on its own line, so the lines only one dump has are
    exactly the rows added, removed or edited.
    """
    keys: dict[str, set[str]] = {}
    schema_changed = False
    for line in set(before.splitlines()) ^ set(after.splitlines()):
        if match := DUMP_ROW.match(line):
            keys.setdefault(match['table'], set()).add(match['key'])
        elif line.strip() and not line.startswith('--'):
            schema_changed = True
    return {table: len(rows) for table, rows in sorted(keys.items())}, schema_changed


def release_level(changes: dict[str, int], schema_changed: bool) -> str | None:
    """The bump a change calls for, or None when only the version row moved."""
    if schema_changed:
        return 'major'
    tables = set(changes) - {'dbversion'}
    if not tables:
        return None
    return 'minor' if tables & PREPROCESS_TABLES else 'patch'


def git(root: Path, *args: str) -> str:
    result = subprocess.run(['git', *args], cwd=root, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise ReleaseError(f'git {args[0]} failed: {result.stderr.strip() or result.stdout.strip()}')
    return result.stdout.rstrip()


def checkout_root() -> Path:
    """The FluMut checkout the open database belongs to."""
    database = Path(DATABASE_PROXY.obj.database).resolve()
    try:
        root = Path(git(database.parent, 'rev-parse', '--show-toplevel')).resolve()
    except ReleaseError as error:
        raise ReleaseError(f'The open database, {database}, is not inside a git checkout of FluMut.') from error
    if root / DATABASE_FILE != database:
        raise ReleaseError(f'The open database is {database}, not the one a release ships ({root / DATABASE_FILE}).')
    return root


def dump(root: Path) -> str:
    """The database as SQL, rendered by the same script that wrote the committed dump."""
    script = runpy.run_path(str(root / DUMP_SCRIPT))
    with script['connect'](root / DATABASE_FILE) as connection:
        return script['render'](connection)


def rewrite(path: Path, pattern: re.Pattern, value: str) -> None:
    """Replace the value ``pattern`` captures in its second group with ``value``."""
    updated, count = pattern.subn(rf'\g<1>{value}', path.read_text(encoding='utf-8'), count=1)
    if count != 1:
        raise ReleaseError(f'Could not find {pattern.pattern!r} in {path}.')
    path.write_text(updated, encoding='utf-8', newline='\n')


def prepare() -> Proposal:
    """Work out the release the open database calls for, without changing anything.

    Raises:
        ReleaseError: If the checkout is not in a state a release can be cut from.
    """
    root = checkout_root()
    problems = []

    if (branch := git(root, 'branch', '--show-current')) != BASE_BRANCH:
        problems.append(
            f'The checkout is on "{branch}", not "{BASE_BRANCH}". '
            f'Switch to {BASE_BRANCH} so the pull request contains only the database change.'
        )

    changed = git(root, 'diff', '--name-only', 'HEAD').splitlines() + git(root, 'ls-files', '--others', '--exclude-standard').splitlines()
    allowed = tuple(path.as_posix() for path in RELEASE_PATHS)
    if unexpected := [path for path in changed if not path.startswith(allowed)]:
        more = f' (and {len(unexpected) - 5} more)' if len(unexpected) > 5 else ''
        problems.append(
            f'The working tree has changes a database release does not cover: {", ".join(unexpected[:5])}{more}. Commit or stash them.'
        )

    if problems:
        raise ReleaseError('\n'.join(problems))

    changes, schema_changed = changed_tables(git(root, 'show', f'HEAD:{DUMP_FILE.as_posix()}'), dump(root))
    level = release_level(changes, schema_changed)
    if level is None:
        raise ReleaseError('The database is unchanged from the last commit, so there is nothing to release.')

    changes.pop('dbversion', None)
    flumut_version = VERSION_ASSIGNMENT.search((root / VERSION_FILE).read_text(encoding='utf-8'))[2]
    return Proposal(root, level, changes, str(DbVersion.get()), flumut_version)


def release(proposal: Proposal, summary: str, snapshot: bool = True) -> None:
    """Commit the release on its own branch and push it.

    The branch is created first, so nothing is committed to the base branch. If a
    step fails, the checkout stays on the new branch with whatever was written.
    """
    root = proposal.root
    git(root, 'switch', '--create', proposal.branch)

    bump_db_version(proposal.level)
    (root / DUMP_FILE).write_text(dump(root), encoding='utf-8', newline='\n')
    rewrite(root / VERSION_FILE, VERSION_ASSIGNMENT, proposal.new_flumut_version)

    if snapshot:
        # A change that alters a marker call legitimately changes the recorded
        # outputs, and that diff is what the reviewer reads to confirm the science.
        command = [sys.executable, '-m', 'pytest', 'tests/e2e', '--snapshot-update', '-q']
        result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise ReleaseError(f'Could not re-record the snapshot:\n{result.stdout}\n{result.stderr}')

    git(root, 'add', '--', *(path.as_posix() for path in RELEASE_PATHS if (root / path).exists()))
    git(root, 'commit', '--message', f'feat (db): {proposal.title}\n\n{description(proposal, summary)}')
    git(root, 'push', '--set-upstream', 'origin', proposal.branch)


def description(proposal: Proposal, summary: str) -> str:
    """What changed, for the commit message and the pull request."""
    return '\n'.join(
        [
            summary or proposal.default_summary,
            '',
            f'Bump: {proposal.level} (FluMutDB {proposal.db_version} -> {proposal.new_db_version})',
            f'FluMut: {proposal.flumut_version} -> {proposal.new_flumut_version}',
        ]
    )


def pull_request_url(proposal: Proposal, summary: str) -> str:
    """GitHub's page for opening the release pull request, with title and description filled in."""
    remote = git(proposal.root, 'remote', 'get-url', 'origin')
    if (match := GITHUB_REMOTE.search(remote)) is None:
        raise ReleaseError(f'The "origin" remote is not a GitHub repository: {remote}')

    query = urlencode({'expand': 1, 'title': proposal.title, 'body': description(proposal, summary) + '\n' + REVIEW_CHECKLIST})
    return f'https://github.com/{match["repo"]}/compare/{BASE_BRANCH}...{proposal.branch}?{query}'
