"""Tests for the release flow: the rules in VERSIONING.md, and what a release changes."""

from pathlib import Path

import pytest

from flumut_db_editor.release import (
    DUMP_FILE,
    GLOBALS_FILE,
    VERSION_FILE,
    ReleaseError,
    bump,
    changed_tables,
    prepare,
    pull_request_url,
    release,
    release_level,
)

from .conftest import run

ADD_EVIDENCE = 'INSERT INTO evidence (id, paper_id) VALUES (99, 1)'
EDIT_REFERENCE = "UPDATE reference SET name = 'A/reference/2' WHERE id = 1"


@pytest.mark.parametrize(
    ('level', 'expected'),
    [('major', '8.0.0'), ('minor', '7.4.0'), ('patch', '7.3.2')],
)
def test_bump_resets_the_lower_parts(level: str, expected: str) -> None:
    assert bump('7.3.1', level) == expected


@pytest.mark.parametrize(
    ('tables', 'schema_changed', 'expected'),
    [
        ({'evidence', 'paper'}, False, 'patch'),
        ({'evidence', 'reference'}, False, 'minor'),
        ({'evidence'}, True, 'major'),
        ({'dbversion'}, False, None),
        (set(), False, None),
    ],
    ids=['analysis', 'preprocessing', 'schema', 'version-only', 'nothing'],
)
def test_the_highest_level_a_change_touches_wins(tables: set[str], schema_changed: bool, expected: str | None) -> None:
    assert release_level(dict.fromkeys(tables, 1), schema_changed) == expected


def test_changed_tables_counts_added_removed_and_edited_rows_once() -> None:
    before = '-- FluMutDB v.7.0.0\nCREATE TABLE "paper" (...);\nINSERT INTO "paper" ("id", "title") VALUES (1, \'A\');\nINSERT INTO "paper" ("id", "title") VALUES (2, \'B\');'
    after = '-- FluMutDB v.7.0.1\nCREATE TABLE "paper" (...);\nINSERT INTO "paper" ("id", "title") VALUES (1, \'A2\');\nINSERT INTO "paper" ("id", "title") VALUES (3, \'C\');'

    assert changed_tables(before, after) == ({'paper': 3}, False)


def test_changed_tables_notices_a_schema_change() -> None:
    assert changed_tables('CREATE TABLE "paper" ("id");', 'CREATE TABLE "paper" ("id", "doi");') == ({}, True)


def test_prepare_proposes_a_patch_for_new_evidence(sandbox: Path, edit) -> None:
    edit(ADD_EVIDENCE)

    proposal = prepare()

    assert proposal.level == 'patch'
    assert proposal.changes == {'evidence': 1}
    assert proposal.new_db_version == '7.0.1'
    assert proposal.new_flumut_version == '1.0.1'
    assert proposal.branch == 'db/v7.0.1'


def test_prepare_proposes_a_minor_for_a_reference(sandbox: Path, edit) -> None:
    edit(EDIT_REFERENCE)

    assert prepare().new_db_version == '7.1.0'


def test_prepare_refuses_when_nothing_changed(sandbox: Path) -> None:
    with pytest.raises(ReleaseError, match='unchanged'):
        prepare()


def test_prepare_refuses_another_branch(sandbox: Path, edit) -> None:
    """Branching off a feature branch would put its commits in the database pull request."""
    edit(ADD_EVIDENCE)
    run('git', 'switch', '--create', 'some-work', cwd=sandbox)

    with pytest.raises(ReleaseError, match='not "main"'):
        prepare()


def test_prepare_refuses_unrelated_changes(sandbox: Path, edit) -> None:
    edit(ADD_EVIDENCE)
    (sandbox / 'notes.txt').write_text('scratch', encoding='utf-8')

    with pytest.raises(ReleaseError, match='notes.txt'):
        prepare()


def test_release_commits_every_version_and_pushes_the_branch(sandbox: Path, edit) -> None:
    edit(ADD_EVIDENCE)

    release(prepare(), 'Add the 2024 H5N1 evidence', snapshot=False)

    assert run('git', 'branch', '--show-current', cwd=sandbox).strip() == 'db/v7.0.1'
    assert run('git', 'status', '--porcelain', cwd=sandbox) == '', 'everything the release touched should be committed'
    assert "__version__ = '1.0.1'" in (sandbox / VERSION_FILE).read_text(encoding='utf-8')
    assert run('git', 'ls-remote', '--heads', 'origin', 'db/v7.0.1', cwd=sandbox)

    commit = run('git', 'show', 'HEAD', '--', DUMP_FILE.as_posix(), cwd=sandbox)
    assert 'feat (db): FluMutDB v.7.0.1' in commit
    assert 'Add the 2024 H5N1 evidence' in commit
    assert '+INSERT INTO "evidence" ("id", "paper_id", "notes") VALUES (99, 1, NULL);' in commit
    assert '+INSERT INTO "dbversion" ("id", "major", "minor", "patch", "notes") VALUES (1, 7, 0, 1, NULL);' in commit


def test_a_major_release_commits_the_hand_edited_db_major_version(sandbox: Path, edit) -> None:
    """DB_MAJOR_VERSION is changed by hand, so the release must accept and ship that edit."""
    edit('ALTER TABLE evidence ADD COLUMN weight REAL')
    (sandbox / GLOBALS_FILE).write_text('DB_MAJOR_VERSION = 8\n', encoding='utf-8')

    proposal = prepare()
    release(proposal, '', snapshot=False)

    assert proposal.new_db_version == '8.0.0'
    assert 'DB_MAJOR_VERSION = 8' in run('git', 'show', 'HEAD', '--', GLOBALS_FILE.as_posix(), cwd=sandbox)


def test_the_pull_request_url_fills_in_the_github_compare_page(sandbox: Path, edit) -> None:
    edit(ADD_EVIDENCE)
    run('git', 'remote', 'set-url', 'origin', 'git@github.com:izsvenezie-virology/FluMut.git', cwd=sandbox)

    url = pull_request_url(prepare(), 'Add the 2024 H5N1 evidence')

    assert url.startswith('https://github.com/izsvenezie-virology/FluMut/compare/main...db/v7.0.1?expand=1&')
    assert 'title=FluMutDB+v.7.0.1' in url


def test_the_pull_request_url_needs_a_github_remote(sandbox: Path, edit) -> None:
    edit(ADD_EVIDENCE)

    with pytest.raises(ReleaseError, match='not a GitHub repository'):
        pull_request_url(prepare(), '')
