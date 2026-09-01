"""Exact-output regression test.

The tests in ``test_pipeline.py`` cannot notice a change that keeps the output
well-formed but alters the science, so this one compares the three TSV outputs
byte for byte against committed files.

The comparison is tied to the database the snapshot was recorded with, so a
FluMutDB update legitimately breaks it. That is deliberate: read the diff,
confirm the new markers or evidence explain it, then re-record with::

    uv run pytest tests/e2e --snapshot-update
"""

import difflib
import json

import pytest

from tests.e2e.conftest import SNAPSHOT_DIR, TSV_OUTPUTS, current_db_version

SNAPSHOT = SNAPSHOT_DIR / 'single_sample'
RECORDED_WITH = SNAPSHOT / 'recorded_with.json'


def diff_against_snapshot(name: str, actual: str) -> str:
    """Return a unified diff against the recorded output, or '' if they match."""
    path = SNAPSHOT / f'{name}.tsv'
    if not path.exists():
        return f'{name}.tsv has never been recorded.\n'
    expected = path.read_text(encoding='utf-8')
    if actual == expected:
        return ''
    diff = difflib.unified_diff(
        expected.splitlines(keepends=True),
        actual.splitlines(keepends=True),
        fromfile=f'recorded/{name}.tsv',
        tofile=f'current/{name}.tsv',
        n=1,
    )
    return ''.join(list(diff)[:40])


def test_single_sample_output_matches_snapshot(default_run, request) -> None:
    assert default_run.exception is None

    if request.config.getoption('--snapshot-update'):
        SNAPSHOT.mkdir(parents=True, exist_ok=True)
        for name in TSV_OUTPUTS:
            (SNAPSHOT / f'{name}.tsv').write_text(default_run.text(name), encoding='utf-8')
        recorded = {'flumut_db_version': current_db_version(), 'input': 'single_sample.fa'}
        RECORDED_WITH.write_text(json.dumps(recorded, indent=2) + '\n', encoding='utf-8')
        pytest.skip(f'Snapshot re-recorded with FluMutDB {current_db_version()}.')

    diffs = [diff for name in TSV_OUTPUTS if (diff := diff_against_snapshot(name, default_run.text(name)))]
    if not diffs:
        return

    recorded_version = json.loads(RECORDED_WITH.read_text(encoding='utf-8'))['flumut_db_version'] if RECORDED_WITH.exists() else 'unknown'
    pytest.fail(
        f'Output differs from the snapshot recorded with FluMutDB {recorded_version} '
        f'(currently loaded: {current_db_version()}).\n'
        'If the database changed, or the change is intended, re-record with:\n'
        '    uv run pytest tests/e2e --snapshot-update\n\n' + '\n'.join(diffs),
        pytrace=False,
    )
