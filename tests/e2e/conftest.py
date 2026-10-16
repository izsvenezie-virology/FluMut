"""Shared fixtures for the end-to-end tests.

These tests drive the real Click command against the bundled FluMutDB and
``examples/single_sample.fa``. A full run takes a couple of seconds, so most
tests share one session-scoped run.
"""

import csv
from dataclasses import dataclass
from pathlib import Path

import pytest
from click.testing import CliRunner

from flumut.cli import cli

REPO_ROOT = Path(__file__).resolve().parents[2]
SINGLE_SAMPLE_FASTA = REPO_ROOT / 'examples' / 'single_sample.fa'

#: CLI flag and file name for every output the command can produce.
OUTPUTS = {
    'markers': ('-m', 'markers.tsv'),
    'mutations': ('-M', 'mutations.tsv'),
    'literature': ('-l', 'literature.tsv'),
    'excel': ('-x', 'excel.xlsm'),
}


@dataclass
class RunResult:
    exit_code: int
    exception: BaseException | None
    paths: dict[str, Path]

    def rows(self, name: str) -> list[dict[str, str]]:
        """Parse one TSV output into a list of row dicts."""
        with self.paths[name].open(encoding='utf-8') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def columns(self, name: str) -> list[str]:
        """Return the header row of one TSV output."""
        with self.paths[name].open(encoding='utf-8') as handle:
            return next(csv.reader(handle, delimiter='\t'))


def run_flumut(directory: Path, fasta: Path = SINGLE_SAMPLE_FASTA, outputs=tuple(OUTPUTS), extra_args=()) -> RunResult:
    """Run the CLI on ``fasta``, writing the requested outputs into ``directory``."""
    paths = {name: directory / OUTPUTS[name][1] for name in outputs}
    args = [arg for name, path in paths.items() for arg in (OUTPUTS[name][0], str(path))]
    result = CliRunner().invoke(cli, [*args, *extra_args, str(fasta)])
    return RunResult(result.exit_code, result.exception, paths)


@pytest.fixture
def flumut(tmp_path):
    """Run the FluMut CLI into a per-test directory. Use for one-off argument combinations."""
    return lambda **kwargs: run_flumut(tmp_path, **kwargs)


@pytest.fixture(scope='session')
def default_run(tmp_path_factory) -> RunResult:
    """One run over ``single_sample.fa`` producing every output, shared by most tests."""
    return run_flumut(tmp_path_factory.mktemp('default'))


@pytest.fixture(scope='session')
def relaxed_run(tmp_path_factory) -> RunResult:
    """The same input analysed with ``--relaxed``."""
    return run_flumut(tmp_path_factory.mktemp('relaxed'), extra_args=('-r',))
