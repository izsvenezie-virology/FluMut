"""Shared fixtures for the end-to-end tests.

These tests drive the real Click command against the bundled FluMutDB and the
FASTA files in ``examples/``, exercising the whole pipeline: FASTA parsing,
reference selection, alignment, translation, quality checks, marker scanning
and output writing.

A full run takes a couple of seconds, so the tests share a small number of
session-scoped runs rather than each paying for their own.
"""

import csv
import itertools
from dataclasses import dataclass
from pathlib import Path

import pytest
from click.testing import CliRunner

from flumut.cli import cli
from flumut.flumutdb import DbVersion

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_DIR = Path(__file__).parent / 'data'

SINGLE_SAMPLE_FASTA = REPO_ROOT / 'examples' / 'single_sample.fa'

#: CLI flag and file extension for every output the command can produce.
OUTPUTS = {
    'markers': ('-m', '.tsv'),
    'mutations': ('-M', '.tsv'),
    'literature': ('-l', '.tsv'),
    'excel': ('-x', '.xlsm'),
}
TSV_OUTPUTS = ('markers', 'mutations', 'literature')
ALL_OUTPUTS = TSV_OUTPUTS + ('excel',)

#: Column schemas produced by the ``get_*_data`` extractors in ``core/io/output.py``.
MARKERS_COLUMNS = ['Sample', 'Marker', 'Mutations in your sample', 'Effect', 'Subtype', 'Literature']
LITERATURE_COLUMNS = ['Short name', 'Title', 'Authors', 'Year', 'Journal', 'Link', 'DOI']


def pytest_addoption(parser) -> None:
    """Register ``--snapshot-update`` for re-recording the pinned outputs."""
    parser.addoption(
        '--snapshot-update',
        action='store_true',
        default=False,
        help='Re-record the end-to-end snapshot instead of comparing against it.',
    )


def current_db_version() -> str:
    """Return the version of the currently loaded FluMutDB, e.g. ``'7.0.0'``."""
    return str(DbVersion.get_or_none())


@dataclass
class RunResult:
    """Outcome of one CLI invocation, with helpers to read the files it wrote.

    Attributes:
        exit_code: Process exit status reported by the Click test runner.
        exception: The exception that escaped the command, or None.
        paths: Map from output name (``markers``, ``excel``, ...) to its path.
    """

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

    def text(self, name: str) -> str:
        """Return the raw text of one TSV output."""
        return self.paths[name].read_text(encoding='utf-8')


def _make_runner(base_dir: Path):
    """Build a callable that runs the CLI into numbered subdirectories of ``base_dir``."""
    counter = itertools.count()

    def run(fasta=SINGLE_SAMPLE_FASTA, outputs=ALL_OUTPUTS, extra_args=()) -> RunResult:
        run_dir = base_dir / f'run{next(counter)}'
        run_dir.mkdir()

        paths = {}
        args = []
        for name in outputs:
            flag, suffix = OUTPUTS[name]
            paths[name] = run_dir / f'{name}{suffix}'
            args += [flag, str(paths[name])]

        args += list(extra_args)
        fastas = [fasta] if isinstance(fasta, (str, Path)) else list(fasta)
        args += [str(path) for path in fastas]

        result = CliRunner().invoke(cli, args)
        return RunResult(result.exit_code, result.exception, paths)

    return run


@pytest.fixture
def flumut(tmp_path):
    """Run the FluMut CLI into a per-test directory. Use for one-off argument combinations."""
    return _make_runner(tmp_path)


@pytest.fixture(scope='session')
def default_run(tmp_path_factory) -> RunResult:
    """One run over ``single_sample.fa`` producing every output, shared by most tests."""
    return _make_runner(tmp_path_factory.mktemp('flumut-default'))()


@pytest.fixture(scope='session')
def relaxed_run(tmp_path_factory) -> RunResult:
    """The same input analysed with ``--relaxed``."""
    return _make_runner(tmp_path_factory.mktemp('flumut-relaxed'))(extra_args=('-r',))


@pytest.fixture(scope='session')
def tiny_fasta(tmp_path_factory) -> Path:
    """A single-record FASTA carved out of ``single_sample.fa``.

    Used by the reproducibility test, which has to pay for a fresh interpreter
    on every run and so wants the smallest input that still exercises the
    whole pipeline.
    """
    records = SINGLE_SAMPLE_FASTA.read_text(encoding='utf-8').split('\n>')
    path = tmp_path_factory.mktemp('flumut-tiny') / 'tiny.fa'
    path.write_text(records[0].rstrip() + '\n', encoding='utf-8')
    return path
