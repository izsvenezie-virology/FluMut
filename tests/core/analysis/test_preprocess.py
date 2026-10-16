import re
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from flumut.core.analysis.models import Analysis
from flumut.core.analysis.preprocess import load_nucleotide_fasta, parse_header

_PATTERN = re.compile(r'(?P<sample>.+)')


# ---------------------------------------------------------------------------
# parse_header
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'header, pattern, expected_sample, expected_segment',
    [
        ('sample1_PB2', r'NOMATCH', 'sample1_PB2', None),
        ('sample2_PB2', r'(?P<sample>.+)_(?P<segment>.+)', 'sample2', 'PB2'),
        ('PB2_sample3', r'(?P<segment>.+)_(?P<sample>.+)', 'sample3', 'PB2'),
        ('sample5_PB2', r'(.+)_(.+)', 'sample5', 'PB2'),
        # One positional group: group(2) raises IndexError, so segment stays None.
        ('sample6_PB2', r'(.+)_', 'sample6', None),
        # No capture groups: group(1) raises IndexError, so the header is used as sample.
        ('sample7_PB2', r'.+_.+', 'sample7_PB2', None),
        # An empty capture is falsy, so the header is used as sample.
        ('_PB2', r'(?P<sample>.*)_(?P<segment>.+)', '_PB2', 'PB2'),
    ],
    ids=[
        'no_match',
        'named_groups',
        'named_groups_inverted',
        'positional_groups',
        'one_positional_group',
        'no_capture_groups',
        'empty_sample_capture',
    ],
)
def test_parse_header(header: str, pattern: str, expected_sample: str, expected_segment: str | None) -> None:
    """Named groups win, positional groups are the fallback, and anything unusable yields the header."""
    assert parse_header(header, re.compile(pattern)) == (expected_sample, expected_segment)


# ---------------------------------------------------------------------------
# load_nucleotide_fasta
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_deps():
    with (
        patch('flumut.core.analysis.preprocess.read_fasta') as read_fasta,
        patch('flumut.core.analysis.preprocess.select_candidate_references') as select_candidates,
        patch('flumut.core.analysis.preprocess.get_best_alignment') as get_alignment,
        patch('flumut.core.analysis.preprocess.get_cds') as get_cds,
        patch('flumut.core.analysis.preprocess.translate') as translate,
    ):
        yield SimpleNamespace(
            read_fasta=read_fasta,
            select_candidates=select_candidates,
            get_alignment=get_alignment,
            get_cds=get_cds,
            translate=translate,
        )


def _make_sequence(header: str) -> MagicMock:
    sequence = MagicMock()
    sequence.header = header
    return sequence


def _make_nt_alignment(num_proteins: int) -> MagicMock:
    nt = MagicMock()
    nt.reference.segment.proteins = [MagicMock() for _ in range(num_proteins)]
    return nt


@pytest.mark.parametrize(
    'headers, num_proteins, expected',
    [
        ([], 1, {}),
        (['sample1'], 1, {'sample1': 1}),
        (['sample1', 'sample2'], 1, {'sample1': 1, 'sample2': 1}),
        (['sample1', 'sample1'], 1, {'sample1': 2}),
        (['sample1'], 3, {'sample1': 3}),
    ],
    ids=['empty_fasta', 'one_sequence', 'two_samples', 'same_sample_twice', 'three_proteins'],
)
def test_load_nucleotide_fasta(mock_deps, headers: list[str], num_proteins: int, expected: dict[str, int]) -> None:
    """One alignment is stored per protein per sequence, grouped under the parsed sample ID."""
    mock_deps.read_fasta.return_value = [_make_sequence(header) for header in headers]
    mock_deps.get_alignment.return_value = _make_nt_alignment(num_proteins)
    analysis = Analysis()

    load_nucleotide_fasta(analysis, MagicMock(), _PATTERN)

    assert {sample_id: len(sample.alignments) for sample_id, sample in analysis.samples.items()} == expected
    for sample in analysis.samples.values():
        assert sample.alignments == [mock_deps.translate.return_value] * len(sample.alignments)


