import re
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from flumut.core.analysis.models import Analysis
from flumut.core.analysis.preprocess import load_nucleotide_fasta, parse_header


# ---------------------------------------------------------------------------
# parse_header tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'header, pattern_str, expected_sample, expected_segment',
    [
        # No match: whole header used as sample, segment is None
        ('sample1_PB2', r'NOMATCH', 'sample1_PB2', None),
        # Named groups
        ('sample2_PB2', r'(?P<sample>.+)_(?P<segment>.+)', 'sample2', 'PB2'),
        # Named groups inverted order
        ('PB2_sample3', r'(?P<segment>.+)_(?P<sample>.+)', 'sample3', 'PB2'),
        # Two positional groups (no named groups)
        ('sample5_PB2', r'(.+)_(.+)', 'sample5', 'PB2'),
        # One positional group: group(2) raises IndexError, segment stays None
        ('sample6_PB2', r'(.+)_', 'sample6', None),
        # No capture groups: group(1) raises IndexError, sample stays None, falls back to header
        ('sample7_PB2', r'.+_.+', 'sample7_PB2', None),
    ],
)
def test_parse_header(
    header: str, pattern_str: str, expected_sample: str, expected_segment: str | None
) -> None:
    sample, segment = parse_header(header, re.compile(pattern_str))
    assert sample == expected_sample
    assert segment == expected_segment


def test_parse_header_empty_sample_falls_back_to_header() -> None:
    # An empty-string capture for sample is falsy, so the full header is returned instead.
    header = '_PB2'
    pattern = re.compile(r'(?P<sample>.*)_(?P<segment>.+)')
    sample, segment = parse_header(header, pattern)
    assert sample == header
    assert segment == 'PB2'


# ---------------------------------------------------------------------------
# load_nucleotide_fasta helpers and fixtures
# ---------------------------------------------------------------------------


def _make_sequence(header: str) -> MagicMock:
    seq = MagicMock()
    seq.header = header
    return seq


def _make_nt_alignment(num_proteins: int = 1) -> MagicMock:
    nt = MagicMock()
    nt.reference.segment.proteins = [MagicMock() for _ in range(num_proteins)]
    return nt


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


_PATTERN = re.compile(r'(?P<sample>.+)')


# ---------------------------------------------------------------------------
# load_nucleotide_fasta tests
# ---------------------------------------------------------------------------


def test_load_single_sequence_creates_sample(mock_deps) -> None:
    mock_deps.read_fasta.return_value = [_make_sequence('sample1')]
    mock_deps.get_alignment.return_value = _make_nt_alignment(1)
    analysis = Analysis()
    load_nucleotide_fasta(analysis, MagicMock(), _PATTERN)
    assert 'sample1' in analysis.samples
    assert len(analysis.samples['sample1'].alignments) == 1


def test_load_empty_fasta_leaves_analysis_unchanged(mock_deps) -> None:
    mock_deps.read_fasta.return_value = []
    analysis = Analysis()
    load_nucleotide_fasta(analysis, MagicMock(), _PATTERN)
    assert analysis.samples == {}


def test_load_multiple_sequences_different_samples(mock_deps) -> None:
    mock_deps.read_fasta.return_value = [_make_sequence('sample1'), _make_sequence('sample2')]
    mock_deps.get_alignment.return_value = _make_nt_alignment(1)
    analysis = Analysis()
    load_nucleotide_fasta(analysis, MagicMock(), _PATTERN)
    assert set(analysis.samples.keys()) == {'sample1', 'sample2'}


def test_load_same_sample_id_reuses_existing_sample(mock_deps) -> None:
    mock_deps.read_fasta.return_value = [_make_sequence('sample1'), _make_sequence('sample1')]
    mock_deps.get_alignment.return_value = _make_nt_alignment(1)
    analysis = Analysis()
    load_nucleotide_fasta(analysis, MagicMock(), _PATTERN)
    assert len(analysis.samples) == 1
    assert len(analysis.samples['sample1'].alignments) == 2


def test_load_multiple_proteins_appended_per_sample(mock_deps) -> None:
    mock_deps.read_fasta.return_value = [_make_sequence('sample1')]
    mock_deps.get_alignment.return_value = _make_nt_alignment(3)
    analysis = Analysis()
    load_nucleotide_fasta(analysis, MagicMock(), _PATTERN)
    assert len(analysis.samples['sample1'].alignments) == 3


def test_load_translate_called_for_each_protein(mock_deps) -> None:
    nt_alignment = _make_nt_alignment(2)
    mock_deps.read_fasta.return_value = [_make_sequence('sample1')]
    mock_deps.get_alignment.return_value = nt_alignment
    load_nucleotide_fasta(Analysis(), MagicMock(), _PATTERN)
    proteins = nt_alignment.reference.segment.proteins
    assert mock_deps.get_cds.call_count == 2
    mock_deps.get_cds.assert_any_call(nt_alignment, proteins[0])
    mock_deps.get_cds.assert_any_call(nt_alignment, proteins[1])
    assert mock_deps.translate.call_count == 2
    mock_deps.translate.assert_any_call(mock_deps.get_cds.return_value)
    mock_deps.translate.assert_any_call(mock_deps.get_cds.return_value)


def test_load_translate_result_stored_in_alignments(mock_deps) -> None:
    aa_alignment = MagicMock()
    mock_deps.translate.return_value = aa_alignment
    mock_deps.read_fasta.return_value = [_make_sequence('sample1')]
    mock_deps.get_alignment.return_value = _make_nt_alignment(1)
    analysis = Analysis()
    load_nucleotide_fasta(analysis, MagicMock(), _PATTERN)
    assert analysis.samples['sample1'].alignments == [aa_alignment]
