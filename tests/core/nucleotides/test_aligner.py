from unittest.mock import MagicMock, patch

import pytest

from flumut.core.nucleotides.aligner import get_best_alignment, select_candidate_references
from flumut.core.globals import GAP_SYMBOL
from flumut.core.nucleotides.models import Nucleotide


def _make_reference(name: str, segment_name: str) -> MagicMock:
    ref = MagicMock()
    ref.name = name
    ref.segment.name = segment_name
    return ref


@pytest.fixture
def mock_all():
    with patch('flumut.core.nucleotides.aligner.Reference') as mock_ref_class:
        yield mock_ref_class.all


# ---------------------------------------------------------------------------
# select_candidate_references tests
# ---------------------------------------------------------------------------


def test_none_hint_returns_all_references(mock_all) -> None:
    refs = [_make_reference('PB2_ref', 'PB2'), _make_reference('PB1_ref', 'PB1')]
    mock_all.return_value = refs
    assert select_candidate_references(None) is refs


def test_empty_string_hint_returns_all_references(mock_all) -> None:
    refs = [_make_reference('PB2_ref', 'PB2'), _make_reference('PB1_ref', 'PB1')]
    mock_all.return_value = refs
    assert select_candidate_references('') is refs


def test_hint_matches_nothing_returns_all_references(mock_all) -> None:
    refs = [_make_reference('PB2_ref', 'PB2'), _make_reference('PB1_ref', 'PB1')]
    mock_all.return_value = refs
    assert select_candidate_references('UNKNOWN') is refs


def test_hint_matches_reference_name(mock_all) -> None:
    pb2 = _make_reference('PB2_ref', 'PB2')
    pb1 = _make_reference('PB1_ref', 'PB1')
    mock_all.return_value = [pb2, pb1]
    assert select_candidate_references('PB2_ref') == [pb2]


def test_hint_matches_segment_name(mock_all) -> None:
    pb2 = _make_reference('PB2_ref', 'PB2')
    pb1 = _make_reference('PB1_ref', 'PB1')
    mock_all.return_value = [pb2, pb1]
    assert select_candidate_references('PB2') == [pb2]


def test_hint_matches_multiple_references_in_same_segment(mock_all) -> None:
    pb2_a = _make_reference('PB2_ref_A', 'PB2')
    pb2_b = _make_reference('PB2_ref_B', 'PB2')
    pb1 = _make_reference('PB1_ref', 'PB1')
    mock_all.return_value = [pb2_a, pb2_b, pb1]
    assert select_candidate_references('PB2') == [pb2_a, pb2_b]


# ---------------------------------------------------------------------------
# get_best_alignment helpers and fixture
# ---------------------------------------------------------------------------


def _make_query(sequence: str) -> MagicMock:
    q = MagicMock()
    q.sequence = sequence
    return q


def _make_seq_reference(sequence: str) -> MagicMock:
    ref = MagicMock()
    ref.sequence = sequence
    return ref


def _make_aln_result(score: float, ref_chars: list[str], query_chars: list[str]) -> MagicMock:
    aln = MagicMock()
    aln.score = score
    aln.__getitem__ = MagicMock(side_effect=lambda i: ref_chars if i == 0 else query_chars)
    return aln


@pytest.fixture
def mock_aligner():
    with patch('flumut.core.nucleotides.aligner._aligner') as mock:
        yield mock


# ---------------------------------------------------------------------------
# get_best_alignment tests
# ---------------------------------------------------------------------------


def test_empty_candidates_raises_value_error(mock_aligner) -> None:
    with pytest.raises(ValueError, match='No reference found'):
        get_best_alignment(_make_query('ACGT'), candidates=[])


def test_single_candidate_returns_nucleotide(mock_aligner) -> None:
    query = _make_query('ACGT')
    reference = _make_seq_reference('ACGT')
    aln = _make_aln_result(score=10.0, ref_chars=['A', 'C'], query_chars=['A', 'C'])
    mock_aligner.align.return_value = [aln]
    result = get_best_alignment(query, [reference])
    assert isinstance(result, Nucleotide)
    assert result.query is query
    assert result.reference is reference


def test_highest_scoring_reference_is_selected(mock_aligner) -> None:
    query = _make_query('ACGT')
    ref_low = _make_seq_reference('AAAA')
    ref_high = _make_seq_reference('ACGT')
    aln_low = _make_aln_result(score=1.0, ref_chars=['A'], query_chars=['A'])
    aln_high = _make_aln_result(score=100.0, ref_chars=['A', 'C'], query_chars=['A', 'C'])
    mock_aligner.align.side_effect = [[aln_low], [aln_high]]
    result = get_best_alignment(query, [ref_low, ref_high])
    assert result.reference is ref_high


def test_gap_symbols_stripped_from_query_before_alignment(mock_aligner) -> None:
    query = _make_query('AC' + GAP_SYMBOL + 'GT')
    reference = _make_seq_reference('ACGT')
    aln = _make_aln_result(score=10.0, ref_chars=['A'], query_chars=['A'])
    mock_aligner.align.return_value = [aln]
    get_best_alignment(query, [reference])
    mock_aligner.align.assert_called_once_with('ACGT', 'ACGT')


def test_alignment_characters_stored_in_result(mock_aligner) -> None:
    ref_chars = ['A', 'C', GAP_SYMBOL, 'T']
    query_chars = ['A', 'C', 'G', 'T']
    aln = _make_aln_result(score=10.0, ref_chars=ref_chars, query_chars=query_chars)
    mock_aligner.align.return_value = [aln]
    result = get_best_alignment(_make_query('ACGT'), [_make_seq_reference('ACGT')])
    assert result.alignment.reference == ref_chars
    assert result.alignment.query == query_chars
