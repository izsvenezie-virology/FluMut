from unittest.mock import MagicMock, patch

import pytest

from flumut.core.globals import GAP_SYMBOL
from flumut.core.nucleotides.aligner import get_best_alignment, select_candidate_references
from flumut.core.nucleotides.models import Nucleotide

#: (reference name, segment name) for the references every selection test sees.
REFERENCES = [('PB2_ref_A', 'PB2'), ('PB2_ref_B', 'PB2'), ('PB1_ref', 'PB1')]
ALL_NAMES = [name for name, _ in REFERENCES]


def _make_reference(name: str, segment_name: str = 'PB2', sequence: str = 'ACGT') -> MagicMock:
    ref = MagicMock()
    ref.name = name
    ref.segment.name = segment_name
    ref.sequence = sequence
    return ref


def _make_query(sequence: str) -> MagicMock:
    query = MagicMock()
    query.sequence = sequence
    return query


def _make_aln_result(score: float, reference: list[str], query: list[str]) -> MagicMock:
    aln = MagicMock()
    aln.score = score
    aln.__getitem__ = MagicMock(side_effect=lambda i: reference if i == 0 else query)
    return aln


@pytest.fixture
def mock_get_references():
    with patch('flumut.core.nucleotides.aligner.loader.get') as mock_get:
        mock_get.return_value = [_make_reference(name, segment) for name, segment in REFERENCES]
        yield mock_get


@pytest.fixture
def mock_aligner():
    with patch('flumut.core.nucleotides.aligner._aligner') as mock:
        yield mock


# ---------------------------------------------------------------------------
# select_candidate_references
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'hint, expected',
    [
        (None, ALL_NAMES),
        ('', ALL_NAMES),
        ('UNKNOWN', ALL_NAMES),
        ('PB2_ref_A', ['PB2_ref_A']),
        ('PB2', ['PB2_ref_A', 'PB2_ref_B']),
        ('PB1', ['PB1_ref']),
    ],
    ids=['no_hint', 'empty_hint', 'unmatched_hint', 'reference_name', 'segment_with_two', 'segment_with_one'],
)
def test_select_candidate_references(mock_get_references, hint: str | None, expected: list[str]) -> None:
    """A hint narrows to a reference or a segment; anything unusable falls back to all references."""
    assert [ref.name for ref in select_candidate_references(hint)] == expected


# ---------------------------------------------------------------------------
# get_best_alignment
# ---------------------------------------------------------------------------


def test_get_best_alignment_no_candidates_raises(mock_aligner) -> None:
    with pytest.raises(ValueError, match='No reference found'):
        get_best_alignment(_make_query('ACGT'), candidates=[])


def test_get_best_alignment_returns_query_reference_and_alignment(mock_aligner) -> None:
    query = _make_query('ACGT')
    reference = _make_reference('PB2_ref_A')
    aligned_reference = ['A', 'C', GAP_SYMBOL, 'T']
    aligned_query = ['A', 'C', 'G', 'T']
    mock_aligner.align.return_value = [_make_aln_result(10.0, aligned_reference, aligned_query)]

    result = get_best_alignment(query, [reference])

    assert isinstance(result, Nucleotide)
    assert result.query is query
    assert result.reference is reference
    assert result.alignment.reference == aligned_reference
    assert result.alignment.query == aligned_query


def test_get_best_alignment_picks_the_highest_scoring_reference(mock_aligner) -> None:
    low = _make_reference('low', sequence='AAAA')
    high = _make_reference('high', sequence='ACGT')
    mock_aligner.align.side_effect = [
        [_make_aln_result(1.0, ['A'], ['A'])],
        [_make_aln_result(100.0, ['A', 'C'], ['A', 'C'])],
    ]

    assert get_best_alignment(_make_query('ACGT'), [low, high]).reference is high


def test_get_best_alignment_strips_gaps_from_the_query(mock_aligner) -> None:
    mock_aligner.align.return_value = [_make_aln_result(10.0, ['A'], ['A'])]

    get_best_alignment(_make_query(f'AC{GAP_SYMBOL}GT'), [_make_reference('PB2_ref_A')])

    mock_aligner.align.assert_called_once_with('ACGT', 'ACGT')
