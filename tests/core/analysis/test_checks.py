from unittest.mock import MagicMock, patch

import pytest

from flumut.core.analysis.checks import (
    check_duplications,
    check_enlongation,
    check_frameshifts,
    check_truncation,
    perform_checks,
)
from flumut.core.analysis.models import Check, DuplicationCheck, EnlongationCheck, FrameshiftCheck, TruncationCheck
from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_alignment(
    protein_name: str,
    query: list[str],
    reference: list[str] | None = None,
    frameshifts: list | None = None,
) -> MagicMock:
    a = MagicMock()
    a.protein.name = protein_name
    a.alignment.query = query
    if reference is not None:
        a.alignment.reference = reference
    if frameshifts is not None:
        a.cds.frameshifts = frameshifts
    else:
        a.cds = None
    return a


def _make_sample(sample_id: str, alignments: list) -> MagicMock:
    s = MagicMock()
    s.id = sample_id
    s.checks = []
    s.alignments = alignments
    return s


def _make_analysis(*samples) -> MagicMock:
    analysis = MagicMock()
    analysis.samples = {s.id: s for s in samples}
    return analysis


# ---------------------------------------------------------------------------
# Check base class and subclass message tests
# ---------------------------------------------------------------------------


def test_check_stores_message() -> None:
    c = Check('test message')
    assert c.message == 'test message'


def test_truncation_check_message() -> None:
    c = TruncationCheck('sampleA', 'PB2', 42)
    assert 'sampleA' in c.message
    assert 'PB2' in c.message
    assert '42' in c.message


def test_enlongation_check_message() -> None:
    c = EnlongationCheck('sampleA', 'PB2')
    assert 'sampleA' in c.message
    assert 'PB2' in c.message


def test_duplication_check_message() -> None:
    c = DuplicationCheck('sampleA', 'PB2')
    assert 'sampleA' in c.message
    assert 'PB2' in c.message


def test_frameshift_check_message_with_end_position() -> None:
    c = FrameshiftCheck('sampleA', 'PB2', 10, 20)
    assert 'sampleA' in c.message
    assert 'PB2' in c.message
    assert '10' in c.message
    assert '20' in c.message


def test_frameshift_check_message_with_no_end_uses_end_literal() -> None:
    c = FrameshiftCheck('sampleA', 'PB2', 10, None)
    assert '10' in c.message
    assert 'end' in c.message


# ---------------------------------------------------------------------------
# check_truncation tests
# ---------------------------------------------------------------------------


def test_check_truncation_no_stop_codon_appends_no_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K', 'V'])])
    check_truncation(sample)
    assert sample.checks == []


def test_check_truncation_stop_codon_in_middle_appends_check() -> None:
    # query[1] = '*', reference[1] != '*' → truncation at 1-based position 2
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', STOP_CODON_SYMBOL, 'V'], reference=['A', 'K', 'V'])])
    check_truncation(sample)
    assert len(sample.checks) == 1
    assert isinstance(sample.checks[0], TruncationCheck)


def test_check_truncation_position_is_one_based() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', STOP_CODON_SYMBOL, 'V'], reference=['A', 'K', 'V'])])
    check_truncation(sample)
    assert '2' in sample.checks[0].message


def test_check_truncation_no_check_when_reference_has_stop_at_same_position() -> None:
    # reference also has stop at that position = correct termination, not truncation
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K', STOP_CODON_SYMBOL], reference=['A', 'K', STOP_CODON_SYMBOL])])
    check_truncation(sample)
    assert sample.checks == []


def test_check_truncation_multiple_stop_codons_multiple_checks() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', [STOP_CODON_SYMBOL, STOP_CODON_SYMBOL, 'V'], reference=['K', 'K', 'V'])])
    check_truncation(sample)
    assert len(sample.checks) == 2


# ---------------------------------------------------------------------------
# check_enlongation tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'last_aa',
    [STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL, GAP_SYMBOL],
    ids=['stop_codon', 'unknown_aa', 'gap'],
)
def test_check_enlongation_no_check_when_last_aa_is_terminal(last_aa: str) -> None:
    # reference ends with stop codon — the query guard is in play
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K', last_aa], reference=['A', STOP_CODON_SYMBOL])])
    check_enlongation(sample)
    assert sample.checks == []


def test_check_enlongation_regular_last_aa_appends_check() -> None:
    # reference ends with stop codon — query is elongated
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K', 'V'], reference=['A', STOP_CODON_SYMBOL])])
    check_enlongation(sample)
    assert len(sample.checks) == 1
    assert isinstance(sample.checks[0], EnlongationCheck)


def test_check_enlongation_no_check_when_reference_does_not_end_with_stop_query_stop() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', STOP_CODON_SYMBOL], reference=['A', 'K'])])
    check_enlongation(sample)
    assert sample.checks == []


def test_check_enlongation_no_check_when_reference_does_not_end_with_stop_query_regular() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'V'], reference=['A', 'K'])])
    check_enlongation(sample)
    assert sample.checks == []


# ---------------------------------------------------------------------------
# check_frameshifts tests
# ---------------------------------------------------------------------------


def test_check_frameshifts_no_cds_no_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K'])])
    check_frameshifts(sample)
    assert sample.checks == []


def test_check_frameshifts_cds_no_frameshifts_no_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K'], frameshifts=[])])
    check_frameshifts(sample)
    assert sample.checks == []


def test_check_frameshifts_single_frameshift_appends_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K'], frameshifts=[(5, 10)])])
    check_frameshifts(sample)
    assert len(sample.checks) == 1
    assert isinstance(sample.checks[0], FrameshiftCheck)


def test_check_frameshifts_multiple_frameshifts_multiple_checks() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M'], frameshifts=[(1, 3), (7, 9)])])
    check_frameshifts(sample)
    assert len(sample.checks) == 2


def test_check_frameshifts_frameshift_with_none_end() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M'], frameshifts=[(5, None)])])
    check_frameshifts(sample)
    assert len(sample.checks) == 1
    assert 'end' in sample.checks[0].message


# ---------------------------------------------------------------------------
# check_duplications tests
# ---------------------------------------------------------------------------


def test_check_duplications_unique_proteins_no_check() -> None:
    sample = _make_sample(
        's1',
        [
            _make_alignment('PB2', ['M']),
            _make_alignment('PB1', ['M']),
        ],
    )
    check_duplications(sample)
    assert sample.checks == []


def test_check_duplications_repeated_protein_appends_check() -> None:
    sample = _make_sample(
        's1',
        [
            _make_alignment('PB2', ['M']),
            _make_alignment('PB2', ['M']),
        ],
    )
    check_duplications(sample)
    assert len(sample.checks) == 1
    assert isinstance(sample.checks[0], DuplicationCheck)


def test_check_duplications_three_occurrences_two_checks() -> None:
    sample = _make_sample(
        's1',
        [
            _make_alignment('PB2', ['M']),
            _make_alignment('PB2', ['M']),
            _make_alignment('PB2', ['M']),
        ],
    )
    check_duplications(sample)
    assert len(sample.checks) == 2


# ---------------------------------------------------------------------------
# perform_checks tests
# ---------------------------------------------------------------------------


@patch('flumut.core.analysis.checks.check_duplications')
@patch('flumut.core.analysis.checks.check_frameshifts')
@patch('flumut.core.analysis.checks.check_enlongation')
@patch('flumut.core.analysis.checks.check_truncation')
def test_perform_checks_calls_all_four_per_sample(mock_truncation, mock_enlongation, mock_frameshifts, mock_duplications) -> None:
    sample = MagicMock()
    analysis = _make_analysis(sample)
    perform_checks(analysis)
    mock_truncation.assert_called_once_with(sample)
    mock_enlongation.assert_called_once_with(sample)
    mock_frameshifts.assert_called_once_with(sample)
    mock_duplications.assert_called_once_with(sample)
