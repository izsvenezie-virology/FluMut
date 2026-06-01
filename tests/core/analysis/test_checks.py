from unittest.mock import MagicMock, patch

import pytest

from flumut.core.analysis.checks import (
    Check,
    DuplicationCheck,
    EnlongationCheck,
    FrameshiftCheck,
    TruncationCheck,
    check_duplications,
    check_enlongation,
    check_frameshifts,
    check_truncation,
    perform_checks,
)
from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_alignment(protein_name: str, query: list[str], frameshifts: list | None = None) -> MagicMock:
    a = MagicMock()
    a.protein.name = protein_name
    a.alignment.query = query
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
    check_truncation(_make_analysis(sample))
    assert sample.checks == []


def test_check_truncation_stop_codon_in_middle_appends_check() -> None:
    # query[1] = '*', position is 1-based → 2
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', STOP_CODON_SYMBOL, 'V'])])
    check_truncation(_make_analysis(sample))
    assert len(sample.checks) == 1
    assert isinstance(sample.checks[0], TruncationCheck)


def test_check_truncation_position_is_one_based() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', STOP_CODON_SYMBOL, 'V'])])
    check_truncation(_make_analysis(sample))
    assert '2' in sample.checks[0].message


def test_check_truncation_stop_codon_at_last_position_no_check() -> None:
    # range(len - 1) never reaches the last index — last stop = elongation, not truncation
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K', STOP_CODON_SYMBOL])])
    check_truncation(_make_analysis(sample))
    assert sample.checks == []


def test_check_truncation_multiple_stop_codons_multiple_checks() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', [STOP_CODON_SYMBOL, STOP_CODON_SYMBOL, 'V'])])
    check_truncation(_make_analysis(sample))
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
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K', last_aa])])
    check_enlongation(_make_analysis(sample))
    assert sample.checks == []


def test_check_enlongation_regular_last_aa_appends_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K', 'V'])])
    check_enlongation(_make_analysis(sample))
    assert len(sample.checks) == 1
    assert isinstance(sample.checks[0], EnlongationCheck)


# ---------------------------------------------------------------------------
# check_frameshifts tests
# ---------------------------------------------------------------------------


def test_check_frameshifts_no_cds_no_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K'])])
    check_frameshifts(_make_analysis(sample))
    assert sample.checks == []


def test_check_frameshifts_cds_no_frameshifts_no_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K'], frameshifts=[])])
    check_frameshifts(_make_analysis(sample))
    assert sample.checks == []


def test_check_frameshifts_single_frameshift_appends_check() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K'], frameshifts=[(5, 10)])])
    check_frameshifts(_make_analysis(sample))
    assert len(sample.checks) == 1
    assert isinstance(sample.checks[0], FrameshiftCheck)


def test_check_frameshifts_multiple_frameshifts_multiple_checks() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M'], frameshifts=[(1, 3), (7, 9)])])
    check_frameshifts(_make_analysis(sample))
    assert len(sample.checks) == 2


def test_check_frameshifts_frameshift_with_none_end() -> None:
    sample = _make_sample('s1', [_make_alignment('PB2', ['M'], frameshifts=[(5, None)])])
    check_frameshifts(_make_analysis(sample))
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
    check_duplications(_make_analysis(sample))
    assert sample.checks == []


def test_check_duplications_repeated_protein_appends_check() -> None:
    sample = _make_sample(
        's1',
        [
            _make_alignment('PB2', ['M']),
            _make_alignment('PB2', ['M']),
        ],
    )
    check_duplications(_make_analysis(sample))
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
    check_duplications(_make_analysis(sample))
    assert len(sample.checks) == 2


# ---------------------------------------------------------------------------
# perform_checks tests
# ---------------------------------------------------------------------------


@patch('flumut.core.analysis.checks.check_duplications')
@patch('flumut.core.analysis.checks.check_frameshifts')
@patch('flumut.core.analysis.checks.check_enlongation')
@patch('flumut.core.analysis.checks.check_truncation')
def test_perform_checks_calls_all(mock_truncation, mock_enlongation, mock_frameshifts, mock_duplications) -> None:
    analysis = MagicMock()
    perform_checks(analysis)
    mock_truncation.assert_called_once_with(analysis)
    mock_enlongation.assert_called_once_with(analysis)
    mock_frameshifts.assert_called_once_with(analysis)
    mock_duplications.assert_called_once_with(analysis)
