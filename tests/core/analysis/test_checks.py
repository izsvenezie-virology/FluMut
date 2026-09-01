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

STOP, UNKNOWN, GAP = STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL, GAP_SYMBOL

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


# ---------------------------------------------------------------------------
# Check messages
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'check_class, args, fragments',
    [
        (Check, ('test message',), ['test message']),
        (TruncationCheck, ('sampleA', 'PB2', 42), ['sampleA', 'PB2', '42']),
        (EnlongationCheck, ('sampleA', 'PB2'), ['sampleA', 'PB2']),
        (DuplicationCheck, ('sampleA', 'PB2'), ['sampleA', 'PB2']),
        (FrameshiftCheck, ('sampleA', 'PB2', 10, 20), ['sampleA', 'PB2', '10', '20']),
        (FrameshiftCheck, ('sampleA', 'PB2', 10, None), ['sampleA', 'PB2', '10', 'end']),
    ],
    ids=['base', 'truncation', 'enlongation', 'duplication', 'frameshift', 'frameshift_open_ended'],
)
def test_check_message_identifies_the_problem(check_class, args: tuple, fragments: list[str]) -> None:
    """Every check message names the sample, the protein and the affected positions."""
    message = check_class(*args).message
    for fragment in fragments:
        assert fragment in message


# ---------------------------------------------------------------------------
# check_truncation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'query, reference, expected_positions',
    [
        (['M', 'K', 'V'], ['A', 'K', 'V'], []),
        (['M', STOP, 'V'], ['A', 'K', 'V'], [2]),
        (['M', 'K', STOP], ['A', 'K', STOP], []),
        ([STOP, STOP, 'V'], ['K', 'K', 'V'], [1, 2]),
    ],
    ids=['no_stop', 'premature_stop', 'stop_matches_reference', 'two_premature_stops'],
)
def test_check_truncation(query: list[str], reference: list[str], expected_positions: list[int]) -> None:
    """A stop codon is a truncation only where the reference has none. Positions are 1-based."""
    # A protein name without digits, so the position assertion cannot pass on the name alone.
    sample = _make_sample('s1', [_make_alignment('NP', query, reference=reference)])
    check_truncation(sample)

    assert len(sample.checks) == len(expected_positions)
    for check, position in zip(sample.checks, expected_positions):
        assert isinstance(check, TruncationCheck)
        assert check.message.endswith(f'at {position}')


# ---------------------------------------------------------------------------
# check_enlongation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'query, reference, expected',
    [
        (['M', 'K', STOP], ['A', STOP], 0),
        (['M', 'K', UNKNOWN], ['A', STOP], 0),
        (['M', 'K', GAP], ['A', STOP], 0),
        (['M', 'K', 'V'], ['A', STOP], 1),
        (['M', STOP], ['A', 'K'], 0),
        (['M', 'V'], ['A', 'K'], 0),
    ],
    ids=[
        'terminates_with_stop',
        'terminates_with_unknown',
        'terminates_with_gap',
        'reads_through',
        'reference_has_no_stop_query_stops',
        'reference_has_no_stop_query_continues',
    ],
)
def test_check_enlongation(query: list[str], reference: list[str], expected: int) -> None:
    """Elongation is reported only when the reference terminates and the query does not."""
    sample = _make_sample('s1', [_make_alignment('PB2', query, reference=reference)])
    check_enlongation(sample)

    assert len(sample.checks) == expected
    assert all(isinstance(check, EnlongationCheck) for check in sample.checks)


# ---------------------------------------------------------------------------
# check_frameshifts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'frameshifts, expected, fragment',
    [
        (None, 0, None),
        ([], 0, None),
        ([(5, 10)], 1, '10'),
        ([(1, 3), (7, 9)], 2, None),
        ([(5, None)], 1, 'end'),
    ],
    ids=['no_cds', 'cds_without_frameshifts', 'single', 'multiple', 'open_ended'],
)
def test_check_frameshifts(frameshifts: list | None, expected: int, fragment: str | None) -> None:
    """One check per frameshift range recorded on the CDS; none when there is no CDS."""
    sample = _make_sample('s1', [_make_alignment('PB2', ['M', 'K'], frameshifts=frameshifts)])
    check_frameshifts(sample)

    assert len(sample.checks) == expected
    assert all(isinstance(check, FrameshiftCheck) for check in sample.checks)
    if fragment:
        assert fragment in sample.checks[0].message


# ---------------------------------------------------------------------------
# check_duplications
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'protein_names, expected',
    [
        (['PB2', 'PB1'], 0),
        (['PB2', 'PB2'], 1),
        (['PB2', 'PB2', 'PB2'], 2),
    ],
    ids=['unique', 'one_repeat', 'two_repeats'],
)
def test_check_duplications(protein_names: list[str], expected: int) -> None:
    """Every occurrence of a protein after the first is reported."""
    sample = _make_sample('s1', [_make_alignment(name, ['M']) for name in protein_names])
    check_duplications(sample)

    assert len(sample.checks) == expected
    assert all(isinstance(check, DuplicationCheck) for check in sample.checks)


# ---------------------------------------------------------------------------
# perform_checks
# ---------------------------------------------------------------------------


@patch('flumut.core.analysis.checks.check_duplications')
@patch('flumut.core.analysis.checks.check_frameshifts')
@patch('flumut.core.analysis.checks.check_enlongation')
@patch('flumut.core.analysis.checks.check_truncation')
def test_perform_checks_runs_every_check_on_every_sample(
    mock_truncation, mock_enlongation, mock_frameshifts, mock_duplications
) -> None:
    sample = MagicMock()
    analysis = MagicMock()
    analysis.samples = {sample.id: sample}

    perform_checks(analysis)

    mock_truncation.assert_called_once_with(sample)
    mock_enlongation.assert_called_once_with(sample)
    mock_frameshifts.assert_called_once_with(sample)
    mock_duplications.assert_called_once_with(sample)
