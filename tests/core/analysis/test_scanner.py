from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from flumut.core.analysis.models import Analysis, ProteinAlignment, Sample
from flumut.core.analysis.scanner import analyse, scan_positions
from flumut.core.globals import GAP_SYMBOL
from flumut.core.nucleotides.models import Alignment
from flumut.flumutdb.models import MutationType

# ---------------------------------------------------------------------------
# scan_positions
# ---------------------------------------------------------------------------


def _make_mapping(position: int, alteration: str) -> MagicMock:
    mapping = MagicMock()
    mapping.position = position
    mapping.alteration = alteration
    mapping.mutation.type = MutationType.SNP.value
    return mapping


def _make_protein_alignment(reference: str, query: str, mappings: list) -> ProteinAlignment:
    reference_model = MagicMock()
    protein = MagicMock()
    # scan_positions reads the mappings pre-grouped per (reference, protein).
    reference_model.mappings_by_protein = {protein: mappings}
    return ProteinAlignment(
        protein=protein,
        reference=reference_model,
        alignment=Alignment(reference=list(reference), query=list(query)),
    )


@pytest.mark.parametrize(
    'reference, query, positions, expected',
    [
        ('ACT', 'KMR', (), ''),
        ('ACT', 'KMR', ((2, 'M'),), 'M'),
        ('ACT', 'KMR', ((1, 'K'), (3, 'R')), 'KR'),
        # A gap column is not numbered, so reference position 2 lands on index 2.
        (f'{GAP_SYMBOL}AC', 'XKR', ((2, 'R'),), 'R'),
    ],
    ids=['no_mappings', 'single_mapping', 'two_mappings', 'gap_does_not_shift_positions'],
)
def test_scan_positions(reference: str, query: str, positions: tuple, expected: str) -> None:
    """Each mapping is resolved to the residue at its 1-based reference position."""
    mappings = [_make_mapping(position, alteration) for position, alteration in positions]
    result = scan_positions(_make_protein_alignment(reference, query, mappings))

    assert [scan.mapping for scan in result] == mappings
    assert ''.join(scan.ammino_acid for scan in result) == expected


# ---------------------------------------------------------------------------
# analyse
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_scanner():
    with (
        patch('flumut.core.analysis.scanner.scan_positions') as mock_scan_positions,
        patch('flumut.core.analysis.scanner.scan_markers') as mock_scan_markers,
    ):
        mock_scan_positions.side_effect = lambda alignment: [MagicMock(is_detected=False)]
        mock_scan_markers.return_value = []
        yield SimpleNamespace(scan_positions=mock_scan_positions, scan_markers=mock_scan_markers)


def _make_analysis(num_samples: int, alignments_per_sample: int) -> Analysis:
    samples = {
        f's{i}': Sample(id=f's{i}', alignments=[MagicMock() for _ in range(alignments_per_sample)])
        for i in range(num_samples)
    }
    return Analysis(samples=samples)


def test_analyse_clears_results_from_a_previous_run(mock_scanner) -> None:
    analysis = Analysis(mutations={MagicMock()}, markers={MagicMock()}, literature={MagicMock()})
    analyse(analysis)
    assert (analysis.mutations, analysis.markers, analysis.literature) == (set(), set(), set())


@pytest.mark.parametrize(
    'num_samples, alignments_per_sample',
    [(0, 0), (1, 2), (2, 1), (2, 3)],
    ids=['no_samples', 'one_sample_two_alignments', 'two_samples', 'two_samples_three_alignments'],
)
def test_analyse_scans_every_alignment_of_every_sample(mock_scanner, num_samples: int, alignments_per_sample: int) -> None:
    analysis = _make_analysis(num_samples, alignments_per_sample)

    analyse(analysis)

    assert mock_scanner.scan_positions.call_count == num_samples * alignments_per_sample
    assert mock_scanner.scan_markers.call_count == num_samples
    for sample in analysis.samples.values():
        assert len(sample.positions) == alignments_per_sample


def test_analyse_collects_mutations_markers_and_literature(mock_scanner) -> None:
    """Only detected positions become mutations; every reported marker contributes its papers."""
    mutation = MagicMock()
    detected = MagicMock(is_detected=True, mutation=mutation)
    undetected = MagicMock(is_detected=False)
    paper = MagicMock()
    marker = MagicMock(evidences=[MagicMock(paper=paper)])
    marker_scan = MagicMock(marker=marker)

    mock_scanner.scan_positions.side_effect = None
    mock_scanner.scan_positions.return_value = [detected, undetected]
    mock_scanner.scan_markers.return_value = [marker_scan]

    sample = Sample(id='s1', alignments=[MagicMock()])
    analysis = Analysis(samples={'s1': sample})
    analyse(analysis)

    assert analysis.mutations == {mutation}
    assert analysis.markers == {marker}
    assert analysis.literature == {paper}
    assert sample.marker_scans == [marker_scan]


@pytest.mark.parametrize(
    'kwargs, expected',
    [({}, True), ({'relaxed': True}, True), ({'relaxed': False}, False)],
    ids=['default', 'explicit_true', 'explicit_false'],
)
def test_analyse_passes_relaxed_through_to_scan_markers(mock_scanner, kwargs: dict, expected: bool) -> None:
    analysis = Analysis(samples={'s1': Sample(id='s1', alignments=[MagicMock()])})

    analyse(analysis, **kwargs)

    assert mock_scanner.scan_markers.call_args.args[1] is expected
