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


def test_analyse_twice_does_not_duplicate_positions(mock_scanner) -> None:
    """Per-sample results are recomputed like the analysis-level ones, not appended to."""
    analysis = _make_analysis(num_samples=2, alignments_per_sample=3)

    analyse(analysis)
    analyse(analysis)

    for sample in analysis.samples.values():
        assert len(sample.positions) == 3


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


