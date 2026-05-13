from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from flumut.core.analysis.models import Analysis, ProteinAlignment, Sample
from flumut.core.analysis.scanner import analyse, scan_positions
from flumut.core.globals import GAP_SYMBOL
from flumut.core.nucleotides.models import Alignment
from flumut.flumutdb.models import MutationType


def _make_mapping(position: int, alteration: str, reference: MagicMock) -> MagicMock:
    mapping = MagicMock()
    mapping.reference = reference
    mapping.position = position
    mapping.alteration = alteration
    mapping.mutation.type = MutationType.SNP.value
    return mapping


def _make_protein_alignment(
    reference_chars: list[str],
    query_chars: list[str],
    mutations: list | None = None,
    reference: MagicMock | None = None,
) -> ProteinAlignment:
    if reference is None:
        reference = MagicMock()
    protein = MagicMock()
    protein.mutations = mutations or []
    return ProteinAlignment(
        protein=protein,
        reference=reference,
        alignment=Alignment(reference=reference_chars, query=query_chars),
    )


# ---------------------------------------------------------------------------
# scan_positions tests
# ---------------------------------------------------------------------------


def test_scan_positions_no_mutations_returns_empty() -> None:
    pa = _make_protein_alignment(['A', 'C', 'T'], ['K', 'M', 'R'])
    assert scan_positions(pa) == []


def test_scan_positions_matching_reference_returns_position_scan() -> None:
    reference = MagicMock()
    # positions: [1, 2, 3] → position 2 is at index 1 → aa 'M'
    mapping = _make_mapping(position=2, alteration='M', reference=reference)
    mutation = MagicMock(mappings=[mapping])
    pa = _make_protein_alignment(
        reference_chars=['A', 'C', 'T'],
        query_chars=['K', 'M', 'R'],
        mutations=[mutation],
        reference=reference,
    )
    result = scan_positions(pa)
    assert len(result) == 1
    assert result[0].mapping is mapping
    assert result[0].ammino_acid == 'M'


def test_scan_positions_nonmatching_reference_excluded() -> None:
    mapping = _make_mapping(position=1, alteration='K', reference=MagicMock())
    mutation = MagicMock(mappings=[mapping])
    # pa.reference is a different object than mapping.reference
    pa = _make_protein_alignment(
        reference_chars=['A'],
        query_chars=['K'],
        mutations=[mutation],
    )
    assert scan_positions(pa) == []


def test_scan_positions_multiple_mutations() -> None:
    reference = MagicMock()
    mapping1 = _make_mapping(position=1, alteration='K', reference=reference)
    mapping2 = _make_mapping(position=3, alteration='R', reference=reference)
    pa = _make_protein_alignment(
        reference_chars=['A', 'C', 'T'],
        query_chars=['K', 'M', 'R'],
        mutations=[MagicMock(mappings=[mapping1]), MagicMock(mappings=[mapping2])],
        reference=reference,
    )
    result = scan_positions(pa)
    assert len(result) == 2
    assert result[0].mapping is mapping1
    assert result[0].ammino_acid == 'K'
    assert result[1].mapping is mapping2
    assert result[1].ammino_acid == 'R'


def test_scan_positions_mixed_references_only_matching_included() -> None:
    reference = MagicMock()
    matching = _make_mapping(position=1, alteration='K', reference=reference)
    excluded = _make_mapping(position=2, alteration='M', reference=MagicMock())
    mutation = MagicMock(mappings=[matching, excluded])
    pa = _make_protein_alignment(
        reference_chars=['A', 'C'],
        query_chars=['K', 'M'],
        mutations=[mutation],
        reference=reference,
    )
    result = scan_positions(pa)
    assert len(result) == 1
    assert result[0].mapping is matching


def test_scan_positions_gap_in_alignment_does_not_shift_positions() -> None:
    # reference: ['-', 'A', 'C'] → get_positions() = [None, 1, 2]
    # position 2 → index 2 → aa 'R' (not the character at index 1)
    reference = MagicMock()
    mapping = _make_mapping(position=2, alteration='R', reference=reference)
    pa = _make_protein_alignment(
        reference_chars=[GAP_SYMBOL, 'A', 'C'],
        query_chars=['X', 'K', 'R'],
        mutations=[MagicMock(mappings=[mapping])],
        reference=reference,
    )
    result = scan_positions(pa)
    assert len(result) == 1
    assert result[0].ammino_acid == 'R'


# ---------------------------------------------------------------------------
# analyse fixture and tests
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_scanner():
    with (
        patch('flumut.core.analysis.scanner.scan_positions') as mock_scan_positions,
        patch('flumut.core.analysis.scanner.scan_markers') as mock_scan_markers,
    ):
        yield SimpleNamespace(
            scan_positions=mock_scan_positions,
            scan_markers=mock_scan_markers,
        )


def test_analyse_clears_existing_state(mock_scanner) -> None:
    analysis = Analysis(
        mutations={MagicMock()},
        markers={MagicMock()},
        literature={MagicMock()},
    )
    mock_scanner.scan_positions.return_value = []
    mock_scanner.scan_markers.return_value = []
    analyse(analysis)
    assert analysis.mutations == set()
    assert analysis.markers == set()
    assert analysis.literature == set()


def test_analyse_no_samples_does_not_call_scan_functions(mock_scanner) -> None:
    analyse(Analysis())
    mock_scanner.scan_positions.assert_not_called()
    mock_scanner.scan_markers.assert_not_called()


def test_analyse_positions_accumulated_from_all_alignments(mock_scanner) -> None:
    sample = Sample(id='s1', alignments=[MagicMock(), MagicMock()])
    analysis = Analysis(samples={'s1': sample})
    pos1, pos2 = MagicMock(is_detected=False), MagicMock(is_detected=False)
    mock_scanner.scan_positions.side_effect = [[pos1], [pos2]]
    mock_scanner.scan_markers.return_value = []
    analyse(analysis)
    assert sample.positions == [pos1, pos2]
    assert mock_scanner.scan_positions.call_count == 2


def test_analyse_scan_markers_called_with_positions_and_relaxed(mock_scanner) -> None:
    sample = Sample(id='s1', alignments=[MagicMock()])
    analysis = Analysis(samples={'s1': sample})
    pos = MagicMock(is_detected=False)
    mock_scanner.scan_positions.return_value = [pos]
    mock_scanner.scan_markers.return_value = []
    analyse(analysis, relaxed=False)
    mock_scanner.scan_markers.assert_called_once_with([pos], False)


def test_analyse_relaxed_defaults_to_true(mock_scanner) -> None:
    sample = Sample(id='s1', alignments=[MagicMock()])
    analysis = Analysis(samples={'s1': sample})
    mock_scanner.scan_positions.return_value = []
    mock_scanner.scan_markers.return_value = []
    analyse(analysis)
    mock_scanner.scan_markers.assert_called_once_with([], True)


def test_analyse_only_detected_positions_added_to_mutations(mock_scanner) -> None:
    sample = Sample(id='s1', alignments=[MagicMock()])
    mutation = MagicMock()
    detected = MagicMock(is_detected=True, mutation=mutation)
    undetected = MagicMock(is_detected=False)
    mock_scanner.scan_positions.return_value = [detected, undetected]
    mock_scanner.scan_markers.return_value = []
    analysis = Analysis(samples={'s1': sample})
    analyse(analysis)
    assert analysis.mutations == {mutation}


def test_analyse_marker_scans_set_on_sample_and_markers_collected(mock_scanner) -> None:
    sample = Sample(id='s1', alignments=[MagicMock()])
    marker = MagicMock(evidences=[])
    ms = MagicMock(marker=marker)
    mock_scanner.scan_positions.return_value = []
    mock_scanner.scan_markers.return_value = [ms]
    analysis = Analysis(samples={'s1': sample})
    analyse(analysis)
    assert sample.marker_scans == [ms]
    assert marker in analysis.markers


def test_analyse_evidence_papers_collected_into_literature(mock_scanner) -> None:
    sample = Sample(id='s1', alignments=[MagicMock()])
    paper = MagicMock()
    marker = MagicMock(evidences=[MagicMock(paper=paper)])
    mock_scanner.scan_positions.return_value = []
    mock_scanner.scan_markers.return_value = [MagicMock(marker=marker)]
    analysis = Analysis(samples={'s1': sample})
    analyse(analysis)
    assert paper in analysis.literature


def test_analyse_multiple_samples_all_processed(mock_scanner) -> None:
    s1 = Sample(id='s1', alignments=[MagicMock()])
    s2 = Sample(id='s2', alignments=[MagicMock()])
    m1, m2 = MagicMock(), MagicMock()
    mock_scanner.scan_positions.side_effect = [
        [MagicMock(is_detected=True, mutation=m1)],
        [MagicMock(is_detected=True, mutation=m2)],
    ]
    mock_scanner.scan_markers.return_value = []
    analyse(Analysis(samples={'s1': s1, 's2': s2}))
    assert mock_scanner.scan_positions.call_count == 2
    assert mock_scanner.scan_markers.call_count == 2
