from unittest.mock import MagicMock

import pytest

from flumut.core.analysis.models import Analysis, MarkerScan, PositionScan, ProteinAlignment, Sample
from flumut.core.nucleotides.models import Alignment
from flumut.flumutdb.models import MutationType

# ---------------------------------------------------------------------------
# PositionScan helpers
# ---------------------------------------------------------------------------


def _make_mapping(alteration: str, mutation_type: str = MutationType.SNP.value) -> MagicMock:
    mutation = MagicMock()
    mutation.type = mutation_type

    mapping = MagicMock()
    mapping.mutation = mutation
    mapping.alteration = alteration
    return mapping


# ---------------------------------------------------------------------------
# PositionScan tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'alteration, ammino_acid, expected',
    [
        ('A', 'A', True),
        ('A', 'G', False),
        ('A', 'AM', True),
        ('A', 'MV', False),
        ('M', 'AM', True),
        ('A', '-', False),
        ('A', '', False),
    ],
)
def test_snp_detection(alteration: str, ammino_acid: str, expected: bool) -> None:
    mapping = _make_mapping(alteration=alteration)
    ps = PositionScan(mapping=mapping, ammino_acid=ammino_acid)
    assert ps.is_detected is expected


def test_unsupported_mutation_type_raises() -> None:
    mapping = _make_mapping(alteration='A', mutation_type='INDEL')
    with pytest.raises(NotImplementedError, match='INDEL'):
        PositionScan(mapping=mapping, ammino_acid='A')


def test_mutation_property_returns_mapping_mutation() -> None:
    mapping = _make_mapping(alteration='A')
    ps = PositionScan(mapping=mapping, ammino_acid='A')
    assert ps.mutation is mapping.mutation


def test_is_detected_is_not_an_init_parameter() -> None:
    mapping = _make_mapping(alteration='A')
    with pytest.raises(TypeError):
        PositionScan(mapping=mapping, ammino_acid='A', is_detected=True)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# MarkerScan helpers
# ---------------------------------------------------------------------------


def _make_position_scan(is_detected: bool) -> PositionScan:
    alteration = 'A'
    ammino_acid = 'A' if is_detected else 'G'
    return PositionScan(mapping=_make_mapping(alteration=alteration), ammino_acid=ammino_acid)


def _make_marker(num_mutations: int) -> MagicMock:
    marker = MagicMock()
    marker.mutations = [MagicMock() for _ in range(num_mutations)]
    return marker


# ---------------------------------------------------------------------------
# MarkerScan tests
# ---------------------------------------------------------------------------


def test_markerscan_all_detected() -> None:
    marker = _make_marker(2)
    positions = [_make_position_scan(True), _make_position_scan(True)]
    ms = MarkerScan(marker=marker, positions=positions)
    assert ms.is_detected is True
    assert ms.is_complete is True
    assert ms.detected_mutations == positions


def test_markerscan_none_detected() -> None:
    marker = _make_marker(2)
    positions = [_make_position_scan(False), _make_position_scan(False)]
    ms = MarkerScan(marker=marker, positions=positions)
    assert ms.is_detected is False
    assert ms.is_complete is False
    assert ms.detected_mutations == []


def test_markerscan_partial_detection() -> None:
    marker = _make_marker(2)
    detected = _make_position_scan(True)
    not_detected = _make_position_scan(False)
    ms = MarkerScan(marker=marker, positions=[detected, not_detected])
    assert ms.is_detected is True
    assert ms.is_complete is False
    assert ms.detected_mutations == [detected]


def test_markerscan_empty_positions_and_mutations() -> None:
    ms = MarkerScan(marker=_make_marker(0), positions=[])
    assert ms.is_detected is False
    assert ms.is_complete is True  # 0 detected == 0 required
    assert ms.detected_mutations == []


def test_markerscan_detected_mutations_excludes_undetected() -> None:
    marker = _make_marker(3)
    p1 = _make_position_scan(True)
    p2 = _make_position_scan(False)
    p3 = _make_position_scan(True)
    ms = MarkerScan(marker=marker, positions=[p1, p2, p3])
    assert p1 in ms.detected_mutations
    assert p2 not in ms.detected_mutations
    assert p3 in ms.detected_mutations


# ---------------------------------------------------------------------------
# ProteinAlignment tests
# ---------------------------------------------------------------------------


def test_proteinalignment_stores_fields() -> None:
    protein = MagicMock()
    reference = MagicMock()
    alignment = Alignment(reference=['A', 'C'], query=['A', 'G'])
    pa = ProteinAlignment(protein=protein, reference=reference, alignment=alignment)
    assert pa.protein is protein
    assert pa.reference is reference
    assert pa.alignment is alignment


def test_proteinalignment_cds_defaults_to_none() -> None:
    pa = ProteinAlignment(protein=MagicMock(), reference=MagicMock(), alignment=Alignment())
    assert pa.cds is None


def test_proteinalignment_cds_can_be_provided() -> None:
    cds = MagicMock()
    pa = ProteinAlignment(protein=MagicMock(), reference=MagicMock(), alignment=Alignment(), cds=cds)
    assert pa.cds is cds


# ---------------------------------------------------------------------------
# Sample tests
# ---------------------------------------------------------------------------


def test_sample_stores_id() -> None:
    s = Sample(id='seq1')
    assert s.id == 'seq1'


def test_sample_alignments_defaults_to_empty_list() -> None:
    s = Sample(id='seq1')
    assert s.alignments == []


def test_sample_alignments_can_be_provided() -> None:
    alignment = MagicMock(spec=ProteinAlignment)
    s = Sample(id='seq1', alignments=[alignment])
    assert s.alignments == [alignment]


def test_sample_positions_starts_empty() -> None:
    s = Sample(id='seq1')
    assert s.positions == []


def test_sample_marker_scans_starts_empty() -> None:
    s = Sample(id='seq1')
    assert s.marker_scans == []


def test_sample_positions_is_not_an_init_parameter() -> None:
    with pytest.raises(TypeError):
        Sample(id='seq1', positions=[])  # type: ignore[call-arg]


def test_sample_marker_scans_is_not_an_init_parameter() -> None:
    with pytest.raises(TypeError):
        Sample(id='seq1', marker_scans=[])  # type: ignore[call-arg]


def test_sample_positions_list_is_independent_per_instance() -> None:
    s1 = Sample(id='seq1')
    s2 = Sample(id='seq2')
    s1.positions.append(MagicMock())
    assert s2.positions == []


def test_sample_marker_scans_list_is_independent_per_instance() -> None:
    s1 = Sample(id='seq1')
    s2 = Sample(id='seq2')
    s1.marker_scans.append(MagicMock())
    assert s2.marker_scans == []


def test_sample_checks_starts_empty() -> None:
    assert Sample(id='seq1').checks == []


def test_sample_checks_list_is_independent_per_instance() -> None:
    s1 = Sample(id='seq1')
    s2 = Sample(id='seq2')
    s1.checks.append(MagicMock())
    assert s2.checks == []


# ---------------------------------------------------------------------------
# Analysis tests
# ---------------------------------------------------------------------------


def test_analysis_defaults_all_empty() -> None:
    a = Analysis()
    assert a.samples == {}
    assert a.mutations == set()
    assert a.markers == set()
    assert a.literature == set()


def test_analysis_samples_dict_is_independent_per_instance() -> None:
    a1 = Analysis()
    a2 = Analysis()
    a1.samples['seq1'] = Sample(id='seq1')
    assert a2.samples == {}


def test_analysis_mutations_set_is_independent_per_instance() -> None:
    a1 = Analysis()
    a2 = Analysis()
    a1.mutations.add(MagicMock())
    assert a2.mutations == set()
