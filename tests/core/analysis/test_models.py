from unittest.mock import MagicMock

import pytest

from flumut.core.analysis.models import MarkerScan, PositionScan
from flumut.flumutdb.models import MutationType

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_mapping(alteration: str, mutation_type: str = MutationType.SNP.value) -> MagicMock:
    mutation = MagicMock()
    mutation.type = mutation_type

    mapping = MagicMock()
    mapping.mutation = mutation
    mapping.alteration = alteration
    return mapping


def _make_position_scan(is_detected: bool) -> PositionScan:
    return PositionScan(mapping=_make_mapping(alteration='A'), ammino_acid='A' if is_detected else 'G')


def _make_marker(num_mutations: int) -> MagicMock:
    marker = MagicMock()
    marker.mutations = [MagicMock() for _ in range(num_mutations)]
    return marker


# ---------------------------------------------------------------------------
# PositionScan
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
    """A SNP is detected when the expected alteration appears in the observed residues."""
    ps = PositionScan(mapping=_make_mapping(alteration=alteration), ammino_acid=ammino_acid)
    assert ps.is_detected is expected


def test_unsupported_mutation_type_raises() -> None:
    mapping = _make_mapping(alteration='A', mutation_type='INDEL')
    with pytest.raises(NotImplementedError, match='INDEL'):
        PositionScan(mapping=mapping, ammino_acid='A')


def test_mutation_property_returns_mapping_mutation() -> None:
    mapping = _make_mapping(alteration='A')
    assert PositionScan(mapping=mapping, ammino_acid='A').mutation is mapping.mutation


# ---------------------------------------------------------------------------
# MarkerScan
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'num_mutations, detections, expected_detected, expected_complete',
    [
        (2, [True, True], True, True),
        (2, [False, False], False, False),
        (2, [True, False], True, False),
        (0, [], False, True),
        (3, [True, False, True], True, False),
    ],
    ids=['all_detected', 'none_detected', 'partial', 'no_positions', 'mixed'],
)
def test_marker_scan_derives_detection_state(
    num_mutations: int, detections: list[bool], expected_detected: bool, expected_complete: bool
) -> None:
    """``is_detected``, ``is_complete`` and ``detected_mutations`` follow from the positions."""
    positions = [_make_position_scan(flag) for flag in detections]
    ms = MarkerScan(marker=_make_marker(num_mutations), positions=positions)

    assert ms.is_detected is expected_detected
    assert ms.is_complete is expected_complete
    assert ms.detected_mutations == [pos for pos, flag in zip(positions, detections) if flag]
