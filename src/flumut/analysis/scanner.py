from flumut.analysis.models import Analysis, MarkerScan, PositionScan
from flumut.core.models import ProteinAlignment
from flumut.flumutdb import Marker


def analyse(analysis: Analysis, relaxed: bool = True) -> None:
    analysis.mutations.clear()
    analysis.markers.clear()
    analysis.literature.clear()

    for sample in analysis.samples.values():
        for alignment in sample.alignments:
            sample.positions += scan_positions(alignment)
        sample.marker_scans = scan_markers(sample.positions, relaxed)

        for position in sample.positions:
            if position.is_detected:
                analysis.mutations.add(position.mutation)

        for scan in sample.marker_scans:
            analysis.markers.add(scan.marker)

    for marker in analysis.markers:
        for evidence in marker.evidences:
            analysis.literature.add(evidence.paper)


def scan_positions(alignment: ProteinAlignment) -> list[PositionScan]:
    positions = alignment.alignment.get_positions()
    result = []

    for mutation in alignment.protein.mutations:
        for mapping in mutation.mappings:
            if not mapping.reference == alignment.reference:
                continue
            index = positions.index(mapping.position)
            aa = alignment.alignment.query[index]
            result.append(PositionScan(mapping=mapping, ammino_acid=aa))
    return result


def scan_markers(positions: list[PositionScan], relaxed: bool) -> list[MarkerScan]:
    markers = []
    mapping = {pos.mutation: pos for pos in positions}

    for marker in Marker.all():
        marker_positions = []
        for mutation in marker.mutations:
            pos = mapping.get(mutation, None)
            if pos:
                marker_positions.append(pos)

        ms = MarkerScan(marker, marker_positions)
        if not ms.is_detected:
            continue
        if not relaxed and not ms.is_complete:
            continue
        markers.append(ms)
    return markers
