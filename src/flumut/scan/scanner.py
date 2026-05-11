from flumut.analysis import ProteinAlignment
from flumut.flumutdb import Marker
from flumut.scan import MarkerScan, PositionScan


def scan_positions(alignment: ProteinAlignment) -> list[PositionScan]:
    positions = alignment.get_positions()
    result = []

    for mutation in alignment.protein.mutations:
        for mapping in mutation.mappings:
            if not mapping.reference == alignment.reference:
                continue
            index = positions.index(mapping.position)
            aa = alignment.aligned_query[index]
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
