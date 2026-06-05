from flumut.core.analysis.models import Analysis, MarkerScan, PositionScan, ProteinAlignment
from flumut.core.logger import LOGGER
from flumut.flumutdb import Marker


def analyse(analysis: Analysis, relaxed: bool = True) -> None:
    """Scan all samples for mutations and markers, populating the analysis results.

    Clears and recomputes ``analysis.mutations``, ``analysis.markers``, and
    ``analysis.literature`` from the current sample alignments. For each sample,
    positions are scanned against the database and matching markers are collected.

    Args:
        analysis: The Analysis object to scan and update.
        relaxed: When True, report markers where at least one mutation is detected.
            When False, only report markers where every mutation is detected.
    """
    analysis.mutations.clear()
    analysis.markers.clear()
    analysis.literature.clear()

    for sample in analysis.samples.values():
        LOGGER.debug(f'Scanning {sample.id}')
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
    """Scan a protein alignment for all known mutation positions.

    For each mutation mapped to the alignment's reference, retrieves the
    observed amino acid at that position and returns a PositionScan result.

    Args:
        alignment: The protein alignment to scan.

    Returns:
        A list of PositionScan objects, one per mutation mapping that matches
        the alignment's reference.
    """
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
    """Evaluate all database markers against a set of scanned positions.

    Args:
        positions: PositionScan results for a single sample.
        relaxed: When True, include markers with at least one detected mutation.
            When False, include only markers where all mutations are detected.

    Returns:
        A list of MarkerScan objects for markers that pass the detection threshold.
    """
    markers = []
    mapping = {pos.mutation: pos for pos in positions}

    for marker in Marker.all():
        marker_positions = []
        for mutation in marker.mutations:
            pos = mapping.get(mutation, None)
            if pos is not None:
                marker_positions.append(pos)

        ms = MarkerScan(marker, marker_positions)
        if not ms.is_detected:
            continue
        if not relaxed and not ms.is_complete:
            continue
        markers.append(ms)
    return markers
