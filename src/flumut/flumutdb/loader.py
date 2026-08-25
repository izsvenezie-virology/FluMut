"""Bulk loading of the database models.

Each loader prefetches a whole branch of the model graph in a handful of
queries and keeps the result, so callers can walk the relations without
issuing one query per row. Loaded data is reused until :func:`clear` drops it.
"""

from peewee import prefetch

from flumut.flumutdb.models import (
    Annotation,
    Effect,
    Evidence,
    Host,
    Mapping,
    Marker,
    Mutation,
    Paper,
    Protein,
    Reference,
    Segment,
    Subtype,
)

_segments: list[Segment] = []
_references: list[Reference] = []
_markers: list[Marker] = []


def load_segments(force_reload: bool = False) -> list[Segment]:
    """Return all Segments, with proteins, references, mutations and mappings attached.

    Args:
        force_reload: Re-fetch from the database even if already loaded.
    """
    global _segments
    if not _segments or force_reload:
        _segments = list(
            prefetch(
                Segment.select(),
                Reference.select(),
                Protein.select(),
                Annotation.select(),
                Mutation.select(),
                Mapping.select(),
            )
        )
        ref_by_id = {ref.get_id(): ref for seg in _segments for ref in seg.references}
        for seg in _segments:
            for protein in seg.proteins:
                for mutation in protein.mutations:
                    for mapping in mutation.mappings:
                        ref_by_id[mapping.reference_id].mappings_by_protein[protein].append(mapping)  # type: ignore[attr-defined]
                for annotation in protein.annotations:
                    ref_by_id[annotation.reference_id].annotations_by_protein[protein].append(annotation)  # type: ignore[attr-defined]
    return _segments


def load_references(force_reload: bool = False) -> list[Reference]:
    """Return all References, grouped by segment.

    Args:
        force_reload: Re-fetch from the database even if already loaded.
    """
    global _references
    if not _references or force_reload:
        _references = [ref for seg in load_segments(force_reload) for ref in seg.references]
    return _references


def load_markers(force_reload: bool = False) -> list[Marker]:
    """Return all Markers, with mutations and evidences attached.

    Args:
        force_reload: Re-fetch from the database even if already loaded.
    """
    global _markers
    if not _markers or force_reload:
        _markers = list(
            prefetch(
                Marker.select(),
                Marker.mutations.through_model.select(),  # type: ignore[union-attr]
                Mutation.select(),
                Evidence.select(),
                Paper.select(),
                Effect.select(),
                Host.select(),
                Subtype.select(),
            )
        )
    return _markers


def load_all(force_reload: bool = False) -> None:
    """Load every branch of the model graph.

    Args:
        force_reload: Re-fetch from the database even if already loaded.
    """
    load_segments(force_reload)
    load_references(force_reload)
    load_markers(force_reload)


def clear() -> None:
    """Drop all loaded data, forcing a reload on next access."""
    global _segments, _references, _markers
    _segments = []
    _references = []
    _markers = []
