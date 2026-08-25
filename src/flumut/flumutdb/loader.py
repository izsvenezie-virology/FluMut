"""In-memory image of the database.

Every table is read once, then all relations are wired by hand: the forward
foreign key on each child and the matching backref on each parent. Callers can
then walk the whole model graph in any direction without touching the database
again. The image is reused until :func:`clear` drops it.
"""

from collections.abc import Iterable
from typing import Any, TypeVar

from peewee import Model

from flumut.flumutdb.models import (
    Annotation,
    Effect,
    Evidence,
    Host,
    Mapping,
    Marker,
    MarkerMutation,
    Mutation,
    Paper,
    Protein,
    Reference,
    Segment,
    Subtype,
)

ModelT = TypeVar('ModelT', bound=Model)

_segments: list[Segment] = []
_references: list[Reference] = []
_markers: list[Marker] = []
_loaded = False


def _by_id(model: type[ModelT]) -> dict[int, ModelT]:
    return {instance.get_id(): instance for instance in model.select()}


def _relate(model: type[Model], field_name: str, children: Iterable[Model], parents: dict[int, Any]) -> None:
    field = model._meta.fields[field_name]  # type: ignore[attr-defined]
    for parent in parents.values():
        setattr(parent, field.backref, [])
    for child in children:
        parent = parents.get(child.__data__.get(field_name))  # type: ignore[attr-defined]
        if parent is None:
            continue
        child.__rel__[field_name] = parent  # type: ignore[attr-defined]
        getattr(parent, field.backref).append(child)


def load(force_reload: bool = False) -> None:
    """Read every table and link the whole model graph.

    Args:
        force_reload: Re-read from the database even if already loaded.
    """
    global _segments, _references, _markers, _loaded
    if _loaded and not force_reload:
        return

    segments = _by_id(Segment)
    proteins = _by_id(Protein)
    references = _by_id(Reference)
    annotations = _by_id(Annotation)
    mutations = _by_id(Mutation)
    mappings = _by_id(Mapping)
    markers = _by_id(Marker)
    evidences = _by_id(Evidence)
    papers = _by_id(Paper)
    effects = _by_id(Effect)
    subtypes = _by_id(Subtype)
    hosts = _by_id(Host)
    marker_mutations = list(MarkerMutation.select())

    _relate(Protein, 'segment', proteins.values(), segments)
    _relate(Reference, 'segment', references.values(), segments)
    _relate(Annotation, 'protein', annotations.values(), proteins)
    _relate(Annotation, 'reference', annotations.values(), references)
    _relate(Mutation, 'protein', mutations.values(), proteins)
    _relate(Mapping, 'mutation', mappings.values(), mutations)
    _relate(Mapping, 'reference', mappings.values(), references)
    _relate(Evidence, 'marker', evidences.values(), markers)
    _relate(Evidence, 'paper', evidences.values(), papers)
    _relate(Evidence, 'effect', evidences.values(), effects)
    _relate(Evidence, 'subtype', evidences.values(), subtypes)
    _relate(Evidence, 'host', evidences.values(), hosts)

    # Linking the through model in both directions is what makes the
    # Marker <-> Mutation many-to-many resolve without a query: the accessor
    # short-circuits as soon as the through backref is a list. The field
    # itself must never be assigned, as that would write to the database.
    _relate(MarkerMutation, 'marker', marker_mutations, markers)
    _relate(MarkerMutation, 'mutation', marker_mutations, mutations)

    _index_by_protein(proteins.values())

    _segments = list(segments.values())
    _references = [reference for segment in _segments for reference in segment.references]
    _markers = list(markers.values())
    _loaded = True


def _index_by_protein(proteins: Iterable[Protein]) -> None:
    """Group each reference's mappings and annotations by protein.

    Scanning and translation look these up per (reference, protein) pair, which
    is an intersection neither backref alone provides.
    """
    for protein in proteins:
        for mutation in protein.mutations:
            for mapping in mutation.mappings:
                mapping.reference.mappings_by_protein[protein].append(mapping)
        for annotation in protein.annotations:
            annotation.reference.annotations_by_protein[protein].append(annotation)


def load_segments(force_reload: bool = False) -> list[Segment]:
    """Return all Segments, with the whole model graph attached.

    Args:
        force_reload: Re-read from the database even if already loaded.
    """
    load(force_reload)
    return _segments


def load_references(force_reload: bool = False) -> list[Reference]:
    """Return all References, grouped by segment.

    Args:
        force_reload: Re-read from the database even if already loaded.
    """
    load(force_reload)
    return _references


def load_markers(force_reload: bool = False) -> list[Marker]:
    """Return all Markers, with mutations and evidences attached.

    Args:
        force_reload: Re-read from the database even if already loaded.
    """
    load(force_reload)
    return _markers


def clear() -> None:
    """Drop the loaded image, forcing a reload on next access."""
    global _segments, _references, _markers, _loaded
    _segments = []
    _references = []
    _markers = []
    _loaded = False
