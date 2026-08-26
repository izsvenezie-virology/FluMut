"""In-memory image of the database.

Every table in :data:`_MODELS` is read once, then all relations are wired: the
forward foreign key on each child and the matching backref on each parent.
Callers ask for any of those models with :func:`load_all` and can then walk the
whole model graph in any direction without touching the database again. The
image is reused until :func:`clear` drops it.
"""

from collections.abc import Iterable
from typing import Any, TypeVar

from peewee import ForeignKeyField, Model

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

_MODELS: tuple[type[Model], ...] = (
    Segment,
    Protein,
    Reference,
    Annotation,
    Mutation,
    Mapping,
    Marker,
    MarkerMutation,
    Evidence,
    Paper,
    Effect,
    Subtype,
    Host,
)

_instances: dict[type[Model], list[Any]] = {}


def _by_id(model: type[ModelT]) -> dict[int, ModelT]:
    return {instance.get_id(): instance for instance in model.select()}


def _relate(field: ForeignKeyField, children: Iterable[Model], parents: dict[int, Any]) -> None:
    """Attach the parent to each child and the children to each parent.

    The foreign key itself must never be assigned, as that would write to the
    database, so the resolved parent goes straight into the relation cache.
    """
    for parent in parents.values():
        setattr(parent, field.backref, [])
    for child in children:
        parent = parents.get(child.__data__.get(field.name))  # type: ignore[attr-defined]
        if parent is None:
            continue
        child.__rel__[field.name] = parent  # type: ignore[attr-defined]
        getattr(parent, field.backref).append(child)


def _index_by_protein(proteins: Iterable[Protein]) -> None:
    """Group each reference's mappings and annotations by protein."""
    for protein in proteins:
        for mutation in protein.mutations:
            for mapping in mutation.mappings:
                mapping.reference.mappings_by_protein[protein].append(mapping)
        for annotation in protein.annotations:
            annotation.reference.annotations_by_protein[protein].append(annotation)


def get(model: type[ModelT]) -> list[ModelT]:
    """Return every instance of ``model``, with the whole model graph attached.

    Args:
        model: Any model listed in :data:`_MODELS`.

    Raises:
        KeyError: If ``model`` is not part of the image.
    """
    if model not in _instances:
        raise KeyError(f'{model.__name__} is not part of the in-memory image.')
    return _instances[model]


def load() -> None:
    """Read every table and link the whole model graph."""
    by_id: dict[type[Model], dict[int, Any]] = {model: _by_id(model) for model in _MODELS}

    for model, instances in by_id.items():
        for field in model._meta.fields.values():  # type: ignore[attr-defined]
            if isinstance(field, ForeignKeyField) and field.rel_model in by_id:
                _relate(field, instances.values(), by_id[field.rel_model])

    _index_by_protein(by_id[Protein].values())

    _instances.clear()
    _instances.update({model: list(instances.values()) for model, instances in by_id.items()})
