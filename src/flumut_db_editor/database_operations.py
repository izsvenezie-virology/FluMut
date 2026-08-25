from dataclasses import dataclass, field
from typing import Callable

from peewee import Model

from flumut.flumutdb import loader
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


class ForeignKeyViolationError(Exception):
    """Raised when deleting a record that still has protected dependents."""

    def __init__(self, item_type: str, item_name: str, violations: dict[str, list]):
        self.item_type = item_type
        self.item_name = item_name
        self.violations = violations
        count = sum(len(rows) for rows in violations.values())
        super().__init__(f"Cannot delete {item_type} '{item_name}': {count} dependent records")


@dataclass(frozen=True)
class DeletePolicy:
    """How a model may be deleted.

    Attributes:
        blocked_by: Maps a human label (e.g. ``'Proteins'``) to the backref
            attribute on the instance (e.g. ``'proteins'``). A non-empty backref
            blocks deletion and is reported as a violation.
        cascade: Callables run (in order) to remove owned children before the
            record itself is deleted.
    """

    blocked_by: dict[str, str] = field(default_factory=dict)
    cascade: tuple[Callable[[Model], None], ...] = ()


def _delete_mappings(mutation: Mutation) -> None:
    Mapping.delete().where(Mapping.mutation == mutation).execute()


def _clear_marker_mutations(marker: Marker) -> None:
    marker.mutations.clear()  # remove the M2M through rows


POLICIES: dict[type[Model], DeletePolicy] = {
    Segment: DeletePolicy(blocked_by={'Proteins': 'proteins', 'References': 'references'}),
    Protein: DeletePolicy(blocked_by={'Annotations': 'annotations', 'Mutations': 'mutations'}),
    Reference: DeletePolicy(blocked_by={'Annotations': 'annotations', 'Mappings': 'mappings'}),
    Mutation: DeletePolicy(blocked_by={'Markers': 'markers'}, cascade=(_delete_mappings,)),
    Marker: DeletePolicy(blocked_by={'Evidences': 'evidences'}, cascade=(_clear_marker_mutations,)),
    Paper: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Effect: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Subtype: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Host: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Annotation: DeletePolicy(),
    Mapping: DeletePolicy(),
    Evidence: DeletePolicy(),
}

def delete(instance: Model) -> None:
    """Delete a record, blocking if protected dependents exist.

    Owned children (e.g. a mutation's mappings) are removed automatically;
    protected dependents raise :class:`ForeignKeyViolationError`.
    """
    policy = POLICIES[type(instance)]
    violations = _blocking_dependents(instance, policy)
    if violations:
        raise ForeignKeyViolationError(type(instance).__name__, str(instance), violations)

    for remove_children in policy.cascade:
        remove_children(instance)
    instance.delete_instance()
    _clear_caches()


def remove_mutation_from_marker(marker: Marker, mutation: Mutation) -> None:
    """Remove a single mutation from a marker (M2M edit, not a row delete)."""
    marker.mutations.remove(mutation)
    _clear_caches()


def _blocking_dependents(instance: Model, policy: DeletePolicy) -> dict[str, list]:
    found = {}
    for label, backref in policy.blocked_by.items():
        rows = list(getattr(instance, backref))
        if rows:
            found[label] = rows
    return found


def _clear_caches() -> None:
    """Invalidate the loaded model graph after any write."""
    loader.clear()
