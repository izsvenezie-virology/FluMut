from dataclasses import dataclass, field

from flumut.flumutdb.models import (
    Annotation,
    BaseModel,
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
    cascade: dict[str, str] = field(default_factory=dict)


class DeleteValidator:
    def __init__(self, instance: BaseModel) -> None:
        self.instance = instance
        self.policy: DeletePolicy = POLICIES[type(instance)]
        self.blocking_items: dict[str, list[BaseModel]] = self._count_dependants(self.policy.blocked_by)
        self.cascade_items: dict[str, list[BaseModel]] = self._count_dependants(self.policy.cascade)

    def can_delete(self) -> bool:
        return not self.blocking_items

    def _count_dependants(self, attributes: dict[str, str]) -> dict[str, list[BaseModel]]:
        result = {}
        for name, attr in attributes.items():
            if items := getattr(self.instance, attr):
                result[name] = items
        return result


POLICIES: dict[type[BaseModel], DeletePolicy] = {
    Segment: DeletePolicy(blocked_by={'Proteins': 'proteins', 'References': 'references'}),
    Protein: DeletePolicy(blocked_by={'Mutations': 'mutations'}, cascade={'Annotations': 'annotations'}),
    Reference: DeletePolicy(cascade={'Annotations': 'annotations', 'Mappings': 'mappings'}),
    Mutation: DeletePolicy(blocked_by={'Markers': 'markers'}, cascade={'Mappings': 'mappings'}),
    Marker: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Paper: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Effect: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Subtype: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Host: DeletePolicy(blocked_by={'Evidences': 'evidences'}),
    Annotation: DeletePolicy(),
    Mapping: DeletePolicy(),
    Evidence: DeletePolicy(),
}
