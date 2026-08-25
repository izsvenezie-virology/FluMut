from collections import defaultdict
from collections.abc import Callable
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
    MutationType,
    Paper,
    Protein,
    Reference,
    Segment,
    Subtype,
)


@dataclass(frozen=True)
class DeletePolicy:
    blocked_by: dict[str, str] = field(default_factory=dict)
    cascade: dict[str, str] = field(default_factory=dict)


DELETION_POLICIES: dict[type[BaseModel], DeletePolicy] = {
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


class DeleteValidator:
    def __init__(self, instance: BaseModel) -> None:
        self.instance = instance
        self.policy: DeletePolicy = DELETION_POLICIES[type(instance)]
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


def validate_sequence(instance: Reference) -> list[str]:
    messages = []
    if unknown_nucleotides := set(instance.sequence).difference({'A', 'T', 'C', 'G'}):
        messages.append(f'sequence contains unknown nucleotides "{", ".join(unknown_nucleotides)}"')
    return messages


def validate_annotation_positions(instance: Annotation) -> list[str]:
    messages = []
    max_value = len(instance.reference.sequence)
    if not 1 <= instance.start <= max_value:
        messages.append(f'start position must be between 1 and {max_value}')
    if not 1 <= instance.end <= max_value:
        messages.append(f'end position must be between 1 and {max_value}')
    if not instance.start < instance.end:
        messages.append('start position must be lower than end position')
    return messages


def validate_mutation_mappings(instance: Mutation) -> list[str]:
    messages = []
    if instance.protein is None:
        return
    valid_references = instance.protein.segment.references
    for mapping in instance.mappings:
        if mapping.reference not in valid_references:
            messages.append(f'mapping has reference {mapping.reference} that is invalid.')
    return messages


def validate_mapping_position(instance: Mapping) -> list[str]:
    messages = []
    max_value = sum([ann.end - ann.start for ann in instance.reference.annotations if ann.protein == instance.mutation.protein]) / 3
    if not 0 <= instance.position <= max_value:
        messages.append(f'position must be between 1 and {max_value}')
    return messages


def validate_mapping_name(instance: Mapping):
    messages = []
    if (dup := instance.get_or_none(instance.mutation_name)) and dup != instance:
        messages.append(f'mutation name "{instance.mutation_name}" is already present into another mapping')
    if Mutation.get_or_none(instance.mutation_name):
        messages.append(f'mutation name "{instance.mutation_name} is already present into a mutation"')
    return messages


def validate_mapping_alteration(instance: Mapping):
    messages = []
    match instance.mutation.type:
        case MutationType.SNP:
            if len(instance.alteration) != 1:
                messages.append('alteration in SNP mutations must be composed by only one letter')
            if instance.alteration not in 'ACDEFGHIKLMNPQRSTVWY':
                messages.append(f'alteration "{instance.alteration}" is not a valid amino acid code')
        case _:
            pass
    return messages


def validate_paper_title(instance: Paper):
    messages = []
    titles = [paper.title.lower() for paper in Paper.select() if paper != instance]
    if instance.title.lower() in titles:
        messages.append('title already exists')
    return messages


@dataclass(frozen=True)
class ValidatePolicy:
    not_null: dict[str, str] = field(default_factory=dict)
    unique: dict[str, str] = field(default_factory=dict)
    not_null_unique: dict[str, str] = field(default_factory=dict)
    extra_validators: dict[str, Callable] = field(default_factory=dict)


VALIDATION_POLICIES: dict[type[BaseModel], ValidatePolicy] = {
    Segment: ValidatePolicy(
        not_null_unique={'Name': 'name'},
    ),
    Protein: ValidatePolicy(
        not_null_unique={'Name': 'name'},
        not_null={'Segment': 'segment'},
    ),
    Reference: ValidatePolicy(
        not_null_unique={'Name': 'name', 'Sequence': 'sequence', 'Source': 'source'},
        not_null={'Segment': 'segment'},
        extra_validators={'Sequence': validate_sequence},
    ),
    Annotation: ValidatePolicy(
        not_null={'Protein': 'protein', 'Reference': 'reference'},
        extra_validators={'Start and end': validate_annotation_positions},
    ),
    Mutation: ValidatePolicy(
        not_null_unique={'Name': 'name'},
        not_null={'Type': 'type', 'Protein': 'protein'},
        extra_validators={'Consistency': validate_mutation_mappings},
    ),
    Mapping: ValidatePolicy(
        not_null={'Mutation': 'mutation', 'Reference': 'reference'},
        extra_validators={'Position': validate_mapping_position, 'Name': validate_mapping_name, 'Alteration': validate_mapping_alteration},
    ),
    Effect: ValidatePolicy(
        not_null_unique={'Name': 'name'},
    ),
    Subtype: ValidatePolicy(
        not_null_unique={'Name': 'name'},
    ),
    Host: ValidatePolicy(
        not_null_unique={'Name': 'name'},
    ),
    Paper: ValidatePolicy(
        not_null_unique={'Short name': 'short_name'},
        not_null={'Title': 'title', 'Authors': 'authors'},
        unique={'DOI': 'doi'},
        extra_validators={'Title': validate_paper_title},
    ),
    Marker: ValidatePolicy(
        unique={'Name': 'name'},
    ),
    Evidence: ValidatePolicy(
        not_null={'Marker': 'marker', 'Paper': 'paper', 'Effect': 'effect', 'Subtype': 'subtype'},
    ),
}


class DataValidator:
    def __init__(self, instance: BaseModel) -> None:
        self.instance = instance
        self.policy: ValidatePolicy = VALIDATION_POLICIES[type(instance)]
        self.errors: dict[str, list[str]] = defaultdict(list)

    def validate(self) -> None:
        self.errors.clear()
        self._validate_not_nulls()
        self._validate_unique()
        self._validate_extra()

    def is_valid(self) -> bool:
        return not self.errors

    def _validate_not_nulls(self) -> None:
        for name, attr in {**self.policy.not_null, **self.policy.not_null_unique}.items():
            if not getattr(self.instance, attr):
                self.errors[name].append('cannot be null.')

    def _validate_unique(self) -> None:
        for name, attr in {**self.policy.unique, **self.policy.not_null_unique}.items():
            value = getattr(self.instance, attr)
            found = type(self.instance).get_or_none(getattr(type(self.instance), attr) == value)
            if found is not None and found != self.instance:
                self.errors[name].append('must be unique.')

    def _validate_extra(self) -> None:
        for name, func in self.policy.extra_validators.items():
            errors = func(self.instance)
            if errors:
                self.errors[name] += errors
