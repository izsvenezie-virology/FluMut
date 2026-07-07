from collections import defaultdict
from enum import Enum
from functools import total_ordering
from typing import List

from peewee import (
    DeferredThroughModel,
    ForeignKeyField,
    IntegerField,
    ManyToManyField,
    Model,
    TextField,
    prefetch,
)

from flumut.core.globals import DATABASE_PROXY, DB_MAJOR_VERSION
from flumut.flumutdb.exceptions import IncompatibleVersionError, MissingVersionError


class MutationType(Enum):
    SNP = 'SNP'


class BaseModel(Model):
    notes: str | None = TextField(null=True)  # type: ignore[assignment]


BaseModel._meta.database = DATABASE_PROXY  # type: ignore[attr-defined]


@total_ordering
class Segment(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    order: int = IntegerField(null=True)  # type: ignore[assignment]

    proteins: list['Protein']
    references: list['Reference']

    def __str__(self) -> str:
        return str(self.name)

    _cache: list['Segment'] = []

    def __lt__(self, other: 'Segment') -> bool:
        return self.order < other.order

    @staticmethod
    def all(force_reload: bool = False) -> list['Segment']:
        """Return all Segment instances, cached after first call.

        Args:
            force_reload: Re-fetch from the database even if already cached.
        """
        if not Segment._cache or force_reload:
            Segment._cache = list(
                prefetch(
                    Segment.select(),
                    Reference.select(),
                    Protein.select(),
                    Annotation.select(),
                    Mutation.select(),
                    Mapping.select(),
                )
            )
            ref_by_id = {ref.get_id(): ref for seg in Segment._cache for ref in seg.references}
            for seg in Segment._cache:
                for protein in seg.proteins:
                    for mutation in protein.mutations:
                        for mapping in mutation.mappings:
                            ref_by_id[mapping.reference_id].mappings_by_protein[protein].append(mapping)  # type: ignore[attr-defined]
                    for annotation in protein.annotations:
                        ref_by_id[annotation.reference_id].annotations_by_protein[protein].append(annotation)  # type: ignore[attr-defined]
        return Segment._cache

    @staticmethod
    def clear_cache():
        """Clear the Segment cache, forcing a reload on next access."""
        Segment._cache = []


@total_ordering
class Protein(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    segment: Segment = ForeignKeyField(Segment, backref='proteins', on_delete='RESTRICT')  # type: ignore[assignment]
    order: int = IntegerField(null=True)  # type: ignore[assignment]

    annotations: list['Annotation']
    mutations: list['Mutation']

    def __str__(self) -> str:
        return f'{self.segment}/{self.name}'

    def __lt__(self, other: 'Protein') -> bool:
        return self.segment <= other.segment and self.order < other.order


@total_ordering
class Reference(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    segment: Segment = ForeignKeyField(Segment, backref='references', on_delete='RESTRICT')  # type: ignore[assignment]
    sequence: str = TextField(unique=True)  # type: ignore[assignment]
    source: str = TextField(unique=True)  # type: ignore[assignment]
    order: int = IntegerField(null=True)  # type: ignore[assignment]

    annotations: list['Annotation']
    mappings: list['Mapping']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mappings_by_protein: dict['Protein', list['Mapping']] = defaultdict(list)
        self.annotations_by_protein: dict['Protein', list['Annotation']] = defaultdict(list)

    def __str__(self) -> str:
        return f'{self.segment}/{self.name}'

    def __lt__(self, other: 'Reference') -> bool:
        return self.segment <= other.segment and self.order < other.order

    _cache: list['Reference'] = []

    @staticmethod
    def all(force_reload: bool = False) -> list['Reference']:
        """Return all Reference instances, cached after first call.

        Args:
            force_reload: Re-fetch from the database even if already cached.
        """
        if not Reference._cache or force_reload:
            Reference._cache = [ref for seg in Segment.all() for ref in seg.references]
        return Reference._cache

    @staticmethod
    def clear_cache():
        """Clear the Reference cache, forcing a reload on next access."""
        Reference._cache = []


@total_ordering
class Annotation(BaseModel):
    protein: Protein = ForeignKeyField(Protein, backref='annotations', on_delete='CASCADE')  # type: ignore[assignment]
    reference: Reference = ForeignKeyField(Reference, backref='annotations', on_delete='CASCADE')  # type: ignore[assignment]
    order: int = IntegerField(null=True)  # type: ignore[assignment]
    start: int = IntegerField()  # type: ignore[assignment]
    end: int = IntegerField()  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.protein} @ {self.reference}: {self.start}-{self.end}'

    def __lt__(self, other: 'Annotation') -> bool:
        return self.protein <= other.protein and self.order < other.order


@total_ordering
class Mutation(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    type: str = TextField(choices=[(t.value, t.name) for t in MutationType])  # type: ignore[assignment]
    protein: Protein = ForeignKeyField(Protein, backref='mutations', on_delete='RESTRICT')  # type: ignore[assignment]
    order: int = IntegerField(null=True)  # type: ignore[assignment]

    mappings: list['Mapping']
    markers: list['Marker']

    def __str__(self) -> str:
        return str(self.name)

    def __lt__(self, other: 'Mutation') -> bool:
        return self.protein <= other.protein and self.order < other.order


class Mapping(BaseModel):
    mutation: Mutation = ForeignKeyField(Mutation, backref='mappings', on_delete='CASCADE')  # type: ignore[assignment]
    reference: Reference = ForeignKeyField(Reference, backref='mappings', on_delete='CASCADE')  # type: ignore[assignment]
    mutation_name: str | None = TextField(null=True)  # type: ignore[assignment]
    position: int = IntegerField()  # type: ignore[assignment]
    alteration: str = TextField()  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.mutation} @ {self.reference} (pos {self.position}, {self.alteration})'

    def is_valid(self) -> bool:
        return self.position > 0


class Effect(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.name)


class Subtype(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.name)


class Host(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.name)


class Paper(BaseModel):
    short_name: str = TextField(unique=True)  # type: ignore[assignment]
    title: str = TextField()  # type: ignore[assignment]
    authors: str = TextField()  # type: ignore[assignment]
    year: int | None = IntegerField(null=True)  # type: ignore[assignment]
    journal: str | None = TextField(null=True)  # type: ignore[assignment]
    url: str | None = TextField(null=True)  # type: ignore[assignment]
    doi: str | None = TextField(null=True, unique=True)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.short_name)


MarkerMutationThrough = DeferredThroughModel()


class Marker(BaseModel):
    name: str | None = TextField(unique=True, null=True)  # type: ignore[assignment]
    mutations: list[Mutation] = ManyToManyField(Mutation, backref='markers', through_model=MarkerMutationThrough)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        mutations = ', '.join(str(m) for m in self.mutations)
        return f'Marker({mutations})'

    _cache: List['Marker'] = []

    @staticmethod
    def all(force_reload: bool = False) -> list['Marker']:
        """Returns a list of all Marker instances, cached after first call.

        Args:
            force_reload: Re-fetch from the database even if already cached.
        """
        if not Marker._cache or force_reload:
            Marker._cache = list(
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
        return Marker._cache

    @staticmethod
    def clear_cache():
        """Clear the Marker cache, forcing a reload on next access."""
        Marker._cache = []


class MarkerMutation(Model):
    """Through model for the Marker <-> Mutation many-to-many relation.

    Asymmetric delete policy: removing a Marker cascades away its links,
    while a Mutation cannot be removed while still referenced by a Marker.
    """

    marker = ForeignKeyField(Marker, on_delete='CASCADE')
    mutation = ForeignKeyField(Mutation, on_delete='RESTRICT')

    class Meta:
        database = DATABASE_PROXY


MarkerMutationThrough.set_model(MarkerMutation)


class Evidence(BaseModel):
    marker: Marker = ForeignKeyField(Marker, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    paper: Paper = ForeignKeyField(Paper, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    effect: Effect = ForeignKeyField(Effect, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    subtype: Subtype = ForeignKeyField(Subtype, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    host: Host | None = ForeignKeyField(Host, backref='evidences', null=True, on_delete='RESTRICT')  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.marker}: {self.effect} in {self.subtype} ({self.paper})'


class DbVersion(BaseModel):
    major: int = IntegerField()  # type: ignore[assignment]
    minor: int = IntegerField()  # type: ignore[assignment]
    patch: int = IntegerField()  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.major}.{self.minor}.{self.patch}'

    @staticmethod
    def is_compatible() -> bool:
        version: DbVersion = DbVersion.get_or_none()
        if version is None:
            raise MissingVersionError()
        if version.major != DB_MAJOR_VERSION:
            raise IncompatibleVersionError(version, DB_MAJOR_VERSION)
        return True
