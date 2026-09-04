from collections import defaultdict
from enum import Enum
from functools import total_ordering

from peewee import (
    DeferredThroughModel,
    ForeignKeyField,
    IntegerField,
    ManyToManyField,
    Model,
    TextField,
)

from flumut.core.globals import DATABASE_PROXY, DB_MAJOR_VERSION
from flumut.flumutdb.exceptions import IncompatibleVersionError, MissingVersionError


class MutationType(Enum):
    SNP = 'SNP'


class BaseModel(Model):
    notes: str | None = TextField(null=True)  # type: ignore[assignment]


BaseModel._meta.database = DATABASE_PROXY  # type: ignore[attr-defined]


@total_ordering
class SortableModel(BaseModel):
    order: int = IntegerField(default=1_000_000)  # type: ignore[assignment]

    @property
    def sort_key(self) -> tuple[int, ...]:
        return (self.order,)

    def __lt__(self, other: 'SortableModel') -> bool:
        return self.sort_key < other.sort_key


class Segment(SortableModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]

    proteins: list['Protein']
    references: list['Reference']

    def __str__(self) -> str:
        return str(self.name)


class Protein(SortableModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    segment: Segment = ForeignKeyField(Segment, backref='proteins', on_delete='RESTRICT')  # type: ignore[assignment]

    annotations: list['Annotation']
    mutations: list['Mutation']

    def __str__(self) -> str:
        return f'{self.segment}/{self.name}'

    @property
    def sort_key(self) -> tuple[int, ...]:
        return self.segment.sort_key + (self.order,)


class Reference(SortableModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    segment: Segment = ForeignKeyField(Segment, backref='references', on_delete='RESTRICT')  # type: ignore[assignment]
    sequence: str = TextField(unique=True)  # type: ignore[assignment]
    source: str = TextField(unique=True)  # type: ignore[assignment]

    annotations: list['Annotation']
    mappings: list['Mapping']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mappings_by_protein: dict[Protein, list[Mapping]] = defaultdict(list)
        self.annotations_by_protein: dict[Protein, list[Annotation]] = defaultdict(list)

    def __str__(self) -> str:
        return f'{self.segment}/{self.name}'

    @property
    def sort_key(self) -> tuple[int, ...]:
        return self.segment.sort_key + (self.order,)


class Annotation(SortableModel):
    protein: Protein = ForeignKeyField(Protein, backref='annotations', on_delete='CASCADE')  # type: ignore[assignment]
    reference: Reference = ForeignKeyField(Reference, backref='annotations', on_delete='CASCADE')  # type: ignore[assignment]
    start: int = IntegerField()  # type: ignore[assignment]
    end: int = IntegerField()  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.protein} @ {self.reference}: {self.start}-{self.end}'

    @property
    def sort_key(self) -> tuple[int, ...]:
        return self.reference.sort_key + (self.protein.order, self.order)


class Mutation(SortableModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    type: str = TextField(choices=[(t.value, t.name) for t in MutationType])  # type: ignore[assignment]
    protein: Protein = ForeignKeyField(Protein, backref='mutations', on_delete='RESTRICT')  # type: ignore[assignment]

    mappings: list['Mapping']
    markers: list['Marker']

    def __str__(self) -> str:
        return str(self.name)

    @property
    def sort_key(self) -> tuple[int, ...]:
        return self.protein.sort_key + (self.order,)


class Mapping(BaseModel):
    mutation: Mutation = ForeignKeyField(Mutation, backref='mappings', on_delete='CASCADE')  # type: ignore[assignment]
    reference: Reference = ForeignKeyField(Reference, backref='mappings', on_delete='CASCADE')  # type: ignore[assignment]
    mutation_name: str | None = TextField(null=True)  # type: ignore[assignment]
    position: int = IntegerField()  # type: ignore[assignment]
    alteration: str = TextField()  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.mutation} @ {self.reference} (pos {self.position}, {self.alteration})'


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


class Target(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.name)


class Paper(BaseModel):
    short_name: str = TextField(unique=True)  # type: ignore[assignment]
    title: str = TextField()  # type: ignore[assignment]
    authors: str = TextField()  # type: ignore[assignment]
    year: int = IntegerField()  # type: ignore[assignment]
    journal: str | None = TextField(null=True)  # type: ignore[assignment]
    url: str | None = TextField(null=True)  # type: ignore[assignment]
    doi: str | None = TextField(null=True, unique=True)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.short_name)


MarkerMutationThrough = DeferredThroughModel()


class Marker(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    mutations: list[Mutation] = ManyToManyField(Mutation, backref='markers', through_model=MarkerMutationThrough)  # type: ignore[assignment]

    evidences: list['Evidence']

    def __str__(self) -> str:
        return self.name


class MarkerMutation(BaseModel):
    """Through model for the Marker <-> Mutation many-to-many relation.

    Asymmetric delete policy: removing a Marker cascades away its links,
    while a Mutation cannot be removed while still referenced by a Marker.
    """

    marker: Marker = ForeignKeyField(Marker, on_delete='CASCADE')  # type: ignore[assignment]
    mutation: Mutation = ForeignKeyField(Mutation, on_delete='RESTRICT')  # type: ignore[assignment]

    class Meta:
        database = DATABASE_PROXY


MarkerMutationThrough.set_model(MarkerMutation)


class Evidence(BaseModel):
    marker: Marker = ForeignKeyField(Marker, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    paper: Paper = ForeignKeyField(Paper, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    effect: Effect = ForeignKeyField(Effect, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    subtype: Subtype = ForeignKeyField(Subtype, backref='evidences', on_delete='RESTRICT')  # type: ignore[assignment]
    host: Host | None = ForeignKeyField(Host, backref='evidences', null=True, on_delete='RESTRICT')  # type: ignore[assignment]
    target: Target | None = ForeignKeyField(Target, backref='evidences', null=True, on_delete='RESTRICT')  # type: ignore[assignment]

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
