from enum import Enum
from typing import List

from peewee import (
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


class Segment(BaseModel):
    name: str = TextField()  # type: ignore[assignment]
    proteins: list['Protein']
    references: list['Reference']

    def __str__(self) -> str:
        return str(self.name)

    _cache: list['Segment'] = []

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
        return Segment._cache

    @staticmethod
    def clear_cache():
        """Clear the Segment cache, forcing a reload on next access."""
        Segment._cache = []


class Protein(BaseModel):
    name: str = TextField()  # type: ignore[assignment]
    segment: Segment = ForeignKeyField(Segment, backref='proteins')  # type: ignore[assignment]
    annotations: list['Annotation']
    mutations: list['Mutation']

    def __str__(self) -> str:
        return f'{self.segment}/{self.name}'


class Reference(BaseModel):
    name: str = TextField()  # type: ignore[assignment]
    segment: Segment = ForeignKeyField(Segment, backref='references')  # type: ignore[assignment]
    sequence: str = TextField()  # type: ignore[assignment]
    source: str = TextField()  # type: ignore[assignment]
    annotations: list['Annotation']
    mappings: list['Mapping']

    def __str__(self) -> str:
        return f'{self.segment}/{self.name}'

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


class Annotation(BaseModel):
    protein: Protein = ForeignKeyField(Protein, backref='annotations')  # type: ignore[assignment]
    reference: Reference = ForeignKeyField(Reference, backref='annotations')  # type: ignore[assignment]
    start: int = IntegerField()  # type: ignore[assignment]
    end: int = IntegerField()  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.protein} @ {self.reference}: {self.start}-{self.end}'


class Mutation(BaseModel):
    name: str = TextField(unique=True)  # type: ignore[assignment]
    type: str = TextField(choices=[(t.value, t.name) for t in MutationType])  # type: ignore[assignment]
    protein: Protein = ForeignKeyField(Protein, backref='mutations')  # type: ignore[assignment]
    default_position: int | None = IntegerField(null=True)  # type: ignore[assignment]
    mappings: list['Mapping']
    markers: list['Marker']

    def __str__(self) -> str:
        return str(self.name)


class Mapping(BaseModel):
    mutation: Mutation = ForeignKeyField(Mutation, backref='mappings')  # type: ignore[assignment]
    reference: Reference = ForeignKeyField(Reference, backref='mappings')  # type: ignore[assignment]
    mutation_name: str | None = TextField(null=True)  # type: ignore[assignment]
    position: int = IntegerField()  # type: ignore[assignment]
    alteration: str = TextField()  # type: ignore[assignment]

    def __str__(self) -> str:
        return f'{self.mutation} @ {self.reference} (pos {self.position}, {self.alteration})'


class Effect(BaseModel):
    name: str = TextField()  # type: ignore[assignment]
    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.name)


class Subtype(BaseModel):
    name: str = TextField()  # type: ignore[assignment]
    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.name)


class Host(BaseModel):
    name: str = TextField()  # type: ignore[assignment]
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
    doi: str | None = TextField(null=True)  # type: ignore[assignment]
    evidences: list['Evidence']

    def __str__(self) -> str:
        return str(self.short_name)


class Marker(BaseModel):
    name: str | None = TextField(unique=True, null=True)  # type: ignore[assignment]
    mutations: list[Mutation] = ManyToManyField(Mutation, backref='markers')  # type: ignore[assignment]
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
                )
            )
        return Marker._cache

    @staticmethod
    def clear_cache():
        """Clear the Marker cache, forcing a reload on next access."""
        Marker._cache = []


class Evidence(BaseModel):
    marker: Marker = ForeignKeyField(Marker, backref='evidences')  # type: ignore[assignment]
    paper: Paper = ForeignKeyField(Paper, backref='evidences')  # type: ignore[assignment]
    effect: Effect = ForeignKeyField(Effect, backref='evidences')  # type: ignore[assignment]
    subtype: Subtype = ForeignKeyField(Subtype, backref='evidences')  # type: ignore[assignment]
    host: Host | None = ForeignKeyField(Host, backref='evidences', null=True)  # type: ignore[assignment]

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
