from dataclasses import dataclass, field

from flumut.flumutdb import Mapping, Marker, Mutation, MutationType, Paper
from flumut.nucleotides.models import Alignment, CDSAlignment, Protein, Reference


@dataclass
class PositionScan:
    mapping: Mapping
    ammino_acid: str
    is_detected: bool = field(init=False)

    def __post_init__(self) -> None:
        match self.mapping.mutation.type:
            case MutationType.SNP.value:
                self.is_detected = self.mapping.alteration in self.ammino_acid
            case _:
                raise NotImplementedError(f'Mutation type {self.mapping.mutation.type} not supported')

    @property
    def mutation(self) -> Mutation:
        return self.mapping.mutation


@dataclass
class MarkerScan:
    marker: Marker
    positions: list[PositionScan]

    detected_mutations: list[PositionScan] = field(default_factory=list, init=False)
    is_detected: bool = field(default=False, init=False)
    is_complete: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        self.detected_mutations = [position for position in self.positions if position.is_detected]
        self.is_detected = len(self.detected_mutations) > 0
        self.is_complete = len(self.detected_mutations) == len(self.marker.mutations)


@dataclass
class ProteinAlignment:
    protein: Protein
    reference: Reference
    alignment: Alignment

    cds: 'CDSAlignment | None' = None


@dataclass
class Sample:
    id: str
    alignments: list[ProteinAlignment] = field(default_factory=list)

    positions: list[PositionScan] = field(default_factory=list, init=False)
    marker_scans: list[MarkerScan] = field(default_factory=list, init=False)


@dataclass
class Analysis:
    samples: dict[str, Sample] = field(default_factory=dict)
    mutations: set[Mutation] = field(default_factory=set)
    markers: set[Marker] = field(default_factory=set)
    literature: set[Paper] = field(default_factory=set)
