from dataclasses import dataclass, field

from flumut.core.logger import LOGGER
from flumut.core.nucleotides.models import CDS, Alignment, Protein, Reference
from flumut.flumutdb import Mapping, Marker, Mutation, MutationType, Paper


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

    cds: CDS | None = None


@dataclass
class Sample:
    id: str
    alignments: list[ProteinAlignment] = field(default_factory=list)

    positions: list[PositionScan] = field(default_factory=list, init=False)
    marker_scans: list[MarkerScan] = field(default_factory=list, init=False)
    checks: list['Check'] = field(default_factory=list, init=False)


@dataclass
class Analysis:
    samples: dict[str, Sample] = field(default_factory=dict)
    mutations: set[Mutation] = field(default_factory=set)
    markers: set[Marker] = field(default_factory=set)
    literature: set[Paper] = field(default_factory=set)


class Check:
    def __init__(self, message: str):
        self.message = message
        LOGGER.warning(self.message)


class TruncationCheck(Check):
    def __init__(self, sample_id: str, protein_name: str, position: int):
        message = f'Premature stop codon detected in sample "{sample_id}" for protein "{protein_name}" at {position}'
        super().__init__(message)


class EnlongationCheck(Check):
    def __init__(self, sample_id: str, protein_name: str):
        message = f'Enlongation detected in sample "{sample_id}" for protein "{protein_name}"'
        super().__init__(message)


class FrameshiftCheck(Check):
    def __init__(self, sample_id: str, protein_name: str, start: int, end: int | None):
        message = f'Frameshift detected in sample "{sample_id}" for protein "{protein_name}" from position {start} to {end or "end"}'
        super().__init__(message)


class DuplicationCheck(Check):
    def __init__(self, sample_id: str, protein_name: str):
        message = f'Duplicate protein "{protein_name}" detected in sample "{sample_id}"'
        super().__init__(message)
