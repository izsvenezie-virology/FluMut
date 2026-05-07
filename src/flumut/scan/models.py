from dataclasses import dataclass, field

from flumutdb import Mapping, Marker, Mutation
from flumutdb.models import MutationType


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
