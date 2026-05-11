from dataclasses import dataclass, field

from flumut.alignment import NucleotideAlignment
from flumut.core.models import Alignment
from flumut.flumutdb import Protein


@dataclass
class CDSAlignment(Alignment):
    alignment: NucleotideAlignment
    protein: Protein

    # From nucleotide alignment properties
    frameshifts: list[tuple[int, int]] = field(default_factory=list)

    def adjust_frame(self) -> None:
        pass
