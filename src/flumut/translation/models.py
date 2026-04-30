from dataclasses import dataclass, field

from flumutdb.models import Protein

from flumut.alignment.models import NucleotideAlignment
from flumut.core.models import Alignment


@dataclass
class ProteinAlignment(Alignment):
    nucleotides: NucleotideAlignment
    protein: Protein
    frameshifts: list[tuple[int, int]] = field(default_factory=list)
