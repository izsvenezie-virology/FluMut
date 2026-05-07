from dataclasses import dataclass, field

from flumutdb import Reference
from flumutdb.models import Protein

from flumut.alignment.models import NucleotideAlignment
from flumut.core.models import Alignment


@dataclass
class ProteinAlignment(Alignment):
    protein: Protein
    reference: Reference

    # From nucleotide alignment properties
    nucleotides: NucleotideAlignment | None = None
    frameshifts: list[tuple[int, int]] = field(default_factory=list)
