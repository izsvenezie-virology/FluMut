from dataclasses import dataclass, field

from Bio.SeqRecord import SeqRecord

from flumut.core.models import Alignment
from flumut.flumutdb import Protein, Reference


@dataclass
class NucleotideAlignment:
    """Stores the alignment of a sequence to a reference."""

    query: SeqRecord
    """Query sequence."""
    reference: Reference
    """Reference sequence."""
    alignment: Alignment


@dataclass
class CDSAlignment:
    nucleotides: NucleotideAlignment
    protein: Protein
    alignment: Alignment

    # From nucleotide alignment properties
    frameshifts: list[tuple[int, int]] = field(default_factory=list)

    def adjust_frame(self) -> None:
        pass
