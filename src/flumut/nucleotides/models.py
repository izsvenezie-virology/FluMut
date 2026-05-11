from dataclasses import dataclass, field

from Bio.Align import Alignment as BioAlignment
from Bio.SeqRecord import SeqRecord

from flumut.core.models import Alignment
from flumut.flumutdb import Protein, Reference


@dataclass
class NucleotideAlignment(Alignment):
    """Stores the alignment of a sequence to a reference."""

    query: SeqRecord
    """Query sequence."""
    reference: Reference
    """Reference sequence."""
    alignment: BioAlignment

    def __post_init__(self) -> None:
        self.aligned_reference = list(self.alignment[0])  # type: ignore
        self.aligned_query = list(self.alignment[1])  # type: ignore


@dataclass
class CDSAlignment(Alignment):
    alignment: NucleotideAlignment
    protein: Protein

    # From nucleotide alignment properties
    frameshifts: list[tuple[int, int]] = field(default_factory=list)

    def adjust_frame(self) -> None:
        pass
