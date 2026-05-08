from dataclasses import dataclass

from Bio.Align import Alignment as BioAlignment
from Bio.SeqRecord import SeqRecord

from flumut.core.models import Alignment
from flumut.flumutdb.models import Reference


@dataclass
class NucleotideAlignment(Alignment):
    """Stores the alignment of a sequence to a reference."""

    query: SeqRecord
    """Query sequence."""
    reference: Reference
    """Reference sequence."""
    alignment: BioAlignment

    def __post_init__(self) -> None:
        self.aligned_reference = list(self.alignment.target)
        self.aligned_query = list(self.alignment.query)
