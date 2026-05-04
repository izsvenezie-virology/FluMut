from dataclasses import dataclass

from Bio.Align import Alignment as BioAlignment
from Bio.SeqRecord import SeqRecord
from flumutdb import Reference

from flumut.core.models import Alignment


@dataclass
class NucleotideAlignment(Alignment):
    """Stores the alignment of a sequence to a reference."""

    query: SeqRecord
    """Query sequence."""
    reference: Reference = None  # type: ignore[assignment]
    """Reference sequence."""

    score: float = -1_000_000

    def set_alignment(self, reference: Reference, alignment: BioAlignment) -> None:
        self.reference = reference
        self.aligned_reference = list(alignment.target)
        self.aligned_query = list(alignment.query)
        self.score = alignment.score  # type: ignore
