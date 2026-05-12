from dataclasses import dataclass, field

from Bio.SeqRecord import SeqRecord

from flumut.core.globals import GAP_SYMBOL
from flumut.flumutdb import Protein, Reference


@dataclass
class Alignment:
    """Stores the alignment of a sequence to a reference."""

    reference: list[str] = field(default_factory=list)
    query: list[str] = field(default_factory=list)

    def get_positions(self) -> list[int]:
        positions = []
        last_position = 0
        for r in self.reference:
            if r == GAP_SYMBOL:
                positions.append(None)
            else:
                last_position += 1
                positions.append(last_position)
        return positions


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
