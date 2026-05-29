from dataclasses import dataclass, field

from flumut.core.globals import GAP_SYMBOL
from flumut.core.io.input import FastaSequence
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
class Nucleotide:
    """Stores the alignment of a sequence to a reference."""

    query: FastaSequence
    """Query sequence."""
    reference: Reference
    """Reference sequence."""
    alignment: Alignment


@dataclass
class CDS:
    nucleotides: Nucleotide
    protein: Protein
    alignment: Alignment

    # From nucleotide alignment properties
    frameshifts: list[tuple[int, int]] = field(default_factory=list)
