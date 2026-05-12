from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from flumut.core.globals import GAP_SYMBOL
from flumut.flumutdb import Protein, Reference

if TYPE_CHECKING:
    from flumut.nucleotides.models import CDSAlignment


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
class ProteinAlignment:
    protein: Protein
    reference: Reference
    alignment: Alignment

    cds: 'CDSAlignment | None' = None
