from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from flumut.core.globals import GAP_SYMBOL
from flumut.flumutdb import Protein, Reference

if TYPE_CHECKING:
    from flumut.nucleotides.models import CDSAlignment


@dataclass(kw_only=True)
class Alignment:
    """Stores the alignment of a sequence to a reference."""

    aligned_reference: list[str] = field(default_factory=list)
    aligned_query: list[str] = field(default_factory=list)

    def get_positions(self) -> list[int]:
        positions = []
        last_position = 0
        for r in self.aligned_reference:
            if r == GAP_SYMBOL:
                positions.append(None)
            else:
                last_position += 1
                positions.append(last_position)
        return positions


@dataclass
class ProteinAlignment(Alignment):
    protein: Protein
    reference: Reference

    cds: 'CDSAlignment | None' = None
