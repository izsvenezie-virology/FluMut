from dataclasses import dataclass

from flumut.core.globals import GAP_SYMBOL


@dataclass(kw_only=True)
class Alignment:
    """Stores the alignment of a sequence to a reference."""

    aligned_reference: str = ''
    aligned_query: str = ''

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
