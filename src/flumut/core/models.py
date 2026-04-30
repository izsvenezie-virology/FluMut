from dataclasses import dataclass, field

from flumut.core.globals import GAP_SYMBOL

# @dataclass
# class Sample:
#     id: str
#     segments: defaultdict[str, list[ProteinAlignment]] = field(default_factory=lambda: defaultdict(list))
#     scan: list[PositionScan] = field(default_factory=list)


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
