from typing import Dict, List


class Mutation:
    def __init__(self, name: str, type: str, ref: str, alt: str, pos: int) -> None:
        self.name: str = name
        self.type: str = type
        self.ref: str = ref
        self.alt: str = alt
        self.pos: int = pos
        self.found: bool = False
        self.samples: Dict[str, List] = {}
