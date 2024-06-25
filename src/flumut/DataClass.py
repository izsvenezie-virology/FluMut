from collections import defaultdict
from typing import Any, Dict, List


class Mutation:
    def __init__(self, name: str, type: str, ref: str, alt: str, pos: int) -> None:
        self.name: str = name
        self.type: str = type
        self.ref: str = ref
        self.alt: str = alt
        self.pos: int = pos
        
        self.protein: str = self.name.split(':')[0]
        self.found: bool = False
        self.samples: Dict[str, List[str]] = {}

class Sample:
    def __init__(self, name: str) -> None:
        self.name: str = name
        self.segments: Dict[str, Dict[str, Any]] = defaultdict(dict)
        self.mutations: List[Mutation] = []
        self.markers: List[Dict[str, str]] = defaultdict(list)
