from dataclasses import dataclass

from flumutdb import Mapping


@dataclass
class PositionScan:
    mapping: Mapping
    ammino_acid: str
