from __future__ import annotations

from typing import Dict, List


class Segment:
    Segments: Dict[str, Segment] = {}

    def __init__(self, name: str) -> None:
        self.name = name
        self.proteins: Dict[str, Protein] = dict()
        self.references: Dict[str, Reference] = dict()

        Segment.Segments[name] = self


class Protein:
    def __init__(self, name: str, segment: Segment) -> None:
        self.name = name
        self.mutations: Dict[str, Mutation] = list()


class Mutation:
    def __init__(self, name: str, type: str, protein: Protein) -> None:
        self.name = name
        self.type = type
        self.protein = protein
        self.mappings: List[MappedMutation] = list()


class Reference:
    References: List[Reference] = []

    def __init__(self, name: str, sequence: str, segment: Segment) -> None:
        self.name = name
        self.sequence = sequence
        self.segment = segment
        self.mapped_proteins: List[MappedProtein] = list()

        Reference.References.append(self)


class MappedProtein:
    def __init__(self, reference: Reference, protein: Protein) -> None:
        self.reference = reference
        self.protein = protein
        self.annotations: List[Annotation] = list()
        self.mutations: List[MappedMutation] = list()


class MappedMutation:
    def __init__(self, mutation: Mutation, ref_seq: str, alt_seq: str, position: int) -> None:
        self.mutation = mutation
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        self.position = position


class Annotation:
    def __init__(self, start: int, end: int) -> None:
        self.start = start
        self.end = end
