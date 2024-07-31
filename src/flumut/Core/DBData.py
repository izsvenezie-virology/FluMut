from __future__ import annotations
from typing import Dict, List

from flumut.Core.DBReader import DBReader


def load_segments() -> Dict[str, Segment]:
    segments = {}
    query_result = DBReader.execute_query(
        """
        SELECT name
        FROM 'segments'
        """
    )
    for row in query_result:
        name = row[0]
        segments[name] = Segment(name)
    return segments


class Segment:
    def __init__(self, name: str) -> None:
        self.name = name
        self.proteins: Dict[str, Protein] = self.load_proteins()
        self.references: Dict[str, Reference] = self.load_references()

    def load_proteins(self) -> Dict[str, Protein]:
        proteins = {}
        query_result = DBReader.execute_query(
            f"""
            SELECT name
            FROM 'proteins'
            WHERE segment_name = '{self.name}'
            """
        )
        for row in query_result:
            name = row[0]
            proteins[name] = Protein(name, self)
        return proteins

    def load_references(self) -> Dict[str, Reference]:
        references = {}
        query_result = DBReader.execute_query(
            f"""
            SELECT name, sequence
            FROM 'references'
            WHERE segment_name = '{self.name}'
            """
        )
        for name, sequence in query_result:
            references[name] = Reference(name, sequence, self)
        return references


class Protein:
    def __init__(self, name: str, segment: Segment) -> None:
        self.name = name
        self.segment = segment
        self.mutations: Dict[str, Mutation] = self.load_mutations()

    def load_mutations(self) -> Dict[str, Mutation]:
        mutations = {}
        query_result = DBReader.execute_query(
            f"""
            SELECT name, type
            FROM 'mutations'
            WHERE protein_name = '{self.name}'
            """
        )
        query_result.rowcount
        for name, mutation_type in query_result:
            mutations[name] = Mutation(name, mutation_type, self)
        return mutations


class Mutation:
    def __init__(self, name: str, type: str, protein: Protein) -> None:
        self.name = name
        self.type = type
        self.protein = protein


class Reference:
    def __init__(self, name: str, sequence: str, segment: Segment) -> None:
        self.name = name
        self.sequence = sequence
        self.segment = segment
        self.mapped_proteins: List[MappedProtein] = self.load_mapped_proteins()

    def load_mapped_proteins(self) -> List[MappedProtein]:
        mapped_proteins = []
        query_result = DBReader.execute_query(
            f"""
            SELECT DISTINCT protein_name
            FROM 'annotations'
            WHERE reference_name = '{self.name}'
            """
        )
        for row in query_result:
            protein = self.segment.proteins[row[0]]
            mapped_proteins.append(MappedProtein(self, protein))
        return mapped_proteins


class MappedProtein:
    def __init__(self, reference: Reference, protein: Protein) -> None:
        self.reference = reference
        self.protein = protein
        self.annotations: List[Annotation] = self.load_annotations()
        self.mutations: List[MappedMutation] = list()

    def load_annotations(self) -> List[Annotation]:
        annotations = []
        query_result = DBReader.execute_query(
            f"""
            SELECT start, end
            FROM 'annotations'
            WHERE protein_name = '{self.protein.name}'
              AND reference_name = '{self.reference.name}'
            """
        )
        for start, end in query_result:
            annotations.append(Annotation(start, end))
        return annotations

    def load_mapped_mutations(self) -> List[MappedMutation]:
        mapped_mutations = []
        query_result = DBReader.execute_query(
            f"""
            SELECT mutation_name, ref_seq, alt_seq, position
            FROM mutation_mappings
            LEFT JOIN mutations ON mutations.name = mutation_mappings.mutation_name
            WHERE mutations.protein_name = '{self.protein.name}'
              AND refrence_name = '{self.reference.name}'
            """
        )
        for mutation_name, ref_seq, alt_seq, position in query_result:
            mutation = self.protein.mutations[mutation_name]
            mapped_mutations.append(MappedMutation(self, mutation, ref_seq, alt_seq, position))
        return mapped_mutations


class MappedMutation:
    def __init__(self, mapped_protein: MappedProtein, mutation: Mutation, ref_seq: str, alt_seq: str, position: int) -> None:
        self.mapped_protein = mapped_protein
        self.mutation = mutation
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        self.position = position


class Annotation:
    def __init__(self, start: int, end: int) -> None:
        self.start = start
        self.end = end
