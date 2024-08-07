from flumut.Core.Reader.DB import DBReader
from flumut.Core.Models.DataModels import Annotation, MappedMutation, MappedProtein, Mutation, Protein, Reference, Segment


def load_db_data() -> None:
    query_result = DBReader.execute_query(
        """
        SELECT name
        FROM 'segments'
        """
    )
    for row in query_result:
        name = row[0]
        segment = Segment(name)
        load_proteins(segment)
        load_references(segment)


def load_proteins(segment: Segment) -> None:
    query_result = DBReader.execute_query(
        f"""
        SELECT name
        FROM 'proteins'
        WHERE segment_name = '{segment.name}'
        """
    )
    for row in query_result:
        name = row[0]
        protein = Protein(row[0], segment)
        segment.proteins[name] = protein
        load_mutations(protein)


def load_mutations(protein: Protein) -> None:
    query_result = DBReader.execute_query(
        f"""
        SELECT name, type
        FROM 'mutations'
        WHERE protein_name = '{protein.name}'
        """
    )
    for name, mutation_type in query_result:
        protein.mutations[name] = Mutation(name, mutation_type, protein)


def load_references(segment: Segment) -> None:
    query_result = DBReader.execute_query(
        f"""
        SELECT name, sequence
        FROM 'references'
        WHERE segment_name = '{segment.name}'
        """
    )
    for name, sequence in query_result:
        reference = Reference(name, sequence, segment)
        segment.references[name] = reference
        load_mapped_proteins(reference)


def load_mapped_proteins(reference: Reference) -> None:
    query_result = DBReader.execute_query(
        f"""
        SELECT DISTINCT protein_name
        FROM 'annotations'
        WHERE reference_name = '{reference.name}'
        """
    )
    for row in query_result:
        protein = reference.segment.proteins[row[0]]
        mapped_protein = MappedProtein(reference, protein)
        reference.mapped_proteins.append(mapped_protein)
        load_annotations(mapped_protein)
        load_mapped_mutations(mapped_protein)


def load_annotations(mapped_protein: MappedProtein) -> None:
    query_result = DBReader.execute_query(
        f"""
        SELECT start, end
        FROM 'annotations'
        WHERE protein_name = '{mapped_protein.protein.name}'
            AND reference_name = '{mapped_protein.reference.name}'
        """
    )
    for start, end in query_result:
        mapped_protein.annotations.append(Annotation(start, end))


def load_mapped_mutations(mapped_protein: MappedProtein) -> None:
    query_result = DBReader.execute_query(
        f"""
        SELECT mutation_name, ref_seq, alt_seq, position
        FROM mutation_mappings
        LEFT JOIN mutations ON mutations.name = mutation_mappings.mutation_name
        WHERE mutations.protein_name = '{mapped_protein.protein.name}'
            AND refrence_name = '{mapped_protein.reference.name}'
        """
    )
    for mutation_name, ref_seq, alt_seq, position in query_result:
        mutation = mapped_protein.protein.mutations[mutation_name]
        mapped_protein.mutations.append(MappedMutation(mapped_protein, mutation, ref_seq, alt_seq, position))
