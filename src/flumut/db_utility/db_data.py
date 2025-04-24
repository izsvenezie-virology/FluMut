from collections import defaultdict
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

from flumut.db_utility.db_connection import execute_query, to_dict
from flumut.db_utility.db_models import Mutation, MutationMapping, Reference


@lru_cache
def references() -> List[Reference]:
    '''
    Retrieve all the references present in the database.

    :return `List[Reference]` references: The list of references
    '''
    res = execute_query("SELECT segment_name, name, sequence FROM 'references'")
    segments = []
    for segment, name, sequence in res:
        segments.append(Reference(segment, name, sequence))
    return segments


@lru_cache
def annotations() -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    '''
    Retrieve all the annotations present in the database.

    :return `Dict[str, Dict[str, List[Tuple[int, int]]]]`: The first dictionary has the reference name as key.
    The second dictionary has the name of the protein as key and a list of start/end values for each part of the coding sequence in the protein.
    '''
    res = execute_query("SELECT reference_name, protein_name, start, end FROM 'annotations'")
    ann = defaultdict(lambda: defaultdict(list))
    for ref, prot, start, end in res:
        ann[ref][prot].append((start, end))
    return ann


@lru_cache
def mutations() -> Dict[str, List[Mutation]]:
    '''
    Retrieve all the mutations from the database, grouped by protein.

    :returns `Dict[str, List[Mutation]]`: The key is the protein name, the value is a list of the mutations in that protein.
    '''
    mutations = defaultdict(list)
    res_mut = execute_query("""SELECT segments.name, segments.number, mutations.name, mutations.type, mutations.protein_name, mutations.default_position FROM 'mutations'
								INNER JOIN proteins ON mutations.protein_name == proteins.name
                                INNER JOIN segments ON proteins.segment_name == segments.name""")
    for segment, segment_number, name, type, protein_name, default_position in res_mut:
        mutation = Mutation(segment, int(segment_number), name, type, protein_name, int(default_position))
        mutations[protein_name].append(mutation)
        res_map = execute_query(
            f"SELECT reference_name, position, ref_seq, alt_seq FROM mutation_mappings WHERE mutation_name = '{name}'")
        for reference_name, position, ref_seq, alt_seq in res_map:
            mutation.mappings[reference_name] = MutationMapping(reference_name, position, alt_seq)
    return mutations


def references_by_segment(segment: Optional[str]) -> List[Reference]:
    '''
    Returns all the references available for a specific segment.
    If no references are available, all the references are returned.

    :param `str`/`None` segment: The name of the segment. If `None` all references are returned.
    :return `List[Reference]`: List of references available for the segment.
    '''
    refs = [ref for ref in references() if ref.segment == segment]
    if not refs:
        refs = references()
    return refs


def annotations_by_reference(ref_name: str) -> Dict[str, List[Tuple[int, int]]]:
    '''
    Returns the annotations for all the proteins present for a reference.

    :param `str` ref_name: Reference name.
    :return `Dict[str, List[Tuple[int, int]]]`: The dictionary has the name of the protein as key and a list of start/end values for each part of the coding sequence as value.
    '''
    return annotations()[ref_name]


def mutations_by_protein(protein_name: str) -> List[Mutation]:
    '''
    Returns all mutations for a specific protein.

    :param `str` protein_name: Name of the protein.
    :return `List[Mutation]`: List of all mutation present in the protein.
    '''
    return mutations()[protein_name]


def markers_by_mutations(mutations: List[Mutation], relaxed: bool) -> List[Dict[str, str]]:
    '''
    Retrieve all markers matching a set of mutations.
    The match can be relaxed (a marker is returned if at least one of its mutations is found)
    or strict (a marker is returned if all of its mutations is found).

    :param `List[Mutation]` mutations: List of the mutations to match with markers.
    :param `bool` relaxed: Selects the result type: relaxed or strict.
    :return `List[Dict[str,str]]`: The list of markers matched. Keys are the name of the property, values are their value.
    '''
    muts_str = ','.join([f"'{mut.name}'" for mut in mutations])
    res = execute_query(f"""
    SELECT mm.marker_id,
        -- Get all mutations composing marker
        (SELECT group_concat(mutation_name, '; ') 
        FROM (
            SELECT mutation_name
            FROM markers_mutations 
            JOIN mutations ON mutation_name = mutations.name
            JOIN proteins ON mutations.protein_name = proteins.name
            JOIN segments ON proteins.segment_name = segments.name
            WHERE marker_id = mm.marker_id 
            ORDER BY segments.number, proteins.name, mutations.default_position)) AS 'Marker',

        -- Get all mutations composing marker found in the samples
        (SELECT group_concat(mutation_name, '; ')
        FROM (
            SELECT DISTINCT mutation_name 
            FROM markers_mutations
            JOIN mutations ON mutation_name = mutations.name
            JOIN proteins ON mutations.protein_name = proteins.name
            JOIN segments ON proteins.segment_name = segments.name
            WHERE mutation_name IN ({muts_str})
                AND marker_id = mm.marker_id
            ORDER BY segments.number, proteins.name, mutations.default_position)) AS 'Mutations in your sample',

        -- Get marker information
        me.effect_name AS 'Effect',
        me.subtype AS 'Subtype',
        (SELECT group_concat(paper_id, '; ')
        FROM (SELECT DISTINCT paper_id
            FROM markers_effects
            WHERE paper_id IN (SELECT paper_id FROM markers_effects WHERE marker_id = mm.marker_id)
            ORDER BY paper_id)) AS 'Literature',

        -- Counts mutations in the marker
        count(DISTINCT mutation_name) AS marker_found, 
        -- Counts mutations composing marker found in the sample
        (SELECT count(*)
            FROM markers_mutations
            WHERE marker_id = mm.marker_id) AS marker_tot

        FROM markers_mutations AS mm
        JOIN markers_effects AS me ON mm.marker_id = me.marker_id
        WHERE mm.mutation_name IN ({muts_str}) 
        GROUP BY mm.marker_id, me.effect_name, me.subtype
        {"HAVING marker_found = marker_tot" if not relaxed else ''};
    """, to_dict)
    return res.fetchall()


def literature_by_ids(ids: List[str]) -> List[Dict[str, str]]:
    '''
    Retrieve all literature records with ID present in a list of IDs.

    :param `List[str]` ids: The list of ids to retrieve.
    :return `List[Dict[str, str]]`: The list of literature records matched. Keys are the name of the property, values are their value.
    '''
    id_list = ','.join(map(lambda x: f"'{x}'", ids))
    res = execute_query(f"""
                        SELECT id AS 'Short name',
                               title AS 'Title',
                               authors AS 'Authors',
                               year AS 'Year',
                               journal AS 'Journal',
                               web_address AS 'Link',
                               doi AS 'DOI'
                        FROM papers
                        WHERE id IN ({id_list})
                        """, to_dict)
    return res.fetchall()
