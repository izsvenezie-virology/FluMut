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
    res_mut = execute_query("SELECT name, type, protein_name FROM 'mutations'")
    for name, type, protein_name in res_mut:
        mutation = Mutation(name, type, protein_name)
        mutations[protein_name].append(mutation)
        res_map = execute_query(
            f"SELECT reference_name, position, ref_seq, alt_seq FROM mutation_mappings WHERE mutation_name = '{name}'")
        for reference_name, position, ref_seq, alt_seq in res_map:
            mutation.mappings[reference_name] = MutationMapping(mutation, reference_name, position, alt_seq)
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
    WITH markers_tbl AS (SELECT marker_id,
                                group_concat(mutation_name) AS found_mutations,
                                count(mutation_name) AS found_mutations_count
                            FROM markers_mutations
                            WHERE mutation_name IN ({muts_str})
                            GROUP BY markers_mutations.marker_id)

    SELECT  markers_summary.all_mutations AS 'Marker',
            markers_tbl.found_mutations AS 'Mutations in your sample',
            markers_effects.effect_name AS 'Effect', 
            group_concat(markers_effects.paper_id, '; ') AS 'Literature', 
            markers_effects.subtype AS 'Subtype'
    FROM markers_effects
    JOIN markers_tbl ON markers_tbl.marker_id = markers_effects.marker_id
    JOIN markers_summary ON markers_summary.marker_id = markers_effects.marker_id
    WHERE markers_effects.marker_id IN (
        SELECT markers_tbl.marker_id 
        FROM markers_tbl) { 'AND markers_summary.all_mutations_count = markers_tbl.found_mutations_count' if not relaxed else '' }
    GROUP BY markers_effects.marker_id, markers_effects.effect_name, markers_effects.subtype
    """, to_dict)
    return res.fetchall()
