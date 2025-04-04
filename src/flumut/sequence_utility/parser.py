from typing import List

from flumut.db_utility.db_data import mutations_by_protein
from flumut.db_utility.db_models import Mutation
from flumut.sequence_utility.models import Sample


def seek_mutations(sample: Sample) -> List[Mutation]:
    '''
    Search for mutations present in the protein sequence.
    '''
    mutations: List[Mutation] = []
    for sequence in sample.sequences:
        for protein in sequence.alignment.proteins:
            for mutation in mutations_by_protein(protein.name):
                mapping = mutation.mappings.get(protein.referecence.name, None)
                if mapping is None:
                    continue
                pos = protein.referecence.convert_position(mapping.position)
                sample.aas[mutation.name] = protein.sequence[pos]
                if mapping.mutation_sequence in protein.sequence[pos]:
                    mutations.append(mutation)
    return mutations
