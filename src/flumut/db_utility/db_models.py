from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict


@dataclass
class Reference:
    '''Reference sequence from DB.'''
    segment: str
    '''Segment name.'''
    name: str
    '''Name of the reference.'''
    sequence: str
    '''Nucleotide sequence.'''


@dataclass
class Mutation:
    '''Mutations data not connected to reference sequences.'''
    name: str
    '''Name of the mutation.'''
    type: str
    '''Type of the mutation.'''
    protein_name: str
    '''Name of protein.'''

    mappings: Dict[str, MutationMapping] = field(default_factory=dict, init=False)
    '''Collection of mappings for different reference sequences.'''
    aas_in_samples: Dict[str, str] = field(default_factory=dict, init=False)
    '''Collection of AAs of each sample at mutation position.'''

    def __hash__(self):
        return hash(self.name)


@dataclass
class MutationMapping:
    '''Mutations data specific of a reference sequence.'''
    mutation: Mutation
    '''Parent mutation.'''
    reference_name: str
    '''Name of the reference sequence.'''
    position: int
    '''Position of the mutation on the reference sequence.'''
    mutation_sequence: str
    '''Sequence of the mutation.'''
