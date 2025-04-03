from dataclasses import dataclass, field
from typing import Dict, List

from flumut.db_utility.db_models import Mutation


@dataclass
class ReferenceSequence:
    segment: str
    '''The segment the reference belongs to.'''
    name: str
    '''Reference name.'''
    sequence: str
    '''Aligned sequence of reference.'''

    def convert_position(self, position: int) -> int:
        '''
        Converts theorical position to actual position based on gaps in the reference sequence.
        Coverts also 1-based positions to 0-based.

        :param `int` position: Theorical 1-based position to convert.
        :return `int`: Actual 0-based converted position.
        '''
        position -= 1    # conversion 1-based to 0-based numeration
        dashes = 0
        adj_pos = position
        while self.sequence.count('-', 0, adj_pos + 1) != dashes:
            dashes = self.sequence.count('-', 0, adj_pos + 1)
            adj_pos = position + dashes
        return adj_pos


@dataclass
class AminoAcidSequence():
    '''Stores aligned amminoacidic sequences.'''
    name: str
    '''Protein name.'''
    sequence: List[str]
    '''Aligned proteic sequence of the sample as list, one element per position. 
        Each position is a string of all AAs that may be present.'''
    referecence: ReferenceSequence
    '''The amino acid reference.'''

    mutations: List[Mutation] = field(default_factory=list, init=False)
    '''All mutations found in the sample.'''


@dataclass
class NucleotideSequence():
    '''Aligned nucleotide sequence.'''
    sequence: str
    '''Aligned sequence.'''
    referecence: ReferenceSequence
    '''The nucleotide aligned reference.'''
    score: float
    '''Score of the alignment, the higher the better.'''

    proteins: List[AminoAcidSequence] = field(default_factory=list, init=False)
    '''List of all proteins translated from the sequence.'''


@dataclass
class FastaSequence():
    '''Sequence from Fasta file.'''
    file: str
    '''Name of Fasta file storing this sequence.'''
    header: str
    '''Header of sequence in Fasta file.'''
    sequence: str = ''
    '''Nucleotidic sequence.'''

    sample: str = field(default='', init=False)
    '''Sample name.'''
    segment: str = field(default='', init=False)
    '''Segment name retrieved in header.'''

    alignment: NucleotideSequence = field(default=None, init=False)
    '''Aligned nucleotide sequence.'''


@dataclass
class Sample:
    name: str
    '''Sample name.'''

    sequences: List[FastaSequence] = field(default_factory=list, init=False)
    '''Sequences related to this sample.'''
    aas: Dict[str, str] = field(default_factory=dict, init=False)
    '''Amino acids present at position of each mutation.
        Keys are the mutation names, values are the amino acids present in the sample.'''
    mutations: List[Mutation] = field(default_factory=list, init=False)
    '''Mutations present in the sample.'''
    markers: List[Dict[str, str]] = field(default_factory=list, init=False)
    '''Markers found for the sample.'''
