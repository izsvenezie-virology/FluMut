import itertools
from copy import deepcopy
from typing import List, Optional, Tuple

from flumut.db_utility.db_data import annotations_by_reference
from flumut.exceptions import UnknownNucleotideException
from flumut.sequence_utility.models import (AminoAcidSequence,
                                            NucleotideSequence)


def translate(alignment: NucleotideSequence) -> List[AminoAcidSequence]:
    '''
    Translates an alignment into all proteins for the segment.

    :param `NucleotideSequence` alignemnt: The alignment to translate.
    :return `List[AminoAcidSequence]`: The sequence with translated proteins.
    '''
    proteins = []
    for protein_name, annotations in annotations_by_reference(alignment.referecence.name).items():
        sample_cds, reference_cds = _get_cds(alignment, annotations)

        reference = deepcopy(alignment.referecence)
        reference.sequence = ''.join(_translate_sequence(reference_cds))

        protein_sequence = _translate_sequence(sample_cds)
        protein = AminoAcidSequence(protein_name, protein_sequence, reference)
        proteins.append(protein)
    return proteins


def _get_cds(alignment: NucleotideSequence, cds: List[Tuple[int, int]]) -> Tuple[str, str]:
    '''
    Cut and assemble the nucleotide sequences based on positions given by the cds.

    :param `NucleotideSequence` alignment: Sequence alignment.
    :param `List[Tuple[int, int]]` cds: The list of starting and ending points of coding sequence.
    :return `str`, `str`: sequence coding sequence, reference coding sequence 
    '''
    cds.sort(key=lambda x: x[0])
    seq_cds = ''
    ref_cds = ''

    for rng in cds:
        start = alignment.referecence.convert_position(rng[0])
        end = alignment.referecence.convert_position(rng[1]) + 1

        ref_cds += alignment.referecence.sequence[start:end]
        seq_cds += alignment.sequence[start:end]
    return seq_cds, ref_cds


def _translate_sequence(seq: str) -> List[str]:
    '''
    Translate a nucleotide sequence in an amino acid sequence.

    :param `str` seq: Aligned nucleotide sequence to translate.
    :return `str`: Amino acid translated sequence.
    '''
    nts = list(seq)
    aas = []
    is_first = True

    for i in range(0, len(nts), 3):
        codon = _get_codon(nts, i, is_first)
        aa = _translate_codon(codon)
        if aa not in ('-', '?'):
            is_first = False
        aas.append(aa)
    return aas


def _get_codon(seq: List[str], start: int, is_first: bool) -> List[str]:
    '''
    Extract the codon from a sequence.

    :param `List[str] seq`: Sequence splitted by nucleotides.
    :param `int` start: Start position of the codon.
    :param `bool` is_first: If `true` codons starting with `-` are returned as is. 
    If `false` a frameshift is created.
    :return `str`, `str`, `str`: The codon splitted by nucleotides. 
    '''
    codon = seq[start:start + 3]
    if codon == ['-', '-', '-']:  # If the codon is a deletion
        return codon
    # If the codon starts from mid codon (to avoid frameshifts in truncated sequences):
    if is_first and '-' in codon:
        return codon
    if 'N' in codon:
        return codon
    codon = [n for n in codon if n != '-']
    while len(codon) < 3:
        next_nucl = _find_next_nucl(seq, start)
        if not next_nucl:
            break
        codon.append(seq[next_nucl])
        seq[next_nucl] = '-'
    return codon


def _translate_codon(codon: List[str]) -> str:
    '''
    Translate a codon into a set of AAs, containing all possible combinations in case of degenerations.

    :param `List[str]` codon: The codon to translate.
    :return `str`: All possible AAs concatenated.
    '''
    if 'N' in codon:
        return '?'
    try:
        undegenerated_codon = [_degeneration_dict[nucl] for nucl in codon]
    except KeyError:
        raise UnknownNucleotideException(''.join(codon)) from None
    codons = list(itertools.product(*undegenerated_codon))
    aas = [_translation_dict.get(''.join(c), '?') for c in codons]
    return ''.join(sorted(set(aas)))


def _find_next_nucl(seq: List[str], start: int) -> Optional[int]:
    '''
    Return the position of the next non deleted nucleotide.

    :param `List[str]` seq: The sequence splitted by nucleotide.
    :param `int` start: The position where to start to search for nucleotide.
    :return `int`|`None`: The position of next non-deleted nucleotide. If no nucleotide is fuound it returns `None`.
    '''
    for i in range(start + 3, len(seq)):
        if not seq[i] == '-':
            return i
    return None


'''Codon to translated amminoacid'''
_translation_dict = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    '---': '-'
}

'''Degenerated nucleotide to list of possible nucleotides'''
_degeneration_dict = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'U': ['T'], '-': ['-'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}
