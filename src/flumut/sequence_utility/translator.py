import itertools
from copy import deepcopy
from typing import List, Optional, Tuple

from flumut.db_utility.db_data import annotations_by_reference
from flumut.sequence_utility.exceptions import UnknownNucleotideException
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
        sequence_aa, referenece_aa, frameshifts = _translate_sequence(sample_cds, reference_cds)

        reference = deepcopy(alignment.referecence)
        reference.sequence = ''.join(referenece_aa)
        protein = AminoAcidSequence(protein_name, sequence_aa, reference)
        protein.frameshifts = frameshifts
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


def _translate_sequence(seq: str, ref: str) -> Tuple[List[str], List[str], List[Tuple[int, int]]]:
    '''
    Translate a nucleotide sequence in an amino acid sequence.

    :param `str` seq: Aligned nucleotide sequence to translate.
    :return `List[str]`: Amino acid translated sequence of the sample.
    :return `List[str]`: Amino acid translated sequence of the reference.
    :return `List[Tuple[int, int]]`: List of frameshifts.
    '''
    ref_nts = list(ref)
    seq_nts = list(seq)

    ref_aas = []
    seq_aas = []

    frameshift_start = None
    frameshifts = []
    is_first = True

    for i in range(0, len(seq_nts), 3):
        codon, seq_frameshift = _get_codon(seq_nts, i, is_first)
        aa = _translate_codon(codon)
        if aa not in ('-', '?'):
            is_first = False
        seq_aas.append(aa)

        codon, ref_frameshift = _get_codon(ref_nts, i)
        aa = _translate_codon(codon)
        ref_aas.append(aa)

        if seq_frameshift is None or ref_frameshift is None:
            continue
        if not seq_frameshift - ref_frameshift == 0 and frameshift_start is None:
            frameshift_start = int(i/3 + 1)
            continue
        if seq_frameshift - ref_frameshift == 0 and frameshift_start is not None:
            frameshift_end = int(i/3)
            if frameshift_end - frameshift_start > 1:
                frameshifts.append((frameshift_start, frameshift_end))
            frameshift_start = None
            continue

    if frameshift_start is not None:
        frameshift_end = len(seq_aas)
        if frameshift_end - frameshift_start > 1:
            frameshifts.append((frameshift_start, frameshift_end))
    return seq_aas, ref_aas, frameshifts


def _get_codon(seq: List[str], start: int, is_first: bool = False) -> Tuple[List[str], int]:
    '''
    Extract the codon from a sequence.

    :param `List[str] seq`: Sequence splitted by nucleotides.
    :param `int` start: Start position of the codon.
    :param `bool` is_first: If `true` codons containing `-` are returned as is. 
    If `false` a frameshift is created.
    :return `str`, `str`, `str` codon: The codon splitted by nucleotides.
    :return `int` frameshift_status: Number of nucleotides moved in order to create the codon.
    '''
    codon = seq[start:start + 3]
    if codon == ['-', '-', '-']:  # If the codon is a deletion
        return codon, None
    # If the codon starts from mid codon (to avoid frameshifts in truncated sequences):
    if is_first and '-' in codon:
        return codon, None
    if 'N' in codon:
        return codon, None
    codon = [n for n in codon if n != '-']
    frameshift_status = 3 - len(codon)
    while len(codon) < 3:
        next_nucl = _find_next_nucl(seq, start)
        if not next_nucl:
            break
        codon.append(seq[next_nucl])
        seq[next_nucl] = '-'
    return codon, frameshift_status


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
'''Codon to translated amminoacid'''

_degeneration_dict = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'U': ['T'], '-': ['-'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}
'''Degenerated nucleotide to list of possible nucleotides'''
