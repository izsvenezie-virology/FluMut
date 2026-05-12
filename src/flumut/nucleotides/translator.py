import itertools

from flumut.core.globals import GAP_SYMBOL, WILDCARD
from flumut.core.models import ProteinAlignment
from flumut.flumutdb import Protein
from flumut.nucleotides import CDSAlignment
from flumut.nucleotides.models import NucleotideAlignment


def translate(alignment: NucleotideAlignment, protein: Protein) -> ProteinAlignment:
    cds = get_cds(alignment, protein)
    cds.adjust_frame()

    aa_reference = []
    aa_query = []

    for i in range(0, len(cds.aligned_query), 3):
        codon_query = cds.aligned_query[i : i + 3]
        codon_reference = cds.aligned_reference[i : i + 3]
        aa_reference += translate_codon(codon_reference)
        aa_query += translate_codon(codon_query)

    return ProteinAlignment(protein, alignment.reference, cds=cds, aligned_reference=aa_reference, aligned_query=aa_query)


def get_cds(alignment: NucleotideAlignment, protein: Protein) -> CDSAlignment:
    annotations = protein.annotations
    nucleotides = alignment
    positions = nucleotides.get_positions()
    cds_reference = []
    cds_query = []
    for annotation in annotations:
        if annotation.reference == nucleotides.reference:
            start = positions.index(annotation.start)
            end = positions.index(annotation.end) + 1
            cds_query += nucleotides.aligned_query[start:end]
            cds_reference += nucleotides.aligned_reference[start:end]
    return CDSAlignment(alignment, protein, aligned_reference=cds_reference, aligned_query=cds_query)


def translate_codon(codon: list[str]) -> str:
    if WILDCARD in codon:
        return '?'
    if len(codon) < 3:
        return '?'

    try:
        undegenerated_codon = [_degeneration_dict[nucl] for nucl in codon]
    except KeyError:
        raise ValueError('Unknown nucleotide found in codon: ' + ''.join(codon))

    codons = list(itertools.product(*undegenerated_codon))
    aas = [_translation_dict.get(''.join(c), '?') for c in codons]
    return ''.join(sorted(set(aas)))


def get_start_position(query: list[str]) -> int:
    for i in range(0, len(query), 3):
        codon = query[i : i + 3]
        if WILDCARD not in codon and GAP_SYMBOL not in codon:
            return i
    raise ValueError('No valid codon found in query sequence')


# fmt: off
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
# fmt: on
