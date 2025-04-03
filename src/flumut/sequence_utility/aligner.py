from Bio.Align import PairwiseAligner, Alignment

from flumut.db_utility.db_data import references_by_segment
from flumut.sequence_utility.models import NucleotideSequence, ReferenceSequence


def align(sequence: str, segment: str) -> NucleotideSequence:
    '''
    Find the best reference and creates the alignment for the given sequence.

    :param `str` sequence: The sequence to align.
    :param `str` segment: The segment of the sequence.
    :return `NucleotideSequence`: The aligned sequence.
    '''
    sequence = sequence.replace('-', '')
    references = references_by_segment(segment)
    best_score = 0

    for reference in references:
        alignment = _pairwise_alignment(reference.sequence, sequence)
        if alignment.score > best_score:
            best_referecence = ReferenceSequence(reference.segment, reference.name, alignment[0])
            best_alignment = NucleotideSequence(alignment[1], best_referecence, alignment.score)
            best_score = alignment.score
    return best_alignment


def _pairwise_alignment(reference: str, sample: str) -> Alignment:
    '''
    Align sequence against a reference.

    :param `str` ref: Reference to align on.
    :param `str` sample: Sequence to align.
    :return `Alignment`: Best alignment.
    '''
    aligner = PairwiseAligner()
    aligner.mismatch_score = 0
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2
    aligner.query_left_open_gap_score = 1
    aligner.query_right_open_gap_score = 1

    return aligner.align(reference, sample)[0]
