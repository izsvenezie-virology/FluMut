import logging
from typing import List, Optional

from Bio.Align import PairwiseAligner
from flumutdb import Reference, Segment

from flumut.sequence_utility.models import NucleotideSequence, ReferenceSequence

WILDCARD = 'N'  # Wildcard character for not known nucleotides
MISMATCH_SCORE = -1  # Score for a mismatch between the reference and the sample sequence
GAP_OPEN_SCORE = -5  # Score for opening a gap in the alignment
GAP_EXTEND_SCORE = -1  # Score for extending a gap in the alignment
GAP_END_SCORE = -0.1  # Score for gaps at the end of the alignment, must be low to allow for partial sequences

_aligner = PairwiseAligner()
_aligner.wildcard = WILDCARD
_aligner.mismatch_score = MISMATCH_SCORE
_aligner.open_gap_score = GAP_OPEN_SCORE
_aligner.extend_gap_score = GAP_EXTEND_SCORE
_aligner.query_end_gap_score = GAP_END_SCORE
_aligner.target_end_gap_score = GAP_END_SCORE


def select_candidate_references(segment_name: Optional[str]) -> List[Reference]:
    candidates: List[Reference] = []
    for segment in Segment.all():
        if segment.name == segment_name:
            for reference in segment.references:
                candidates.append(reference)
    if not candidates:
        for segment in Segment.all():
            for reference in segment.references:
                candidates.append(reference)
    return candidates


def align(sequence: str, candidates: List[Reference]) -> NucleotideSequence:
    """
    Find the best reference and creates the alignment for the given sequence.

    :param `str` sequence: The sequence to align.
    :param `str` segment: The segment of the sequence.
    :return `NucleotideSequence`: The aligned sequence.
    """
    query_sequence = sequence.replace('-', '')
    best_score = -1_000_000
    best_reference: Optional[ReferenceSequence] = None
    best_alignment: Optional[NucleotideSequence] = None

    for reference in candidates:
        reference_sequence = str(reference.sequence)
        alignment = _aligner.align(reference_sequence, query_sequence)[0]
        if alignment.score > best_score:  # type: ignore
            best_score = alignment.score  # type: ignore
            best_reference = ReferenceSequence(str(reference.segment), str(reference.name), str(alignment[0]))
            best_alignment = NucleotideSequence(str(alignment[1]), best_reference, alignment.score)  # type: ignore

    if best_reference is None or best_alignment is None:
        raise ValueError('No candidates provided for alignment')

    logging.debug(f'Best reference: {best_reference.name} (segment {best_reference.segment})')
    logging.debug(f'reference: {best_reference.sequence}')
    logging.debug(f'sample:    {best_alignment.sequence}')
    return best_alignment
