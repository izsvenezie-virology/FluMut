import logging

from Bio.Align import PairwiseAligner
from flumutdb import Reference

from flumut.sequence_utility.models import NucleotideSequence, ReferenceSequence

WILDCARD = 'N'  # Wildcard character for not known nucleotides
GAP_SYMBOL = '-'  # Symbol for gaps in the alignment
MISMATCH_SCORE = -1  # Score for a mismatch between the reference and the sample sequence
GAP_OPEN_SCORE = -5  # Score for opening a gap in the alignment
GAP_EXTEND_SCORE = -1  # Score for extending a gap in the alignment
GAP_END_SCORE = -0.1  # Score for gaps at the end of the alignment, must be low to allow for partial sequences

_aligner = PairwiseAligner()
_aligner.wildcard = WILDCARD
_aligner.mismatch_score = MISMATCH_SCORE
_aligner.open_gap_score = GAP_OPEN_SCORE
_aligner.extend_gap_score = GAP_EXTEND_SCORE
_aligner.end_deletion_score = GAP_END_SCORE
_aligner.end_insertion_score = GAP_END_SCORE


def select_candidate_references(candidate_hint: str | None) -> list[Reference]:
    if not candidate_hint:
        return Reference.all()

    candidates: list[Reference] = []

    for reference in Reference.all():
        if _is_candidate_reference(reference, candidate_hint):
            candidates.append(reference)

    if not candidates:
        candidates = Reference.all()

    return candidates


def find_best_reference(sequence: str, candidates: list[Reference]) -> NucleotideSequence:
    """
    Find the best reference and creates the alignment for the given sequence.

    :param `str` sequence: The sequence to align.
    :param `str` segment: The segment of the sequence.
    :return `NucleotideSequence`: The aligned sequence.
    """
    query_sequence = sequence.replace(GAP_SYMBOL, '')
    best_score = -1_000_000
    best_reference: ReferenceSequence | None = None
    best_alignment: NucleotideSequence | None = None

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


def _is_candidate_reference(reference: Reference, candidate_hint: str) -> bool:
    if reference.name == candidate_hint:
        return True
    if reference.segment.name == candidate_hint:
        return True
    return False
