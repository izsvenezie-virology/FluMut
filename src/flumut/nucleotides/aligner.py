import logging

from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord

from flumut.core.globals import GAP_END_SCORE, GAP_EXTEND_SCORE, GAP_OPEN_SCORE, GAP_SYMBOL, MISMATCH_SCORE, WILDCARD
from flumut.flumutdb import Reference
from flumut.nucleotides.models import Alignment, NucleotideAlignment

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


def get_best_alignment(query: SeqRecord, candidates: list[Reference]) -> NucleotideAlignment:
    """
    Find the best reference and creates the alignment for the given sequence.

    :param `str` sequence: The sequence to align.
    :param `str` segment: The segment of the sequence.
    :return `NucleotideSequence`: The aligned sequence.
    """
    best_alignment = None
    best_score = -1_000_000
    query_sequence = str(query.seq).replace(GAP_SYMBOL, '')

    for reference in candidates:
        reference_sequence = str(reference.sequence)
        alignment = _aligner.align(reference_sequence, query_sequence)[0]
        if best_score < alignment.score:  # type: ignore
            aln = Alignment(list(alignment[0]), list(alignment[1]))  # type: ignore
            best_alignment = NucleotideAlignment(query, reference, aln)
            best_score = alignment.score  # type: ignore

    if not best_alignment:
        raise ValueError(f'No reference found for sequence {query.name}')

    logging.debug(f'Best reference: {best_alignment.reference.name} (segment {best_alignment.reference.segment})')
    return best_alignment


def _is_candidate_reference(reference: Reference, candidate_hint: str) -> bool:
    if reference.name == candidate_hint:
        return True
    if reference.segment.name == candidate_hint:
        return True
    return False
