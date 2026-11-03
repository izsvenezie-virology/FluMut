from Bio.Align import PairwiseAligner

from flumut.core.globals import GAP_END_SCORE, GAP_EXTEND_SCORE, GAP_OPEN_SCORE, GAP_SYMBOL, MISMATCH_SCORE, WILDCARD
from flumut.core.io.input import FastaSequence
from flumut.core.logger import LOGGER
from flumut.core.nucleotides.exceptions import NoReferenceFoundError
from flumut.core.nucleotides.models import Alignment, Nucleotide
from flumut.flumutdb import Reference, loader

_aligner = PairwiseAligner()
_aligner.wildcard = WILDCARD
_aligner.mismatch_score = MISMATCH_SCORE
_aligner.open_gap_score = GAP_OPEN_SCORE
_aligner.extend_gap_score = GAP_EXTEND_SCORE
_aligner.end_deletion_score = GAP_END_SCORE
_aligner.end_insertion_score = GAP_END_SCORE


def select_candidate_references(candidate_hint: str | None) -> list[Reference]:
    """Return the subset of references matching a segment or reference name hint.

    If ``candidate_hint`` is None or no reference matches, all known references
    are returned as candidates.

    Args:
        candidate_hint: A reference name or segment name to filter by, or None
            to use all references.

    Returns:
        A non-empty list of Reference objects to consider during alignment.
    """
    references = loader.get(Reference)
    if not candidate_hint:
        return references

    candidates: list[Reference] = []

    for reference in references:
        if _is_candidate_reference(reference, candidate_hint):
            candidates.append(reference)
    return candidates or references


def get_best_alignment(query: FastaSequence, candidates: list[Reference]) -> Nucleotide:
    """Align a query sequence against all candidate references and return the best match.

    Uses pairwise alignment scores to select the highest-scoring reference.
    The returned Nucleotide bundles the query, the winning reference, and the
    alignment.

    Args:
        query: The FASTA sequence to align.
        candidates: Reference sequences to align against.

    Returns:
        A Nucleotide containing the query, the best-scoring Reference, and the
        resulting Alignment.

    Raises:
        NoReferenceFoundError: If ``candidates`` is empty and no alignment can be produced.
    """
    if not candidates:
        raise NoReferenceFoundError(f'No reference found for sequence {query.header} in file {query.file}')

    scores: dict[Reference, float] = {}
    query_sequence = str(query.sequence).replace(GAP_SYMBOL, '')

    LOGGER.debug(f'Aligning "{query.header}" against {len(candidates)} candidate reference(s)...')

    for reference in candidates:
        scores[reference] = _aligner.score(reference.sequence, query_sequence)

    best_reference: Reference = max(scores, key=scores.get)
    alignment = _aligner.align(best_reference.sequence, query_sequence)[0]
    aln = Alignment(list(alignment[0]), list(alignment[1]))  # type: ignore
    LOGGER.debug(f'Best reference: {best_reference.name} (segment {best_reference.segment})')

    return Nucleotide(query, best_reference, aln)


def _is_candidate_reference(reference: Reference, candidate_hint: str) -> bool:
    """Check whether a reference matches a given name or segment hint.

    Args:
        reference: The Reference to test.
        candidate_hint: A reference name or segment name to match against.

    Returns:
        True if the reference name or its segment name equals ``candidate_hint``.
    """
    if reference.name == candidate_hint:
        return True
    return reference.segment.name == candidate_hint
