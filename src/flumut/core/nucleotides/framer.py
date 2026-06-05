from flumut.core.globals import GAP_SYMBOL, WILDCARD
from flumut.core.nucleotides.models import CDS, Alignment


def adjust_frame(cds: CDS) -> None:
    """Resolve frameshift indels in a CDS alignment by sliding gaps out of codons.

    Iterates over the alignment in codon-sized steps. When a gap causes a
    frameshift within a codon, the gap is shifted to the nearest codon boundary.
    Paired reference/query gaps are removed entirely. Frameshift ranges that
    could not be closed are recorded in ``cds.frameshifts``.

    Args:
        cds: The CDS whose alignment will be modified in-place.
    """
    aln = cds.alignment

    frameshifts = set()
    non_closing_frameshift = False

    pos = get_start_position(aln.query)

    while pos < len(aln.reference):
        codon_reference = aln.reference[pos : pos + 3]
        codon_query = aln.query[pos : pos + 3]
        moved = False
        if is_frameshift(codon_reference) or is_frameshift(codon_query):
            for i in range(3):
                if codon_reference[i] == GAP_SYMBOL and codon_query[i] == GAP_SYMBOL:
                    compensate_gap(aln, pos + i)
                    moved = True
                    break

                if codon_reference[i] == GAP_SYMBOL:
                    moved = move_gap(aln.reference, pos + i)
                    break
                if codon_query[i] == GAP_SYMBOL:
                    moved = move_gap(aln.query, pos + i)
                    break
            if moved:
                frameshifts.add(pos + i + 1)  # type: ignore[reportPossiblyUnboundVariable]
            else:
                non_closing_frameshift = True
                pos += 3
        else:
            pos += 3
    cds.frameshifts = _collapse_to_ranges(sorted(frameshifts), non_closing_frameshift)


def _collapse_to_ranges(positions: list[int], non_closing_frameshift: bool) -> list[tuple[int, int | None]]:
    """Collapse a sorted list of frameshift nucleotide positions into contiguous ranges.

    Consecutive positions are merged into a single range; a gap of more than
    9 nucleotides starts a new range. The last range's end position is set to
    None when ``non_closing_frameshift`` is True, indicating it extends to the
    end of the sequence.

    Args:
        positions: Sorted list of 1-based nucleotide positions with a frameshift.
        non_closing_frameshift: When True, the final range has no closing position.

    Returns:
        A list of ``(start, end)`` tuples; ``end`` is None for an open-ended frameshift.
    """
    if not positions:
        return []
    result = []
    start = prev = positions[0]
    for p in positions[1:]:
        if p > prev + 9:
            result.append((start, prev))
            start = p
        prev = p
    last = None if non_closing_frameshift else prev
    result.append((start, last))
    return result


def move_gap(sequence: list[str], pos: int) -> bool:
    """Slide the gap symbol at ``pos`` forward to the next non-gap position.

    Replaces ``sequence[pos]`` with the first non-gap nucleotide found after it,
    then places the gap symbol at that nucleotide's original position.

    Args:
        sequence: Mutable list of nucleotide characters (modified in-place).
        pos: Index of the gap symbol to move.

    Returns:
        True if a non-gap nucleotide was found and the swap was performed;
        False if no non-gap nucleotide exists after ``pos``.
    """
    for i, nt in enumerate(sequence[pos + 1 :]):
        if nt != GAP_SYMBOL:
            sequence[pos] = nt
            sequence[pos + 1 + i] = GAP_SYMBOL
            return True
    return False


def compensate_gap(alignment: Alignment, pos: int) -> None:
    """Remove a paired reference/query gap at the given position.

    Both ``alignment.reference`` and ``alignment.query`` have their element at
    ``pos`` removed in-place, shortening both lists by one.

    Args:
        alignment: The Alignment whose lists are modified in-place.
        pos: Index of the gap column to remove from both lists.
    """
    alignment.reference.pop(pos)
    alignment.query.pop(pos)


def is_frameshift(codon: list[str]) -> bool:
    """Return True if the number of gaps in a codon is not a multiple of three.

    A gap count that is not divisible by 3 indicates an insertion or deletion
    that would shift the reading frame.

    Args:
        codon: A list of up to three nucleotide characters.

    Returns:
        True if the gap count in ``codon`` is not divisible by 3.
    """
    return codon.count(GAP_SYMBOL) % 3 != 0


def get_start_position(query: list[str]) -> int:
    """Find the first codon-aligned position that contains no wildcards or gaps.

    Scans the query in 3-nucleotide steps and returns the offset of the first
    complete, unambiguous codon.

    Args:
        query: List of nucleotide characters representing the query alignment.

    Returns:
        The 0-based index of the first character of the first valid codon.

    Raises:
        ValueError: If no valid codon is found in the entire query sequence.
    """
    for i in range(0, len(query), 3):
        codon = query[i : i + 3]
        if WILDCARD not in codon and GAP_SYMBOL not in codon:
            return i
    raise ValueError('No valid codon found in query sequence')
