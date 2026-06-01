from flumut.core.globals import GAP_SYMBOL, WILDCARD
from flumut.core.nucleotides.models import CDS, Alignment


def adjust_frame(cds: CDS) -> None:
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
    for i, nt in enumerate(sequence[pos + 1 :]):
        if nt != GAP_SYMBOL:
            sequence[pos] = nt
            sequence[pos + 1 + i] = GAP_SYMBOL
            return True
    return False


def compensate_gap(alignment: Alignment, pos: int) -> None:
    alignment.reference.pop(pos)
    alignment.query.pop(pos)


def is_frameshift(codon: list[str]) -> bool:
    return codon.count(GAP_SYMBOL) % 3 != 0


def get_start_position(query: list[str]) -> int:
    for i in range(0, len(query), 3):
        codon = query[i : i + 3]
        if WILDCARD not in codon and GAP_SYMBOL not in codon:
            return i
    raise ValueError('No valid codon found in query sequence')
