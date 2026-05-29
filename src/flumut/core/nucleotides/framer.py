from flumut.core.globals import GAP_SYMBOL, WILDCARD
from flumut.core.nucleotides.models import CDS, Alignment


def adjust_frame(cds: CDS) -> None:
    aln = cds.alignment

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
            if not moved:
                pos += 3
        else:
            pos += 3


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
