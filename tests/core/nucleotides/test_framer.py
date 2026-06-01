from unittest.mock import MagicMock

from flumut.core.nucleotides.framer import adjust_frame
from flumut.core.nucleotides.models import CDS, Alignment


def _make_cds(reference_chars: list[str], query_chars: list[str]) -> CDS:
    return CDS(
        nucleotides=MagicMock(),
        protein=MagicMock(),
        alignment=Alignment(reference=reference_chars, query=query_chars),
    )


# ---------------------------------------------------------------------------
# No frameshifts — alignment must be left unchanged
# ---------------------------------------------------------------------------


def test_adjust_frame_without_frameshifts() -> None:
    ref = list('ATGAAA')
    qry = list('GTGCAT')
    cds = _make_cds(ref, qry)
    adjust_frame(cds)
    assert cds.alignment.reference == ref
    assert cds.alignment.query == qry
    assert cds.frameshifts == []


# ---------------------------------------------------------------------------
# Opening gaps in query — must not be moved (treated as leading Ns)
# ---------------------------------------------------------------------------


def test_adjust_frame_opening_gaps_in_query_are_left_untouched() -> None:
    # Query does not cover the first codon; those gaps are NOT frameshifts
    ref = list('ATGAAGGGA')
    qry = list('---AAGGGA')
    cds = _make_cds(ref, qry)
    adjust_frame(cds)
    assert cds.alignment.reference == list('ATGAAGGGA')
    assert cds.alignment.query == list('---AAGGGA')
    assert cds.frameshifts == []


# ---------------------------------------------------------------------------
# Unresolved insertions in query (gap in reference)
# Gap moves forward until end of sequence.
# Frameshift is recorded at the last codon where the gap lands (aa 3).
# ---------------------------------------------------------------------------


def test_adjust_frame_insertion_at_codon_boundary_moves_to_end() -> None:
    # 1-nt insertion: ref gap at codon boundary (pos 3)
    cds = _make_cds(list('ATG-AAGGG'), list('ATGCAAGGG'))
    adjust_frame(cds)
    assert cds.alignment.reference == list('ATGAAGGG-')
    assert cds.alignment.query == list('ATGCAAGGG')
    assert cds.frameshifts == [(4, None)]


def test_adjust_frame_insertion_mid_codon_moves_to_end() -> None:
    # 1-nt insertion: ref gap inside first codon (pos 2)
    cds = _make_cds(list('AT-GAAGGG'), list('ATCGAAGGG'))
    adjust_frame(cds)
    assert cds.alignment.reference == list('ATGAAGGG-')
    assert cds.alignment.query == list('ATCGAAGGG')
    assert cds.frameshifts == [(3, None)]


# ---------------------------------------------------------------------------
# Unresolved deletions in query (gap in query)
# Gap moves forward until end of sequence.
# ---------------------------------------------------------------------------


def test_adjust_frame_deletion_at_codon_boundary_moves_to_end() -> None:
    # 1-nt deletion: qry gap at codon boundary (pos 3)
    cds = _make_cds(list('ATGAAGGGC'), list('ATG-AGGGC'))
    adjust_frame(cds)
    assert cds.alignment.reference == list('ATGAAGGGC')
    assert cds.alignment.query == list('ATGAGGGC-')
    assert cds.frameshifts == [(4, None)]


def test_adjust_frame_two_nt_deletion_moves_to_end() -> None:
    # 2-nt deletion: two qry gaps spread across codons → both move to end
    cds = _make_cds(list('ATGAAGGGATTT'), list('ATGAA--GATTT'))
    adjust_frame(cds)
    assert cds.alignment.reference == list('ATGAAGGGATTT')
    assert cds.alignment.query == list('ATGAAGATTT--')
    assert cds.frameshifts == [(6, None)]


# ---------------------------------------------------------------------------
# Compensation — one ref gap and one qry gap cancel each other
# The column where both are '-' is removed; both sequences shrink by 1.
# ---------------------------------------------------------------------------


def test_adjust_frame_compensation_removes_column_and_shrinks_sequences() -> None:
    # ref gap at pos 3 (+1 frame), qry gap at pos 5 (-1 frame) → cancel
    cds = _make_cds(list('ATG-AAGGG'), list('ATGCA-GGG'))
    adjust_frame(cds)
    assert cds.alignment.reference == list('ATGAAGGG')  # shrank from 9 → 8
    assert cds.alignment.query == list('ATGCAGGG')  # shrank from 9 → 8
    assert cds.frameshifts == [(4, 6)]


# ---------------------------------------------------------------------------
# Multiple gaps forming a complete deletion codon
# Three spread gaps collapse into a single '---' codon.
# ---------------------------------------------------------------------------


def test_adjust_frame_three_nt_deletion_collapses_to_gap_codon() -> None:
    # qry has 3 gaps at positions 3, 5, 7 → collapse to codon slot 6:9
    cds = _make_cds(list('ATGAAGGGATTT'), list('ATG-A-G-ATTT'))
    adjust_frame(cds)
    assert cds.alignment.reference == list('ATGAAGGGATTT')
    assert cds.alignment.query == list('ATGAGA---TTT')
    assert cds.frameshifts == [(4, 6)]
