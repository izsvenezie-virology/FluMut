from unittest.mock import MagicMock

import pytest

from flumut.core.nucleotides.framer import adjust_frame
from flumut.core.nucleotides.models import CDS, Alignment


def _make_cds(reference: str, query: str) -> CDS:
    return CDS(
        nucleotides=MagicMock(),
        protein=MagicMock(),
        alignment=Alignment(reference=list(reference), query=list(query)),
    )


@pytest.mark.parametrize(
    'reference, query, expected_reference, expected_query, expected_frameshifts',
    [
        # Nothing to correct: the alignment is left exactly as it was.
        ('ATGAAA', 'GTGCAT', 'ATGAAA', 'GTGCAT', []),
        # Query does not cover the first codon; leading gaps are not a frameshift.
        ('ATGAAGGGA', '---AAGGGA', 'ATGAAGGGA', '---AAGGGA', []),
        # Unresolved 1-nt insertion (gap in reference) at a codon boundary: gap travels to the end.
        ('ATG-AAGGG', 'ATGCAAGGG', 'ATGAAGGG-', 'ATGCAAGGG', [(4, None)]),
        # Same insertion, but starting mid-codon.
        ('AT-GAAGGG', 'ATCGAAGGG', 'ATGAAGGG-', 'ATCGAAGGG', [(3, None)]),
        # Unresolved 1-nt deletion (gap in query): gap travels to the end.
        ('ATGAAGGGC', 'ATG-AGGGC', 'ATGAAGGGC', 'ATGAGGGC-', [(4, None)]),
        # Unresolved 2-nt deletion: both gaps travel to the end together.
        ('ATGAAGGGATTT', 'ATGAA--GATTT', 'ATGAAGGGATTT', 'ATGAAGATTT--', [(6, None)]),
        # Compensation: a reference gap and a query gap cancel, so the column is dropped.
        ('ATG-AAGGG', 'ATGCA-GGG', 'ATGAAGGG', 'ATGCAGGG', [(4, 6)]),
        # Three spread query gaps collapse into one complete '---' codon.
        ('ATGAAGGGATTT', 'ATG-A-G-ATTT', 'ATGAAGGGATTT', 'ATGAGA---TTT', [(4, 6)]),
    ],
    ids=[
        'no_frameshift',
        'leading_gaps_ignored',
        'insertion_at_codon_boundary',
        'insertion_mid_codon',
        'deletion_at_codon_boundary',
        'two_nt_deletion',
        'compensated_indel',
        'three_nt_deletion_collapses',
    ],
)
def test_adjust_frame(
    reference: str, query: str, expected_reference: str, expected_query: str, expected_frameshifts: list
) -> None:
    """Unresolved indels are pushed to the end of the CDS and recorded as frameshifts."""
    cds = _make_cds(reference, query)
    adjust_frame(cds)

    assert cds.alignment.reference == list(expected_reference)
    assert cds.alignment.query == list(expected_query)
    assert cds.frameshifts == expected_frameshifts
