from unittest.mock import MagicMock

import pytest

from flumut.core.globals import GAP_SYMBOL
from flumut.core.nucleotides.models import CDS, Alignment, Nucleotide

G = GAP_SYMBOL


# ---------------------------------------------------------------------------
# Alignment dataclass tests
# ---------------------------------------------------------------------------


def test_alignment_defaults_to_empty_lists() -> None:
    a = Alignment()
    assert a.reference == []
    assert a.query == []


def test_alignment_stores_provided_values() -> None:
    a = Alignment(reference=['A', 'C'], query=['K', 'M'])
    assert a.reference == ['A', 'C']
    assert a.query == ['K', 'M']


def test_alignment_lists_independent_per_instance() -> None:
    a1 = Alignment()
    a2 = Alignment()
    a1.reference.append('A')
    assert a2.reference == []


# ---------------------------------------------------------------------------
# get_positions tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'reference, expected',
    [
        ([], []),
        (['A'], [1]),
        ([G], [None]),
        (['A', 'C', 'T'], [1, 2, 3]),
        ([G, 'A', 'C'], [None, 1, 2]),
        (['A', 'C', G], [1, 2, None]),
        (['A', G, 'C'], [1, None, 2]),
        (['A', G, G, 'C'], [1, None, None, 2]),
        ([G, G, G], [None, None, None]),
    ],
    ids=[
        'empty',
        'single_residue',
        'single_gap',
        'no_gaps',
        'gap_at_start',
        'gap_at_end',
        'gap_in_middle',
        'consecutive_gaps',
        'all_gaps',
    ],
)
def test_get_positions(reference: list[str], expected: list) -> None:
    assert Alignment(reference=reference).get_positions() == expected


# ---------------------------------------------------------------------------
# Nucleotide tests
# ---------------------------------------------------------------------------


def test_nucleotide_stores_all_fields() -> None:
    query = MagicMock()
    reference = MagicMock()
    alignment = Alignment()
    n = Nucleotide(query=query, reference=reference, alignment=alignment)
    assert n.query is query
    assert n.reference is reference
    assert n.alignment is alignment


# ---------------------------------------------------------------------------
# CDS tests
# ---------------------------------------------------------------------------


def test_cds_stores_required_fields() -> None:
    nucleotides = MagicMock()
    protein = MagicMock()
    alignment = Alignment()
    cds = CDS(nucleotides=nucleotides, protein=protein, alignment=alignment)
    assert cds.nucleotides is nucleotides
    assert cds.protein is protein
    assert cds.alignment is alignment


def test_cds_frameshifts_defaults_to_empty_list() -> None:
    cds = CDS(nucleotides=MagicMock(), protein=MagicMock(), alignment=Alignment())
    assert cds.frameshifts == []


def test_cds_frameshifts_can_be_provided() -> None:
    shifts = [(1, 3), (7, 9)]
    cds = CDS(nucleotides=MagicMock(), protein=MagicMock(), alignment=Alignment(), frameshifts=shifts)
    assert cds.frameshifts == shifts


def test_cds_frameshifts_list_independent_per_instance() -> None:
    cds1 = CDS(nucleotides=MagicMock(), protein=MagicMock(), alignment=Alignment())
    cds2 = CDS(nucleotides=MagicMock(), protein=MagicMock(), alignment=Alignment())
    cds1.frameshifts.append((1, 3))
    assert cds2.frameshifts == []
