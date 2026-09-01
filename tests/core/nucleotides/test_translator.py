from unittest.mock import MagicMock, patch

import pytest

from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL
from flumut.core.nucleotides.exceptions import UnknownNucleotideError
from flumut.core.nucleotides.models import CDS, Alignment
from flumut.core.nucleotides.translator import get_cds, translate, translate_codon

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_cds(reference: str, query: str) -> CDS:
    return CDS(
        nucleotides=MagicMock(),
        protein=MagicMock(),
        alignment=Alignment(reference=list(reference), query=list(query)),
    )


def _make_nucleotides(reference: str, query: str, annotations: tuple = ()) -> tuple[MagicMock, MagicMock]:
    nucleotides = MagicMock()
    nucleotides.alignment = Alignment(reference=list(reference), query=list(query))
    protein = MagicMock()
    # get_cds reads the annotations pre-grouped per (reference, protein).
    nucleotides.reference.annotations_by_protein = {
        protein: [MagicMock(start=start, end=end) for start, end in annotations]
    }
    return nucleotides, protein


@pytest.fixture
def mock_adjust_frame():
    with patch('flumut.core.nucleotides.translator.adjust_frame') as mock:
        yield mock


# ---------------------------------------------------------------------------
# translate
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'reference, query, expected_reference, expected_query',
    [
        # ATG->M AAA->K ; GTG->V CAT->H
        ('ATGAAA', 'GTGCAT', ['M', 'K'], ['V', 'H']),
        # R = A|G, so RTG translates to both M and V, stored as one element.
        ('ATG', 'RTG', ['M'], ['MV']),
    ],
    ids=['plain_codons', 'degenerate_codon_stays_one_element'],
)
def test_translate(reference: str, query: str, expected_reference: list[str], expected_query: list[str]) -> None:
    """Each codon becomes exactly one element, even when it resolves to several amino acids."""
    result = translate(_make_cds(reference, query))
    assert result.alignment.reference == expected_reference
    assert result.alignment.query == expected_query


def test_translate_result_carries_over_the_cds_fields() -> None:
    cds = _make_cds('ATG', 'GTG')
    result = translate(cds)
    assert result.protein is cds.protein
    assert result.reference is cds.nucleotides.reference
    assert result.cds is cds


# ---------------------------------------------------------------------------
# get_cds
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'reference, query, annotations, expected_reference, expected_query',
    [
        ('ACT', 'KMR', (), '', ''),
        ('ACTG', 'KMRV', ((2, 3),), 'CT', 'MR'),
        ('ACTGA', 'KMRVH', ((1, 1), (4, 5)), 'AGA', 'KVH'),
        # A leading gap shifts the columns but not the reference positions.
        (f'{GAP_SYMBOL}ACT', 'XKMR', ((2, 3),), 'CT', 'MR'),
    ],
    ids=['no_annotations', 'single_annotation', 'annotations_concatenated', 'gap_before_cds'],
)
def test_get_cds_slices_by_reference_position(
    mock_adjust_frame, reference: str, query: str, annotations: tuple, expected_reference: str, expected_query: str
) -> None:
    """Annotations select alignment columns by 1-based reference position, gaps excluded."""
    nucleotides, protein = _make_nucleotides(reference, query, annotations)
    cds = get_cds(nucleotides, protein)
    assert cds.alignment.reference == list(expected_reference)
    assert cds.alignment.query == list(expected_query)


def test_get_cds_wraps_its_inputs_and_adjusts_the_frame(mock_adjust_frame) -> None:
    nucleotides, protein = _make_nucleotides('ACT', 'KMR')
    cds = get_cds(nucleotides, protein)
    assert cds.nucleotides is nucleotides
    assert cds.protein is protein
    mock_adjust_frame.assert_called_once_with(cds)


# ---------------------------------------------------------------------------
# translate_codon
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'codon, expected',
    [
        (['N', 'T', 'G'], UNKNOWN_AA_SYMBOL),
        (['A', 'N', 'G'], UNKNOWN_AA_SYMBOL),
        ([], UNKNOWN_AA_SYMBOL),
        (['A', 'T'], UNKNOWN_AA_SYMBOL),
        (['A', 'T', 'G'], 'M'),
        (['T', 'T', 'T'], 'F'),
        (['T', 'A', 'A'], STOP_CODON_SYMBOL),
        ([GAP_SYMBOL, GAP_SYMBOL, GAP_SYMBOL], GAP_SYMBOL),
        (['R', 'T', 'G'], 'MV'),
        (['T', 'C', 'Y'], 'S'),
        (['A', 'T', GAP_SYMBOL], UNKNOWN_AA_SYMBOL),
    ],
    ids=[
        'wildcard_at_start',
        'wildcard_in_middle',
        'empty',
        'two_bases',
        'methionine',
        'phenylalanine',
        'stop_codon',
        'gap_codon',
        'degenerate_two_aas',
        'degenerate_same_aa',
        'partial_gap',
    ],
)
def test_translate_codon(codon: list[str], expected: str) -> None:
    assert translate_codon(codon) == expected


def test_translate_codon_unknown_nucleotide_raises() -> None:
    with pytest.raises(UnknownNucleotideError):
        translate_codon(['X', 'T', 'G'])
