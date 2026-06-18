from unittest.mock import MagicMock, patch

import pytest

from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL
from flumut.core.nucleotides.exceptions import UnknownNucleotideError
from flumut.core.nucleotides.models import CDS, Alignment
from flumut.core.nucleotides.translator import get_cds, translate, translate_codon

# ---------------------------------------------------------------------------
# translate helpers and tests
# ---------------------------------------------------------------------------


def _make_cds(reference_chars: list[str], query_chars: list[str]) -> CDS:
    return CDS(
        nucleotides=MagicMock(),
        protein=MagicMock(),
        alignment=Alignment(reference=reference_chars, query=query_chars),
    )


def test_translate_produces_correct_amino_acid_sequences() -> None:
    # ATG → M, AAA → K  (reference)
    # GTG → V, CAT → H  (query)
    cds = _make_cds(
        reference_chars=['A', 'T', 'G', 'A', 'A', 'A'],
        query_chars=['G', 'T', 'G', 'C', 'A', 'T'],
    )
    result = translate(cds)
    assert result.alignment.reference == ['M', 'K']
    assert result.alignment.query == ['V', 'H']


def test_translate_each_codon_is_one_list_element() -> None:
    # append() stores each translated codon as one element, even multi-char ones
    cds = _make_cds(
        reference_chars=['A', 'T', 'G'],  # R=A|G → ATG=M or GTG=V → 'MV'
        query_chars=['R', 'T', 'G'],
    )
    result = translate(cds)
    assert result.alignment.reference == ['M']
    assert result.alignment.query == ['MV']


def test_translate_result_stores_fields_from_cds() -> None:
    cds = _make_cds(['A', 'T', 'G'], ['G', 'T', 'G'])
    result = translate(cds)
    assert result.protein is cds.protein
    assert result.reference is cds.nucleotides.reference
    assert result.cds is cds


# ---------------------------------------------------------------------------
# get_cds fixture, helpers, and tests
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_adjust_frame():
    with patch('flumut.core.nucleotides.translator.adjust_frame') as mock:
        yield mock


def _make_nucleotides(
    reference_chars: list[str],
    query_chars: list[str],
    annotations: list | None = None,
) -> tuple[MagicMock, MagicMock]:
    n = MagicMock()
    n.alignment = Alignment(reference=reference_chars, query=query_chars)
    protein = MagicMock()
    # get_cds reads the annotations pre-grouped per (reference, protein)
    n.reference.annotations_by_protein = {protein: annotations or []}
    return n, protein


def _make_annotation(start: int, end: int) -> MagicMock:
    ann = MagicMock()
    ann.start = start
    ann.end = end
    return ann


def test_get_cds_no_annotations_returns_empty_alignment(mock_adjust_frame) -> None:
    n, protein = _make_nucleotides(['A', 'C', 'T'], ['K', 'M', 'R'])
    cds = get_cds(n, protein)
    assert cds.alignment.reference == []
    assert cds.alignment.query == []


def test_get_cds_matching_annotation_extracts_correct_slice(mock_adjust_frame) -> None:
    # positions: [1, 2, 3, 4] — annotation covers positions 2–3 → indices 1:3
    ann = _make_annotation(start=2, end=3)
    n, protein = _make_nucleotides(['A', 'C', 'T', 'G'], ['K', 'M', 'R', 'V'], annotations=[ann])
    cds = get_cds(n, protein)
    assert cds.alignment.reference == ['C', 'T']
    assert cds.alignment.query == ['M', 'R']


def test_get_cds_multiple_annotations_concatenated(mock_adjust_frame) -> None:
    # positions: [1, 2, 3, 4, 5]
    ann1 = _make_annotation(start=1, end=1)  # idx 0:1 → ['A']/['K']
    ann2 = _make_annotation(start=4, end=5)  # idx 3:5 → ['G','A']/['V','H']
    n, protein = _make_nucleotides(['A', 'C', 'T', 'G', 'A'], ['K', 'M', 'R', 'V', 'H'], annotations=[ann1, ann2])
    cds = get_cds(n, protein)
    assert cds.alignment.reference == ['A', 'G', 'A']
    assert cds.alignment.query == ['K', 'V', 'H']


def test_get_cds_result_wraps_nucleotides_and_protein(mock_adjust_frame) -> None:
    n, protein = _make_nucleotides(['A', 'C', 'T'], ['K', 'M', 'R'])
    cds = get_cds(n, protein)
    assert cds.nucleotides is n
    assert cds.protein is protein


def test_get_cds_gap_in_alignment_maps_positions_correctly(mock_adjust_frame) -> None:
    # positions: [None, 1, 2, 3] — annotation covers positions 2–3 → indices 2:4
    ann = _make_annotation(start=2, end=3)
    n, protein = _make_nucleotides([GAP_SYMBOL, 'A', 'C', 'T'], ['X', 'K', 'M', 'R'], annotations=[ann])
    cds = get_cds(n, protein)
    assert cds.alignment.reference == ['C', 'T']
    assert cds.alignment.query == ['M', 'R']


def test_get_cds_calls_adjust_frame_with_cds(mock_adjust_frame) -> None:
    n, protein = _make_nucleotides(['A', 'C', 'T'], ['K', 'M', 'R'])
    cds = get_cds(n, protein)
    mock_adjust_frame.assert_called_once_with(cds)


# ---------------------------------------------------------------------------
# translate_codon tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'codon, expected',
    [
        (['N', 'T', 'G'], UNKNOWN_AA_SYMBOL),  # wildcard at start
        (['A', 'N', 'G'], UNKNOWN_AA_SYMBOL),  # wildcard in middle
        ([], UNKNOWN_AA_SYMBOL),  # empty codon
        (['A', 'T'], UNKNOWN_AA_SYMBOL),  # one base short
        (['A', 'T', 'G'], 'M'),  # ATG = Met
        (['T', 'T', 'T'], 'F'),  # TTT = Phe
        (['T', 'A', 'A'], STOP_CODON_SYMBOL),  # stop codon
        (['-', '-', '-'], '-'),  # gap codon (--- → -)
        (['R', 'T', 'G'], 'MV'),  # R=A|G → ATG=M, GTG=V → sorted 'MV'
        (['T', 'C', 'Y'], 'S'),  # Y=C|T → TCC=TCS=S, dedup → 'S'
        (['A', 'T', '-'], UNKNOWN_AA_SYMBOL),  # mixed gap → 'AT-' not in dict
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
