from unittest.mock import MagicMock, patch

import pytest

from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL
from flumut.core.nucleotides.models import CDS, Alignment
from flumut.core.nucleotides.translator import get_cds, translate, translate_codon


@pytest.fixture
def mock_get_cds():
    with patch('flumut.core.nucleotides.translator.get_cds') as mock:
        yield mock


def _make_cds(reference_chars: list[str], query_chars: list[str]) -> CDS:
    cds = CDS(
        nucleotides=MagicMock(),
        protein=MagicMock(),
        alignment=Alignment(reference=reference_chars, query=query_chars),
    )
    return cds


# ---------------------------------------------------------------------------
# translate tests
# ---------------------------------------------------------------------------


def test_translate_produces_correct_amino_acid_sequences(mock_get_cds) -> None:
    # ATG → M, AAA → K  (reference)
    # GTG → V, CAT → H  (query)
    cds = _make_cds(
        reference_chars=['A', 'T', 'G', 'A', 'A', 'A'],
        query_chars=['G', 'T', 'G', 'C', 'A', 'T'],
    )
    mock_get_cds.return_value = cds
    result = translate(MagicMock(), MagicMock())
    assert result.alignment.reference == ['M', 'K']
    assert result.alignment.query == ['V', 'H']


def test_translate_calls_adjust_frame_on_cds(mock_get_cds) -> None:
    cds = _make_cds(['A', 'T', 'G'], ['G', 'T', 'G'])
    mock_get_cds.return_value = cds
    with patch('flumut.core.nucleotides.translator.adjust_frame') as mock_adjust:
        translate(MagicMock(), MagicMock())
        mock_adjust.assert_called_once_with(cds)


def test_translate_passes_arguments_to_get_cds(mock_get_cds) -> None:
    nucleotides = MagicMock()
    protein = MagicMock()
    mock_get_cds.return_value = _make_cds(['A', 'T', 'G'], ['G', 'T', 'G'])
    translate(nucleotides, protein)
    mock_get_cds.assert_called_once_with(nucleotides, protein)


def test_translate_result_stores_protein_reference_and_cds(mock_get_cds) -> None:
    nucleotides = MagicMock()
    protein = MagicMock()
    cds = _make_cds(['A', 'T', 'G'], ['G', 'T', 'G'])
    mock_get_cds.return_value = cds
    result = translate(nucleotides, protein)
    assert result.protein is protein
    assert result.reference is nucleotides.reference
    assert result.cds is cds


# ---------------------------------------------------------------------------
# get_cds helpers and tests
# ---------------------------------------------------------------------------


def _make_nucleotides(reference_chars: list[str], query_chars: list[str]) -> MagicMock:
    n = MagicMock()
    n.alignment = Alignment(reference=reference_chars, query=query_chars)
    return n


def _make_annotation(reference: MagicMock, start: int, end: int) -> MagicMock:
    ann = MagicMock()
    ann.reference = reference
    ann.start = start
    ann.end = end
    return ann


def test_get_cds_no_annotations_returns_empty_alignment() -> None:
    n = _make_nucleotides(['A', 'C', 'T'], ['K', 'M', 'R'])
    cds = get_cds(n, MagicMock(annotations=[]))
    assert cds.alignment.reference == []
    assert cds.alignment.query == []


def test_get_cds_matching_annotation_extracts_correct_slice() -> None:
    n = _make_nucleotides(['A', 'C', 'T', 'G'], ['K', 'M', 'R', 'V'])
    # positions: [1, 2, 3, 4] — annotation covers positions 2–3 → indices 1:3
    ann = _make_annotation(reference=n.reference, start=2, end=3)
    cds = get_cds(n, MagicMock(annotations=[ann]))
    assert cds.alignment.reference == ['C', 'T']
    assert cds.alignment.query == ['M', 'R']


def test_get_cds_nonmatching_annotation_excluded() -> None:
    n = _make_nucleotides(['A', 'C', 'T'], ['K', 'M', 'R'])
    ann = _make_annotation(reference=MagicMock(), start=1, end=2)
    cds = get_cds(n, MagicMock(annotations=[ann]))
    assert cds.alignment.reference == []
    assert cds.alignment.query == []


def test_get_cds_multiple_annotations_concatenated() -> None:
    n = _make_nucleotides(['A', 'C', 'T', 'G', 'A'], ['K', 'M', 'R', 'V', 'H'])
    # positions: [1, 2, 3, 4, 5]
    ann1 = _make_annotation(reference=n.reference, start=1, end=1)  # idx 0:1 → ['A']/['K']
    ann2 = _make_annotation(reference=n.reference, start=4, end=5)  # idx 3:5 → ['G','A']/['V','H']
    cds = get_cds(n, MagicMock(annotations=[ann1, ann2]))
    assert cds.alignment.reference == ['A', 'G', 'A']
    assert cds.alignment.query == ['K', 'V', 'H']


def test_get_cds_result_wraps_nucleotides_and_protein() -> None:
    n = _make_nucleotides(['A', 'C', 'T'], ['K', 'M', 'R'])
    protein = MagicMock(annotations=[])
    cds = get_cds(n, protein)
    assert cds.nucleotides is n
    assert cds.protein is protein


def test_get_cds_gap_in_alignment_maps_positions_correctly() -> None:
    n = _make_nucleotides([GAP_SYMBOL, 'A', 'C', 'T'], ['X', 'K', 'M', 'R'])
    # positions: [None, 1, 2, 3] — annotation covers positions 2–3 → indices 2:4
    ann = _make_annotation(reference=n.reference, start=2, end=3)
    cds = get_cds(n, MagicMock(annotations=[ann]))
    assert cds.alignment.reference == ['C', 'T']
    assert cds.alignment.query == ['M', 'R']


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
    with pytest.raises(ValueError, match='Unknown nucleotide'):
        translate_codon(['X', 'T', 'G'])
