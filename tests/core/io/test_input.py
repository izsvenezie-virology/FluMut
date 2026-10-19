from io import StringIO

import pytest

from flumut.core.io.input import read_fasta, sanitize_header, sanitize_sequence


def test_read_fasta_keeps_the_whole_description_as_header() -> None:
    """Parsed against the real SeqIO: the header is not truncated at the first space."""
    fasta = StringIO('>ABC DEF\nATCG')
    fasta.name = 'seqs.fasta'
    assert read_fasta(fasta)[0].header == 'ABC DEF'  # type: ignore[arg-type]


@pytest.mark.parametrize(
    'header, expected',
    [
        ('ABC def', 'ABC def'),
        ('ABC\tDEF', 'ABC DEF'),
    ],
    ids=['case_preserved', 'tab_becomes_space'],
)
def test_sanitize_header(header: str, expected: str) -> None:
    assert sanitize_header(header) == expected


@pytest.mark.parametrize(
    'sequence, expected',
    [
        ('ATG CTG', 'ATGCTG'),
        ('atgctg', 'ATGCTG'),
    ],
    ids=['spaces_removed', 'uppercased'],
)
def test_sanitize_sequence(sequence: str, expected: str) -> None:
    assert sanitize_sequence(sequence) == expected
