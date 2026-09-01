from io import StringIO
from unittest.mock import MagicMock, patch

import pytest

from flumut.core.io.input import read_fasta, sanitize_header, sanitize_sequence


@pytest.fixture
def mock_seqio():
    with patch('flumut.core.io.input.SeqIO') as mock:
        yield mock


def _make_seq_record(description: str, sequence: str) -> MagicMock:
    record = MagicMock()
    record.description = description
    record.seq = sequence
    return record


@pytest.mark.parametrize(
    'records',
    [
        [],
        [('seq1', 'ACGT')],
        [('seq1', 'ACGT'), ('seq2', 'TGCA')],
    ],
    ids=['empty', 'single', 'multiple'],
)
def test_read_fasta(mock_seqio, records: list[tuple[str, str]]) -> None:
    """Every record becomes a FastaSequence carrying the source file name."""
    fasta = MagicMock()
    fasta.name = 'seqs.fasta'
    mock_seqio.parse.return_value = [_make_seq_record(name, seq) for name, seq in records]

    result = read_fasta(fasta)

    mock_seqio.parse.assert_called_once_with(fasta, 'fasta')
    assert [(r.header, r.sequence, r.file) for r in result] == [(name, seq, 'seqs.fasta') for name, seq in records]


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
