from io import StringIO
from unittest.mock import MagicMock, patch

import pytest

from flumut.core.io.input import FastaSequence, read_fasta

# ---------------------------------------------------------------------------
# FastaSequence tests
# ---------------------------------------------------------------------------


def test_fasta_sequence_stores_all_fields() -> None:
    fs = FastaSequence(header='seq1', sequence='ACGT', file='test.fasta')
    assert fs.header == 'seq1'
    assert fs.sequence == 'ACGT'
    assert fs.file == 'test.fasta'


# ---------------------------------------------------------------------------
# read_fasta tests
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_seqio():
    with patch('flumut.core.io.input.SeqIO') as mock:
        yield mock


def _make_seq_record(name: str, seq: str) -> MagicMock:
    record = MagicMock()
    record.description = name
    record.seq = seq
    return record


def test_read_fasta_empty_file_returns_empty_list(mock_seqio) -> None:
    fasta = MagicMock(name='empty.fasta')
    mock_seqio.parse.return_value = []
    assert read_fasta(fasta) == []


def test_read_fasta_single_sequence_returns_one_fasta_sequence(mock_seqio) -> None:
    fasta = MagicMock()
    fasta.name = 'seqs.fasta'
    mock_seqio.parse.return_value = [_make_seq_record('seq1', 'ACGT')]
    result = read_fasta(fasta)
    assert len(result) == 1
    assert result[0].header == 'seq1'
    assert result[0].sequence == 'ACGT'
    assert result[0].file == 'seqs.fasta'


def test_read_fasta_multiple_sequences_all_returned(mock_seqio) -> None:
    fasta = MagicMock()
    fasta.name = 'seqs.fasta'
    mock_seqio.parse.return_value = [
        _make_seq_record('seq1', 'ACGT'),
        _make_seq_record('seq2', 'TGCA'),
    ]
    result = read_fasta(fasta)
    assert len(result) == 2
    assert result[0].header == 'seq1'
    assert result[1].header == 'seq2'


def test_read_fasta_space_in_header_keeps_whole_header() -> None:
    fasta = StringIO('>ABC DEF\nATCG')
    fasta.name = 'seqs.fasta'
    result = read_fasta(fasta)  # type: ignore
    assert result[0].header == 'ABC DEF'


def test_read_fasta_passes_file_to_seqio(mock_seqio) -> None:
    fasta = MagicMock()
    fasta.name = 'seqs.fasta'
    mock_seqio.parse.return_value = []
    read_fasta(fasta)
    mock_seqio.parse.assert_called_once_with(fasta, 'fasta')
