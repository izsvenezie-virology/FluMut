from dataclasses import dataclass
from io import TextIOWrapper

from Bio import SeqIO

from flumut.core.logger import LOGGER


@dataclass
class FastaSequence:
    """A single sequence record parsed from a FASTA file.

    Attributes:
        header: Full sequence header.
        sequence: Nucleotide or amino acid sequence string.
        file: Path to the source FASTA file.
    """

    header: str
    sequence: str
    file: str


def read_fasta(fasta: TextIOWrapper) -> list[FastaSequence]:
    """Parse all sequence records from a FASTA file.

    Args:
        fasta: Open file handle for the FASTA file.

    Returns:
        A list of FastaSequence objects, one per record in the file.
    """
    LOGGER.debug(f'Reading FASTA file {fasta.name}')
    result = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        header = sanitize_header(seq.description)
        sequence = sanitize_sequence(str(seq.seq))
        result.append(FastaSequence(header, sequence, fasta.name))
    LOGGER.debug(f'Read {len(result)} sequences from {fasta.name}')
    return result


def sanitize_header(header: str) -> str:
    return ' '.join(header.split())


def sanitize_sequence(sequence: str) -> str:
    return ''.join(sequence.split()).upper()
