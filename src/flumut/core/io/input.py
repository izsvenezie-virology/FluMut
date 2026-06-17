from dataclasses import dataclass
from io import TextIOWrapper

from Bio import SeqIO

from flumut.core.logger import LOGGER


@dataclass
class FastaSequence:
    """A single sequence record parsed from a FASTA file.

    Attributes:
        header: Sequence identifier (the part after ``>`` up to the first space).
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
        result.append(FastaSequence(seq.name, str(seq.seq).upper(), fasta.name))
    LOGGER.debug(f'Read {len(result)} sequences from {fasta.name}')
    return result
