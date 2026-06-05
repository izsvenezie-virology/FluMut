from dataclasses import dataclass
from io import TextIOWrapper

from Bio import SeqIO

from flumut.core.logger import LOGGER


@dataclass
class FastaSequence:
    header: str
    sequence: str
    file: str


def read_fasta(fasta: TextIOWrapper) -> list[FastaSequence]:
    LOGGER.info(f'Reading FASTA file {fasta.name}')
    result = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        result.append(FastaSequence(seq.name, str(seq.seq), fasta.name))
    LOGGER.debug(f'Read {len(result)} sequences from {fasta.name}')
    return result
