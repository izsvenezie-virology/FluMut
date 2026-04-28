import logging
from io import TextIOWrapper

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def collect_sequences(fastas: tuple[TextIOWrapper]) -> list[SeqRecord]:
    logging.info(f'Collecting sequences from {len(fastas)} Fasta files')
    sequences: list[SeqRecord] = []
    for fasta in fastas:
        sequences += read_fasta(fasta)
    logging.info(f'Collected {len(sequences)} sequences')
    return sequences


def read_fasta(fasta: TextIOWrapper) -> list[SeqRecord]:
    return list(SeqIO.parse(fasta, 'fasta'))
