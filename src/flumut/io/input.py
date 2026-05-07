import logging
from io import TextIOWrapper

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_fasta(fasta: TextIOWrapper) -> list[SeqRecord]:
    logging.debug(f'Reading FASTA file {fasta.name}')
    return list(SeqIO.parse(fasta, 'fasta'))
