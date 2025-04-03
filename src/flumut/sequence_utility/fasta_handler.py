from io import TextIOWrapper
import re
from typing import List, Tuple

from flumut.Exceptions import MalformedFastaException
from flumut.sequence_utility.models import FastaSequence

_header_pattern: re.Pattern = None
'''RegEx used to parse Fasta headers.'''


def set_header_pattern(pattern: str) -> None:
    '''
    Set the pattern used for header parsing.

    :param `str` pattern: RegEx string pattern.
    '''
    global _header_pattern
    _header_pattern = re.compile(pattern)


def get_header_pattern() -> str:
    '''
    Return the pattern used for header parsing.

    :return `str`: The RegEx string pattern.
    '''
    return _header_pattern.pattern


def read_fasta(fasta: TextIOWrapper) -> List[FastaSequence]:
    '''
    Read all the sequences in the Fasta file.

    :param `TextIOWrapper` fasta: The opened Fasta file.
    :return `List[FastaSequence]`: All the sequences found in the Fasta file.
    '''
    sequences = []
    for line in fasta:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            sequence = FastaSequence(fasta.name, line[1:])
            sequences.append(sequence)
            continue
        try:
            sequence.sequence += line.upper()
        except UnboundLocalError:
            raise MalformedFastaException() from None
    return sequences


def parse_header(header: str) -> Tuple[str, str]:
    '''
    Parse the header with the header pattern.
    It searches for two groups in order to retrieve sample and segment from the header.
    The header pattern can be set with set_header_pattern function.

    :param `str` header: The FASTA header.
    :return `str` sample: The sample name.
    :return `str`segment: The segment name. 
    '''
    match = _header_pattern.match(header)
    sample = match.groupdict().get('sample', match.group(1))
    segment = match.groupdict().get('segment', match.group(2))
    return sample, segment
