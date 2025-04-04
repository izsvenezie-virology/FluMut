from io import TextIOWrapper
import logging
import re
from typing import List, Optional, Tuple

from flumut.exceptions import MalformedFastaException
from flumut.sequence_utility.models import FastaSequence

_header_pattern: re.Pattern = None
'''RegEx used to parse Fasta headers.'''


def set_header_pattern(pattern: str) -> None:
    '''
    Set the pattern used for header parsing.

    :param `str` pattern: RegEx string pattern.
    '''
    logging.debug(f'Set "{pattern}" as header pattern.')
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
    logging.debug(f'Reading file {fasta.name}.')
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
            raise MalformedFastaException(fasta.name) from None
    logging.debug(f'Found {len(sequences)} sequences.')
    return sequences


def parse_header(header: str) -> Tuple[Optional[str], Optional[str]]:
    '''
    Parse the header with the header pattern.
    It searches for two groups in order to retrieve sample and segment from the header.
    The header pattern can be set with set_header_pattern function.

    :param `str` header: The FASTA header.
    :param `str` allow_unmatching_headers: If true, continues the execution in 
    :return `Optional[str]` sample: The sample name. Returns `None` if  not found.
    :return `Optional[str]`segment: The segment name. Returns `None` if not found.
    '''
    logging.debug(f'Parsing "{header}" with "{get_header_pattern()}"')
    sample = None
    segment = None

    match = _header_pattern.match(header)
    try:
        sample = match.groupdict().get('sample', match.group(1))
        segment = match.groupdict().get('segment', match.group(2))
    except IndexError as e:
        pass

    logging.debug(f'Sample:  "{sample}"')
    logging.debug(f'Segment: "{segment}"')
    return sample, segment
