from difflib import SequenceMatcher
from io import TextIOWrapper
import re
import sys
from typing import Generator, Tuple

from Bio.Align import PairwiseAligner
from flumut.Exceptions import MalformedFastaException, UnmatchHeaderException


PRINT_ALIGNMENT = False


class ProteinSequence:
    def __init__(self, sequence, cds) -> None:
        self.sequence: FastaSequence = sequence
        pass


class FastaSequence:
    def __init__(self, header: str, sequence: str) -> None:
        self.header: str = header
        self.sequence: str = sequence

    def parse_header(self, pattern: re.Pattern) -> None:
        '''Get sample and segment information from FASTA header.'''
        match = pattern.match(self.header)
        try:
            self.sample = match.groupdict().get('sample', match.group(1))
            self.segment = match.groupdict().get('segment', match.group(2))
        except (IndexError, AttributeError):
            raise UnmatchHeaderException(self.header, pattern.pattern) from None

    def align(self, references):
        higher_similarity = -1
        aligner = PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1

        for ref_name, ref_seq in references.items():
            alignment = aligner.align(ref_seq, self.sequence)[0]
            similarity = SequenceMatcher(None, alignment[0], alignment[1], False).ratio()
            if similarity > higher_similarity:
                higher_similarity = similarity
                self.reference_name = ref_name
                self.reference_aligned = alignment[0]
                self.sequence_aligned = alignment[1]

    def clean(self):
        del self.sequence
        del self.reference_aligned
        del self.sequence_aligned


def read(fasta_file: TextIOWrapper) -> Generator[FastaSequence, None, None]:
    '''Create a Fasta reading a file in Fasta format'''
    header = None
    for raw_line in fasta_file:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield FastaSequence(header, ''.join(seq).upper())
            header = line[1:]
            seq = []
        else:
            try:
                seq.append(line)
            except UnboundLocalError:
                raise MalformedFastaException() from None
    if header is not None:
        yield header, ''.join(seq).upper()
