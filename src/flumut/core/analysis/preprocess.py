import logging
import re
from io import TextIOWrapper

from flumut.core.analysis.models import Analysis, Sample
from flumut.core.io.input import read_fasta
from flumut.core.nucleotides.aligner import get_best_alignment, select_candidate_references
from flumut.core.nucleotides.translator import translate


def parse_header(header: str, pattern: re.Pattern) -> tuple[str, str | None]:
    """
    Parse the header with the header pattern.
    It searches for two groups in order to retrieve sample and segment from the header.
    The header pattern can be set with set_header_pattern function.

    :param `str` header: The FASTA header.
    :return `Optional[str]` sample: The sample name. Returns `None` if  not found.
    :return `Optional[str]`segment: The segment name. Returns `None` if not found.
    """
    logging.debug(f'Parsing "{header}" with "{pattern.pattern}"')
    sample, segment = None, None
    match = pattern.match(header)

    if not match:
        logging.info(f'Cannot extract sample ID from "{header}". The whole header is used as sample ID.')
        return header, None

    try:
        sample = match.groupdict().get('sample', match.group(1))
        segment = match.groupdict().get('segment', match.group(2))
    except IndexError:
        pass

    logging.debug(f'Sample:  "{sample}" - Segment: "{segment}"')
    return sample or header, segment


def load_nucleotide_fasta(analysis: Analysis, fasta: TextIOWrapper, header_pattern: re.Pattern) -> None:
    for sequence in read_fasta(fasta):
        sample, candidate_hint = parse_header(sequence.header, header_pattern)
        candidates = select_candidate_references(candidate_hint)
        nt_alignment = get_best_alignment(sequence, candidates)

        if sample not in analysis.samples:
            analysis.samples[sample] = Sample(sample)

        for protein in nt_alignment.reference.segment.proteins:
            aa_alignment = translate(nt_alignment, protein)
            analysis.samples[sample].alignments.append(aa_alignment)
