import re
from io import TextIOWrapper

from flumut.core.analysis.models import Analysis, Sample
from flumut.core.io.input import read_fasta
from flumut.core.logger import LOGGER
from flumut.core.nucleotides.aligner import get_best_alignment, select_candidate_references
from flumut.core.nucleotides.translator import get_cds, translate


def parse_header(header: str, pattern: re.Pattern) -> tuple[str, str | None]:
    """Extract sample ID and segment name from a FASTA header string.

    Applies the given regex pattern to the header. The pattern is expected to
    define named groups ``sample`` and ``segment``, or positional groups 1 and 2
    respectively. If the pattern does not match, the raw header is returned as
    the sample ID with no segment.

    Args:
        header: Raw FASTA header string (without the leading ``>``).
        pattern: Compiled regex pattern used to extract sample and segment.

    Returns:
        A tuple of ``(sample_id, segment_name)`` where ``segment_name`` is
        ``None`` when not found.
    """
    LOGGER.debug(f'Parsing "{header}" with "{pattern.pattern}"')
    sample, segment = None, None
    match = pattern.match(header)

    if not match:
        LOGGER.info(f'Cannot extract sample ID from "{header}". The whole header is used as sample ID.')
        return header, None

    try:
        sample = match.groupdict().get('sample', match.group(1))
        segment = match.groupdict().get('segment', match.group(2))
    except IndexError:
        pass

    LOGGER.debug(f'Sample: "{sample}" - Segment: "{segment}"')
    return sample or header, segment


def load_nucleotide_fasta(analysis: Analysis, fasta: TextIOWrapper, header_pattern: re.Pattern) -> None:
    """Read a nucleotide FASTA file and populate the analysis with aligned protein sequences.

    For each sequence in the FASTA file, parses the header to determine the sample
    ID and segment hint, selects candidate references, computes the best nucleotide
    alignment, and translates the CDS for every protein on the matched reference
    segment. Results are appended to the corresponding Sample in ``analysis.samples``.

    Args:
        analysis: The Analysis object to populate with the loaded sequences.
        fasta: Open file handle for the input FASTA file.
        header_pattern: Compiled regex used to extract sample and segment from headers.
    """
    for sequence in read_fasta(fasta):
        sample, candidate_hint = parse_header(sequence.header, header_pattern)
        candidates = select_candidate_references(candidate_hint)
        nt_alignment = get_best_alignment(sequence, candidates)

        if sample not in analysis.samples:
            analysis.samples[sample] = Sample(sample)

        for protein in nt_alignment.reference.segment.proteins:
            cds = get_cds(nt_alignment, protein)
            aa_alignment = translate(cds)
            analysis.samples[sample].alignments.append(aa_alignment)
