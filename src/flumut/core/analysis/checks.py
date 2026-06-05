from typing import TYPE_CHECKING

from flumut.core.analysis.models import DuplicationCheck, EnlongationCheck, FrameshiftCheck, TruncationCheck
from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL
from flumut.core.logger import LOGGER

if TYPE_CHECKING:
    from flumut.core.analysis.models import Analysis, Sample


def perform_checks(analysis: 'Analysis') -> None:
    """Run all quality checks on every sample in the analysis.

    Args:
        analysis: The Analysis container holding all samples to validate.
    """
    for sample in analysis.samples.values():
        LOGGER.debug(f'Checking {sample.id}')
        check_duplications(sample)
        check_frameshifts(sample)
        check_truncation(sample)
        check_enlongation(sample)


def check_truncation(sample: 'Sample') -> None:
    """Detect premature stop codons in a sample's protein alignments.

    A TruncationCheck is appended for every stop codon symbol found at any
    position where the reference amino acid is not a stop codon.

    Args:
        sample: The sample whose protein alignments are inspected.
    """
    for alignment in sample.alignments:
        aln = alignment.alignment
        for pos in range(len(aln.query)):
            if STOP_CODON_SYMBOL in aln.query[pos] and aln.reference[pos] != STOP_CODON_SYMBOL:
                sample.checks.append(TruncationCheck(sample.id, alignment.protein.name, pos + 1))


def check_enlongation(sample: 'Sample') -> None:
    """Detect elongation beyond the expected stop codon.

    An EnlongationCheck is appended when the reference alignment ends with a
    stop codon but the query contains a valid non-stop amino
    acid in that final position.

    Args:
        sample: The sample whose protein alignments are inspected.
    """
    for alignment in sample.alignments:
        aln = alignment.alignment
        if aln.reference[-1] != STOP_CODON_SYMBOL:
            continue
        last_aa = aln.query[-1]
        for aa in last_aa:
            if aa not in (STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL, GAP_SYMBOL):
                sample.checks.append(EnlongationCheck(sample.id, alignment.protein.name))
                break


def check_frameshifts(sample: 'Sample') -> None:
    """Detect frameshift indels in a sample's CDS alignments.

    Relies on frameshift ranges pre-computed by ``adjust_frame`` during translation.
    A FrameshiftCheck is appended for every frameshift range found.

    Args:
        sample: The sample whose CDS objects are inspected.
    """
    for alignment in sample.alignments:
        if alignment.cds:
            for frameshift in alignment.cds.frameshifts:
                sample.checks.append(FrameshiftCheck(sample.id, alignment.protein.name, frameshift[0], frameshift[1]))


def check_duplications(sample: 'Sample') -> None:
    """Detect duplicate protein entries within a sample.

    A DuplicationCheck is appended for every protein name that appears more
    than once among the sample's alignments.

    Args:
        sample: The sample whose protein alignments are inspected.
    """
    seen = set()
    for alignment in sample.alignments:
        if alignment.protein.name in seen:
            sample.checks.append(DuplicationCheck(sample.id, alignment.protein.name))
        else:
            seen.add(alignment.protein.name)
