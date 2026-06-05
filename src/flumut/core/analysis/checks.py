from typing import TYPE_CHECKING

from flumut.core.analysis.models import DuplicationCheck, EnlongationCheck, FrameshiftCheck, TruncationCheck
from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL

if TYPE_CHECKING:
    from flumut.core.analysis.models import Analysis, Sample


def perform_checks(analysis: 'Analysis') -> None:
    for sample in analysis.samples.values():
        check_duplications(sample)
        check_frameshifts(sample)
        check_truncation(sample)
        check_enlongation(sample)


def check_truncation(sample: 'Sample') -> None:
    for alignment in sample.alignments:
        for pos in range(len(alignment.alignment.query) - 1):
            if STOP_CODON_SYMBOL in alignment.alignment.query[pos]:
                sample.checks.append(TruncationCheck(sample.id, alignment.protein.name, pos + 1))


def check_enlongation(sample: 'Sample') -> None:
    for alignment in sample.alignments:
        if alignment.alignment.reference[-1] != STOP_CODON_SYMBOL:
            continue
        last_aa = alignment.alignment.query[-1]
        for aa in last_aa:
            if aa not in (STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL, GAP_SYMBOL):
                sample.checks.append(EnlongationCheck(sample.id, alignment.protein.name))
                break


def check_frameshifts(sample: 'Sample') -> None:
    for alignment in sample.alignments:
        if alignment.cds:
            for frameshift in alignment.cds.frameshifts:
                sample.checks.append(FrameshiftCheck(sample.id, alignment.protein.name, frameshift[0], frameshift[1]))


def check_duplications(sample: 'Sample') -> None:
    seen = set()
    for alignment in sample.alignments:
        if alignment.protein.name in seen:
            sample.checks.append(DuplicationCheck(sample.id, alignment.protein.name))
        else:
            seen.add(alignment.protein.name)
