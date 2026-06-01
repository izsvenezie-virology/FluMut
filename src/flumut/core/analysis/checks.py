import logging
from typing import TYPE_CHECKING

from flumut.core.globals import GAP_SYMBOL, STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL

if TYPE_CHECKING:
    from flumut.core.analysis.models import Analysis


def perform_checks(analysis: 'Analysis') -> None:
    check_truncation(analysis)
    check_enlongation(analysis)
    check_frameshifts(analysis)
    check_duplications(analysis)


def check_truncation(analysis: 'Analysis') -> None:
    for sample in analysis.samples.values():
        for alignment in sample.alignments:
            for pos in range(len(alignment.alignment.query) - 1):
                if STOP_CODON_SYMBOL in alignment.alignment.query[pos]:
                    sample.checks.append(TruncationCheck(sample.id, alignment.protein.name, pos + 1))


def check_enlongation(analysis: 'Analysis') -> None:
    for sample in analysis.samples.values():
        for alignment in sample.alignments:
            last_aa = alignment.alignment.query[-1]
            if last_aa not in (STOP_CODON_SYMBOL, UNKNOWN_AA_SYMBOL, GAP_SYMBOL):
                sample.checks.append(EnlongationCheck(sample.id, alignment.protein.name))


def check_frameshifts(analysis: 'Analysis') -> None:
    for sample in analysis.samples.values():
        for alignment in sample.alignments:
            if alignment.cds:
                for frameshift in alignment.cds.frameshifts:
                    sample.checks.append(FrameshiftCheck(sample.id, alignment.protein.name, frameshift[0], frameshift[1]))


def check_duplications(analysis: 'Analysis') -> None:
    for sample in analysis.samples.values():
        seen = set()
        for alignment in sample.alignments:
            if alignment.protein.name in seen:
                sample.checks.append(DuplicationCheck(sample.id, alignment.protein.name))
            else:
                seen.add(alignment.protein.name)


class Check:
    def __init__(self, message: str):
        self.message = message
        logging.warning(self.message)


class TruncationCheck(Check):
    def __init__(self, sample_id: str, protein_name: str, position: int):
        message = f'Premature stop codon detected in sample "{sample_id}" for protein "{protein_name}" at {position}'
        super().__init__(message)


class EnlongationCheck(Check):
    def __init__(self, sample_id: str, protein_name: str):
        message = f'Enlongation detected in sample "{sample_id}" for protein "{protein_name}"'
        super().__init__(message)


class FrameshiftCheck(Check):
    def __init__(self, sample_id: str, protein_name: str, start: int, end: int | None):
        message = f'Frameshift detected in sample "{sample_id}" for protein "{protein_name}" from position {start} to {end or "end"}'
        super().__init__(message)


class DuplicationCheck(Check):
    def __init__(self, sample_id: str, protein_name: str):
        message = f'Duplicate protein "{protein_name}" detected in sample "{sample_id}"'
        super().__init__(message)
