import re
from io import TextIOWrapper

from flumut.core.analysis.checks import perform_checks
from flumut.core.analysis.models import Analysis
from flumut.core.analysis.preprocess import load_nucleotide_fasta
from flumut.core.analysis.scanner import analyse
from flumut.core.io.output import write_outputs
from flumut.core.logger import LOGGER


def whole_workflow(
    fasta_files: tuple[TextIOWrapper, ...],
    relaxed: bool,
    name_regex: str,
    markers_output: TextIOWrapper | None,
    mutations_output: TextIOWrapper | None,
    literature_output: TextIOWrapper | None,
    excel_output: TextIOWrapper | None,
) -> None:
    if not fasta_files:
        raise Exception("Missing argument 'FASTA_FILES'")
    analysis = Analysis()
    pattern = re.compile(name_regex)

    LOGGER.info(f'Reading {len(fasta_files)} FASTA file(s)...')
    for fasta in fasta_files:
        load_nucleotide_fasta(analysis, fasta, pattern)

    LOGGER.info(f'Checking {len(analysis.samples)} sample(s)...')
    perform_checks(analysis)

    LOGGER.info(f'Scanning {len(analysis.samples)} sample(s)...')
    analyse(analysis, relaxed)

    LOGGER.info('Writing outputs...')
    write_outputs(analysis, markers_output, mutations_output, literature_output, excel_output)
