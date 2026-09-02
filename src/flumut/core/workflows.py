import re

from flumut.core.analysis.checks import perform_checks
from flumut.core.analysis.models import Analysis
from flumut.core.analysis.preprocess import load_nucleotide_fasta
from flumut.core.analysis.scanner import analyse
from flumut.core.io.output import write_outputs
from flumut.core.logger import LOGGER
from flumut.core.options import FluMutOptions
from flumut.flumutdb import close_database, open_database


class MissingFastaFilesError(ValueError):
    """Raised when no FASTA files are provided to the workflow."""


class NoOutputRequestedError(ValueError):
    """Raised when the workflow is asked to run without a single output file."""


def whole_workflow(options: FluMutOptions) -> None:
    """Run the complete analysis described by ``options``, from FASTA to reports.

    Args:
        options: Every option of the run, grouped by concern.

    Raises:
        MissingFastaFilesError: If no input FASTA file is given.
        NoOutputRequestedError: If no output file is requested, which would run
            the whole analysis and write nothing.
    """
    if not options.input.fasta_files:
        raise MissingFastaFilesError("Missing argument 'FASTA_FILES'")
    if not options.output.any_requested():
        raise NoOutputRequestedError('No output file requested, nothing would be written')

    LOGGER.info('Initializing analysis...')
    analysis = Analysis()
    pattern = re.compile(options.input.name_regex)

    try:
        open_database(options.database)

        LOGGER.info(f'Reading {len(options.input.fasta_files)} FASTA file(s)...')
        for fasta in options.input.fasta_files:
            load_nucleotide_fasta(analysis, fasta, pattern)

        LOGGER.info(f'Checking {len(analysis.samples)} sample(s)...')
        perform_checks(analysis)

        LOGGER.info(f'Scanning {len(analysis.samples)} sample(s)...')
        analyse(analysis, options.analysis.relaxed)

        LOGGER.info('Writing outputs...')
        write_outputs(analysis, options.output)
    finally:
        close_database()
