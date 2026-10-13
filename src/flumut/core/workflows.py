import os
import re

from flumut.core.analysis.checks import perform_checks
from flumut.core.analysis.models import Analysis
from flumut.core.analysis.preprocess import load_nucleotide_fasta
from flumut.core.analysis.scanner import analyse
from flumut.core.exceptions import (
    InvalidNameRegexError,
    MissingFastaFilesError,
    NoOutputRequestedError,
    UnwritableOutputError,
)
from flumut.core.io.output import write_outputs
from flumut.core.logger import LOGGER
from flumut.core.options import FluMutOptions
from flumut.flumutdb import close_database, open_database


def check_output_path(path: str) -> None:
    """Make sure a file can be created at ``path``, before any work is spent on producing it.

    Args:
        path: Where the output file is going to be written.

    Raises:
        UnwritableOutputError: If the directory of ``path`` does not exist or is not writable.
    """
    directory = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(directory):
        raise UnwritableOutputError(path, 'the directory does not exist')
    if not os.access(directory, os.W_OK):
        raise UnwritableOutputError(path, 'the directory is not writable')


def whole_workflow(options: FluMutOptions) -> None:
    """Run the complete analysis described by ``options``, from FASTA to reports.

    Args:
        options: Every option of the run, grouped by concern.

    Raises:
        MissingFastaFilesError: If no input FASTA file is given.
        NoOutputRequestedError: If no output file is requested, which would run
            the whole analysis and write nothing.
        UnwritableOutputError: If the Excel report cannot be written where requested.
        InvalidNameRegexError: If the name regex does not compile.
    """
    if not options.input.fasta_files:
        raise MissingFastaFilesError("Missing argument 'FASTA_FILES'")
    if not options.output.any_requested():
        raise NoOutputRequestedError('No output file requested, nothing would be written')
    if options.output.excel_output:
        check_output_path(options.output.excel_output)

    LOGGER.info('Initializing analysis...')
    analysis = Analysis()
    try:
        pattern = re.compile(options.input.name_regex)
    except re.error as e:
        raise InvalidNameRegexError(options.input.name_regex, e) from e

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
