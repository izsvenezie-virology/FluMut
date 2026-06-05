import logging
import re
from io import TextIOWrapper

import click
from click import File

from flumut import __author__, __contact__, __version__
from flumut.core import logger
from flumut.core.analysis.checks import perform_checks
from flumut.core.analysis.models import Analysis
from flumut.core.analysis.preprocess import load_nucleotide_fasta
from flumut.core.analysis.scanner import analyse
from flumut.core.io.output import write_outputs
from flumut.core.logger import LEVELS, LOGGER
from flumut.flumutdb import initialize
from flumut.flumutdb.models import DbVersion


def update_db(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    raise NotImplementedError()


def print_all_versions(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(f'FluMut v.{__version__}; FluMutDB v.{DbVersion.get_or_none()}')
    ctx.exit()


def set_dbfile(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    initialize(value, read_only=True)


def set_verbosity(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    LOGGER.setLevel(LEVELS[value])


def print_errors(error: Exception) -> None:
    LOGGER.critical(f'{type(error).__name__}: {error}')
    if LOGGER.level == logging.DEBUG:
        raise error


@click.command()
# Help and versions
@click.help_option('-h', '--help')
@click.version_option(__version__, '--version', message=f'%(prog)s, v.%(version)s, by {__author__} ({__contact__})')
# Database selection, must be eager since it must be parsed before update and all-versions
@click.option('-D', '--db-file', type=str, callback=set_dbfile, expose_value=False, is_eager=True, help='Set source database.')
# Options that exits from the workflow
@click.option('--all-versions', is_flag=True, callback=print_all_versions, expose_value=False, help='Prints all versions and exit.')
@click.option('--update', is_flag=True, callback=update_db, expose_value=False, help='Update the database to the latest version and exit.')
# Advanced options
@click.option('-r', '--relaxed', is_flag=True, help='Report markers of which at least one mutation is found.')
@click.option(
    '-n',
    '--name-regex',
    type=str,
    default=r'(?P<sample>.+)_(?P<segment>.+)',
    show_default=True,
    help='Set regular expression to parse sequence name.',
)
@click.option(
    '--loglevel',
    type=click.Choice(logger.LEVELS.keys(), case_sensitive=False),
    callback=set_verbosity,
    expose_value=False,
    default='wrn',
    show_default=True,
    help='Verbosity of the logging messages',
)
# Output files
@click.option('-m', '--markers-output', type=File('w', 'utf-8'), default=None, help='TSV markers output file.')
@click.option('-M', '--mutations-output', type=File('w', 'utf-8'), default=None, help='TSV mutations output file.')
@click.option('-l', '--literature-output', type=File('w', 'utf-8'), default=None, help='TSV literature output file.')
@click.option('-x', '--excel-output', type=File('w', lazy=False), default=None, help='Excel complete report.')
# Input files
@click.argument('fasta-files', type=File('r'), nargs=-1)
def cli(
    fasta_files: tuple[TextIOWrapper, ...],
    relaxed: bool,
    name_regex: str,
    markers_output: TextIOWrapper | None,
    mutations_output: TextIOWrapper | None,
    literature_output: TextIOWrapper | None,
    excel_output: TextIOWrapper | None,
) -> None:
    try:
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
    except Exception as e:
        print_errors(e)


if __name__ == '__main__':
    cli()
