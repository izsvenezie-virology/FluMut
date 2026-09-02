import logging
from io import TextIOWrapper

import click
from click import File

from flumut import __author__, __contact__, __version__
from flumut.core import logger
from flumut.core.logger import LEVELS, LOGGER
from flumut.core.workflows import whole_workflow


def set_verbosity(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    LOGGER.setLevel(LEVELS[value])


def print_errors(error: Exception) -> None:
    LOGGER.critical(f'{type(error).__name__}: {error}')
    if LOGGER.level == logging.DEBUG:
        raise error


@click.command()
# Eager options
@click.help_option('-h', '--help')
@click.version_option(__version__, '--version', message=f'%(prog)s, v.%(version)s, by {__author__} ({__contact__})')
@click.option(
    '--loglevel',
    type=click.Choice(logger.LEVELS.keys(), case_sensitive=False),
    callback=set_verbosity,
    expose_value=False,
    default='wrn',
    show_default=True,
    is_eager=True,
    help='Verbosity of the logging messages',
)
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
        whole_workflow(fasta_files, relaxed, name_regex, markers_output, mutations_output, literature_output, excel_output)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print_errors(e)


if __name__ == '__main__':
    cli()
