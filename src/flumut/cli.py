import logging
from collections.abc import Callable
from typing import NoReturn

import click
from click import File

from flumut import __author__, __contact__, __version__
from flumut.core import logger
from flumut.core.logger import LEVELS, LOGGER
from flumut.core.options import DEFAULT_NAME_REGEX, DatabaseOptions, FluMutOptions
from flumut.core.workflows import whole_workflow
from flumut.flumutdb import initialize
from flumut.flumutdb.models import DbVersion


def set_verbosity(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    LOGGER.setLevel(LEVELS[value])


def print_errors(error: Exception) -> NoReturn:
    LOGGER.critical(f'{type(error).__name__}: {error}')
    if LOGGER.level == logging.DEBUG:
        raise error
    raise SystemExit(1)


def print_all_versions(database: DatabaseOptions) -> None:
    """Print the version of FluMut and of the database it would analyse with."""
    initialize(database.path, database.read_only)
    click.echo(f'FluMut v.{__version__}; FluMutDB v.{DbVersion.get_or_none()}')


def option_group(*options: Callable) -> Callable:
    """Bundle Click options into one reusable decorator.

    Options are listed in the order they should appear in ``--help``, matching
    the order they would have as stacked decorators.
    """

    def decorator(func: Callable) -> Callable:
        for option in reversed(options):
            func = option(func)
        return func

    return decorator


general_options = option_group(
    click.help_option('-h', '--help'),
    click.version_option(__version__, '--version', message=f'%(prog)s, v.%(version)s, by {__author__} ({__contact__})'),
    click.option('--all-versions', is_flag=True, help='Print the FluMut and FluMutDB versions and exit.'),
    click.option(
        '--loglevel',
        type=click.Choice(logger.LEVELS.keys(), case_sensitive=False),
        callback=set_verbosity,
        expose_value=False,
        default='wrn',
        show_default=True,
        is_eager=True,
        help='Verbosity of the logging messages',
    ),
)

database_options = option_group(
    click.option('-D', '--db-file', 'path', type=str, default=None, help='Set source database.'),
)

analysis_options = option_group(
    click.option('-r', '--relaxed', is_flag=True, help='Report markers of which at least one mutation is found.'),
)

output_options = option_group(
    click.option('-m', '--markers-output', type=File('w', 'utf-8'), default=None, help='TSV markers output file.'),
    click.option('-M', '--mutations-output', type=File('w', 'utf-8'), default=None, help='TSV mutations output file.'),
    click.option('-l', '--literature-output', type=File('w', 'utf-8'), default=None, help='TSV literature output file.'),
    click.option('-x', '--excel-output', type=File('w', lazy=False), default=None, help='Excel complete report.'),
)

input_options = option_group(
    click.option(
        '-n',
        '--name-regex',
        type=str,
        default=DEFAULT_NAME_REGEX,
        show_default=True,
        help='Set regular expression to parse sequence name.',
    ),
    click.argument('fasta-files', type=File('r'), nargs=-1),
)


@click.command()
@general_options
@database_options
@analysis_options
@output_options
@input_options
def cli(all_versions: bool, **options) -> None:
    # Every option above must be named after a field of one option group, or from_flat rejects it.
    run_options = FluMutOptions.from_flat(options)
    try:
        if all_versions:
            print_all_versions(run_options.database)
            return
        whole_workflow(run_options)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print_errors(e)


if __name__ == '__main__':
    cli()
