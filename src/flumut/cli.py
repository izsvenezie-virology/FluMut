import logging
import traceback
import click
from click import File

from flumut import __author__, __contact__, __version__, initialize_logging
from flumut.analysis import analyze
from flumut.db_utility.db_connection import DBConnection
from flumut.db_utility.db_update import update
from flumut.output import set_output_file
from flumut.sequence_utility.fasta_handler import set_header_pattern


def update_db(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    old_version = DBConnection().version_string
    update()
    new_version = DBConnection().version_string
    if old_version == new_version:
        print(f'Already using latest flumutdb version ({new_version})')
    else:
        print(f'Updated flumutdb to version {new_version}')
    ctx.exit()


def print_all_versions(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(f'flumut: {__version__}')
    major, minor, date = DBConnection().version
    print(f'flumutdb: {DBConnection().version_string}')
    ctx.exit()


def set_name_regex(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    set_header_pattern(value)


def set_output(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    set_output_file(param.name, value)


def set_dbfile(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    DBConnection().db_file = value


def set_verbosity(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    initialize_logging(value)


def print_errors(error: Exception) -> None:
    if logging.root.level == logging.DEBUG:
        traceback.print_exc()
    else:
        logging.critical(f'{type(error).__name__}: {error}')


@click.command()
# Help and versions
@click.help_option('-h', '--help')
@click.version_option(__version__, '-v', '--version', message=f'%(prog)s, version %(version)s, by {__author__} ({__contact__})')
# Database selection, must be eager since it must be parsed before update and all-versions
@click.option('-D', '--db-file', type=str, callback=set_dbfile, expose_value=False, is_eager=True, help='Set source database.')
# Options that exits from the workflow
@click.option('-V', '--all-versions', is_flag=True, callback=print_all_versions, expose_value=False, help='Prints all versions and exit.')
@click.option('--update', is_flag=True, callback=update_db, expose_value=False, help='Update the database to the latest version and exit.')
# Advanced options
@click.option('--allow-unmatching-headers', is_flag=True, default=False, help='Uses the whole header if this does not match the regular expression pattern.')
@click.option('-r', '--relaxed', is_flag=True, help='Report markers of which at least one mutation is found.')
@click.option('-n', '--name-regex', type=str, callback=set_name_regex, expose_value=False, default=r'(?P<sample>.+)_(?P<segment>.+)', show_default=True,
              help='Set regular expression to parse sequence name.')
@click.option('--loglevel', type=click.Choice(['dbg', 'inf', 'wrn', 'err'], case_sensitive=False), callback=set_verbosity, expose_value=False, default='wrn', show_default=True,
              help='Verbosity of the logging messages')
# Output files
@click.option('-m', '--markers-output', callback=set_output, expose_value=False, type=File('w', 'utf-8'), default=None, help='TSV markers output file.')
@click.option('-M', '--mutations-output', callback=set_output, expose_value=False, type=File('w', 'utf-8'), default=None, help='TSV mutations output file.')
@click.option('-l', '--literature-output', callback=set_output, expose_value=False, type=File('w', 'utf-8'), default=None, help='TSV literature output file.')
@click.option('-x', '--excel-output', callback=set_output, expose_value=False, type=File('w', lazy=False), default=None, help='Excel complete report.')
# Input files
@click.argument('fasta-files', type=File('r'), nargs=-1)
def cli(fasta_files: File, relaxed: bool, allow_unmatching_headers: bool) -> None:
    try:
        analyze(fasta_files, relaxed, allow_unmatching_headers)
    except Exception as e:
        print_errors(e)


if __name__ == '__main__':
    cli()
