import click
from click import File

import flumut
from flumut import __version__, __author__, __contact__


def update(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    flumut.update()
    ctx.exit()

def versions(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    for package, version in flumut.versions().items():
        print(f'{package}: {version}')
    ctx.exit()

@click.command()
@click.help_option('-h', '--help')
@click.version_option(__version__, '-v', '--version', message=f'%(prog)s, version %(version)s, by {__author__} ({__contact__})')
@click.option('-V', '--all-versions', is_flag=True, callback=versions, expose_value=False, is_eager=True, help='Prints all versions and exit.')
@click.option('--update', is_flag=True, callback=update, expose_value=False, is_eager=True, help='Update the database to the latest version and exit.')
@click.option('--skip-unmatch-names', is_flag=True, default=False, help='Skip sequences with name that does not match the regular expression pattern.')
@click.option('--skip-unknown-segments', is_flag=True, default=False, help='Skip sequences with segment not present in the database.')
@click.option('-r', '--relaxed', is_flag=True, help='Report markers of which at least one mutation is found.')
@click.option('-n', '--name-regex', type=str, default=r'(?P<sample>.+)_(?P<segment>.+)', show_default=True, help='Set regular expression to parse sequence name.')
@click.option('-D', '--db-file', type=str, default=None, help='Set source database.')
@click.option('-m', '--markers-output', type=File('w', 'utf-8'), default=None, help='TSV markers output file.')
@click.option('-M', '--mutations-output', type=File('w', 'utf-8'), default=None, help='TSV mutations output file.')
@click.option('-l', '--literature-output', type=File('w', 'utf-8'), default=None, help='TSV literature output file.')
@click.option('-x', '--excel-output', type=str, default=None, help='Excel complete report.')
@click.argument('fasta-file', type=File('r'))
def cli(name_regex: str, fasta_file: File, db_file: str, 
         markers_output: File, mutations_output: File, literature_output: File, excel_output: str,
         relaxed: bool, skip_unmatch_names: bool, skip_unknown_segments: bool) -> None:
    flumut.analyze(**locals())


if __name__ == '__main__':
    cli()
