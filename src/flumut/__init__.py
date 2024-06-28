__version__ = '0.6.0'
__author__ = 'Edoardo Giussani'
__contact__ = 'egiussani@izsvenezie.it'

from io import TextIOWrapper
from typing import Dict


def versions() -> Dict[str, str]:
    """Collect versions of FluMut and FluMutDB packages.
    
    Returns:
    versions : `Dict[str, str]`
            Versions for FluMut and FluMutDB packages.
    """
    from flumut.DbReader import get_db_version
    major, minor, date = get_db_version()
    return {
        'FluMut': __version__,
        'FluMutDB': f'{major}.{minor}, released on {date}'
    }


def update() -> None:
    """Updates FluMutDB package using Pip.
    """
    import flumut.Updater 
    flumut.Updater.update()


def analyze(name_regex: str, fasta_file: TextIOWrapper, db_file: str, 
         markers_output: TextIOWrapper, mutations_output: TextIOWrapper, literature_output: TextIOWrapper, excel_output: str,
         relaxed: bool = False, skip_unmatch_names: bool = False, skip_unknown_segments: bool = False) -> None:
    """Runs the FluMut analysis.

    Args:
        name_regex: `str`
            The regex used to parse FASTA header
        fasta_file: `TextIOWrapper`
            The opened FASTA file to analyze
        db_file: `str|None`
            Path of the db file. If `None` the file from FluMutDB is used
        markers_output: `TextIOWrapper|None`
            The opened file where write Markers output
        mutations_output: `TextIOWrapper|None`
            The opened file where write Mutations output
        literature_output: `TextIOWrapper|None`
            The opened file where write Literature output
        excel_output: `str`
            The complete path to the Excel output file
        relaxed: `bool`
            When `True`, all markers where at least one mutation is detected
        skip_unmatch_names: `bool`
            When `True`, unmatching names do not raise exceptions
        skip_unknown_segments: `bool`
            When `True`, unknown segments do not raise exceptions
    
    """

    if name_regex is None:
        name_regex = r'(?P<sample>.+)_(?P<segment>.+)'
    if db_file is not None:
        from flumut.DbReader import set_db_file
        set_db_file(db_file)

    import flumut.flumut
    flumut.flumut.analyze(name_regex=name_regex, 
                          fasta_file=fasta_file,
                          db_file=db_file,
                          markers_output=markers_output,
                          mutations_output=mutations_output,
                          literature_output=literature_output,
                          excel_output=excel_output,
                          relaxed=relaxed,
                          skip_unmatch_names=skip_unmatch_names,
                          skip_unknown_segments=skip_unknown_segments)
