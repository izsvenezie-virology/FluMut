import subprocess
import sys
from urllib.request import urlretrieve

from flumut.db_utility.db_connection import DBConnection


def update() -> None:
    '''Update FluMutDB.'''
    db = DBConnection()
    if db.is_pip_package:
        _pip_update()
    else:
        _file_update()


def _pip_update() -> None:
    '''Update flumutdb Python package using pip.'''
    _ = subprocess.check_output([sys.executable, '-m', 'pip', 'install',
                                '--upgrade', 'flumutdb'], stderr=subprocess.PIPE)
    pass


def _file_update() -> None:
    '''Download and replace the database file, the Python package is not updated.'''
    url = 'https://github.com/izsvenezie-virology/FluMutDB/releases/latest/download/flumut_db.sqlite'
    db_path = DBConnection().db_file
    _ = urlretrieve(url, db_path)
