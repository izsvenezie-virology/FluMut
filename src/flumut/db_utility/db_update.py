import sys
import subprocess
from urllib.request import urlretrieve
from flumut.db_utility.db_connection import DBConnection


def update() -> None:
    db = DBConnection()
    if db.is_pip_package:
        pip_update()
    else:
        file_update()


def pip_update() -> None:
    '''Update flumutdb Python package using pip.'''
    _ = subprocess.check_output([sys.executable, '-m', 'pip', 'install',
                                '--upgrade', 'flumutdb'], stderr=subprocess.PIPE)
    pass


def file_update() -> None:
    '''Download and replace the database file, the Python package is not updated.'''
    url = 'https://github.com/izsvenezie-virology/FluMutDB/releases/latest/download/flumut_db.sqlite'
    db_path = DBConnection().db_file
    _ = urlretrieve(url, db_path)
