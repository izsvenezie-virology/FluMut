import sys
import subprocess
from urllib.request import urlretrieve
from flumut.db_utility.db_connection import DBConnection


def update():
    output = subprocess.check_output([sys.executable, '-m', 'pip', 'install', '-U', 'flumutdb'], stderr=subprocess.PIPE)


def update_db_file() -> None:
    url = 'https://github.com/izsvenezie-virology/FluMutDB/releases/latest/download/flumut_db.sqlite'
    db_path = DBConnection().db_file
    _ = urlretrieve(url, db_path)
