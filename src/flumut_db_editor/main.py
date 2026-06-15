import sys

import qdarktheme
from PySide6.QtWidgets import QApplication

from flumut.flumutdb.initializer import initialize
from flumut_db_editor.gui.main_window import MainWindow


def main():
    db = None
    if len(sys.argv) > 1:
        db = sys.argv[1]
    initialize(read_only=False, path=db)

    app = QApplication(sys.argv)
    qdarktheme.setup_theme()

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()
