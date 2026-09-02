import sys

import qdarktheme
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QApplication

from flumut.core.options import DatabaseOptions
from flumut.flumutdb import open_database
from flumut_db_editor.gui.main_window import MainWindow

FONT_FAMILIES = ['Inter', 'Segoe UI', 'SF Pro Text', 'Noto Sans', 'Cantarell', 'DejaVu Sans']
FONT_SIZE = 11


def ui_font() -> QFont:
    """The first of `FONT_FAMILIES` installed on this system, a couple of points above the platform default."""
    font = QFont()
    font.setFamilies(FONT_FAMILIES)
    font.setPointSize(FONT_SIZE)
    return font


def main():
    db = None
    if len(sys.argv) > 1:
        db = sys.argv[1]
    open_database(DatabaseOptions(path=db, read_only=False))

    app = QApplication(sys.argv)
    qdarktheme.setup_theme()
    app.setFont(ui_font())

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()
