import sys

from PySide6.QtWidgets import QApplication

from flumut.flumutdb.initializer import initialize
from flumut_db_editor.gui.main_window import MainWindow


def main():
    initialize()

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()
