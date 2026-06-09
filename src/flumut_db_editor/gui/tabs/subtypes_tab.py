from PySide6.QtWidgets import QTableWidgetItem

from flumut.flumutdb.models import Subtype
from flumut_db_editor.gui.base import BaseTableTab


class SubtypesTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        header = ['Name', 'Notes']
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        subtypes: list[Subtype] = Subtype.select()
        self.table.setRowCount(len(subtypes))
        for row, subtype in enumerate(subtypes):
            self.table.setItem(row, 0, QTableWidgetItem(subtype.name))
            self.table.setItem(row, 1, QTableWidgetItem(subtype.notes or ''))
        self.table.resizeColumnsToContents()
