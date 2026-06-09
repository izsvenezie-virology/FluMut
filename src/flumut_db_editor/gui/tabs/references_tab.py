from PySide6.QtWidgets import QTableWidgetItem

from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.models import Reference


class ReferencesTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        header = ['Name', 'Segment', 'Source', 'Notes', 'Sequence']
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        references: list[Reference] = Reference.select()
        self.table.setRowCount(len(references))
        for row, ref in enumerate(references):
            self.table.setItem(row, 0, QTableWidgetItem(ref.name))
            self.table.setItem(row, 1, QTableWidgetItem(ref.segment.name))
            self.table.setItem(row, 2, QTableWidgetItem(ref.source))
            self.table.setItem(row, 3, QTableWidgetItem(ref.notes or ''))
            self.table.setItem(row, 4, QTableWidgetItem(ref.sequence))
        self.table.resizeColumnsToContents()
