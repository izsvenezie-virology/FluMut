from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.models import Reference


class ReferencesTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(['ID', 'Name', 'Segment'])
        references = Reference.select()
        self.table.setRowCount(len(references))
        for row, ref in enumerate(references):
            self.table.setItem(row, 0, self.__class__._create_item(str(ref.id)))
            self.table.setItem(row, 1, self.__class__._create_item(ref.name))
            self.table.setItem(row, 2, self.__class__._create_item(ref.segment.name))

    @staticmethod
    def _create_item(text):
        from PySide6.QtWidgets import QTableWidgetItem

        return QTableWidgetItem(text)
