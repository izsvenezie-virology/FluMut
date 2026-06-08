from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.models import Marker


class MarkersTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(['ID', 'Name'])
        markers = Marker.select()
        self.table.setRowCount(len(markers))
        for row, marker in enumerate(markers):
            self.table.setItem(row, 0, self.__class__._create_item(str(marker.id)))
            self.table.setItem(row, 1, self.__class__._create_item(marker.name))

    @staticmethod
    def _create_item(text):
        from PySide6.QtWidgets import QTableWidgetItem

        return QTableWidgetItem(text)
