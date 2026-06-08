from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.models import Evidence


class EvidencesTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(['ID', 'Marker', 'Paper', 'Effect', 'Subtype'])
        evidences = Evidence.select()
        self.table.setRowCount(len(evidences))
        for row, evidence in enumerate(evidences):
            self.table.setItem(row, 0, self.__class__._create_item(str(evidence.id)))
            self.table.setItem(row, 1, self.__class__._create_item(evidence.marker.name if evidence.marker else ''))
            self.table.setItem(row, 2, self.__class__._create_item(evidence.paper.short_name if evidence.paper else ''))
            self.table.setItem(row, 3, self.__class__._create_item(evidence.effect.name if evidence.effect else ''))
            self.table.setItem(row, 4, self.__class__._create_item(evidence.subtype.name if evidence.subtype else ''))

    @staticmethod
    def _create_item(text):
        from PySide6.QtWidgets import QTableWidgetItem

        return QTableWidgetItem(text)
