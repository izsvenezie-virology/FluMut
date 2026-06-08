from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.models import Mutation


class MutationsTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(['ID', 'Name', 'Type', 'Protein'])
        mutations = Mutation.select()
        self.table.setRowCount(len(mutations))
        for row, mutation in enumerate(mutations):
            self.table.setItem(row, 0, self.__class__._create_item(str(mutation.id)))
            self.table.setItem(row, 1, self.__class__._create_item(mutation.name))
            self.table.setItem(row, 2, self.__class__._create_item(mutation.type))
            self.table.setItem(row, 3, self.__class__._create_item(mutation.protein.name))

    @staticmethod
    def _create_item(text):
        from PySide6.QtWidgets import QTableWidgetItem

        return QTableWidgetItem(text)
