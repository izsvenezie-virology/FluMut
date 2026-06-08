from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.models import Paper


class PapersTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(['ID', 'Short Name', 'Title', 'Year'])
        papers = Paper.select()
        self.table.setRowCount(len(papers))
        for row, paper in enumerate(papers):
            self.table.setItem(row, 0, self.__class__._create_item(str(paper.id)))
            self.table.setItem(row, 1, self.__class__._create_item(paper.short_name))
            self.table.setItem(row, 2, self.__class__._create_item(paper.title))
            self.table.setItem(row, 3, self.__class__._create_item(str(paper.year)))

    @staticmethod
    def _create_item(text):
        from PySide6.QtWidgets import QTableWidgetItem

        return QTableWidgetItem(text)
