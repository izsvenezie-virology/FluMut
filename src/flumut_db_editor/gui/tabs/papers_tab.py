from PySide6.QtWidgets import QTableWidgetItem

from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.models import Paper


class PapersTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        header = ['Short Name', 'Authors', 'Year', 'Title', 'Journal', 'DOI', 'URL', 'Notes']
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        papers: list[Paper] = Paper.select()
        self.table.setRowCount(len(papers))
        for row, paper in enumerate(papers):
            self.table.setItem(row, 0, QTableWidgetItem(paper.short_name))
            self.table.setItem(row, 1, QTableWidgetItem(paper.authors))
            self.table.setItem(row, 2, QTableWidgetItem(str(paper.year)))
            self.table.setItem(row, 3, QTableWidgetItem(paper.title))
            self.table.setItem(row, 4, QTableWidgetItem(paper.journal or ''))
            self.table.setItem(row, 5, QTableWidgetItem(paper.doi or ''))
            self.table.setItem(row, 6, QTableWidgetItem(paper.url or ''))
            self.table.setItem(row, 7, QTableWidgetItem(paper.notes or ''))
        self.table.resizeColumnsToContents()
