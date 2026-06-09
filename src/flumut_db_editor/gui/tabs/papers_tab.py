from PySide6.QtWidgets import QTableWidgetItem

from flumut.flumutdb.models import Paper
from flumut_db_editor.database_operations import DatabaseOperations
from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.gui.crud_mixin import delete_with_confirmation
from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.paper_form import PaperForm


class PapersTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.new_button.clicked.connect(self.handle_new)
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

    def handle_new(self):
        form = PaperForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Paper created successfully.')
            self.load_data()

    def handle_edit(self):
        row = self.get_selected_item()
        if row is None:
            return
        papers = list(Paper.select())
        paper = papers[row]
        form = PaperForm(self, paper)
        if form.exec():
            SuccessNotification.show_success(self, 'Paper updated successfully.')
            self.load_data()

    def handle_delete(self):
        row = self.get_selected_item()
        if row is None:
            return
        papers = list(Paper.select())
        paper = papers[row]
        if delete_with_confirmation(
            self, 'Paper', paper.id, DatabaseOperations.delete_paper
        ):
            self.load_data()
