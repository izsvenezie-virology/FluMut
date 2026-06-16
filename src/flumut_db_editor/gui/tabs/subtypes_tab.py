from PySide6.QtWidgets import QTableWidgetItem

from flumut.flumutdb.models import Subtype
from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.subtype_form import SubtypeForm


class SubtypesTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.new_button.clicked.connect(self.handle_new)
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

    def handle_new(self):
        form = SubtypeForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Subtype created successfully.')
            self.load_data()

    def handle_edit(self):
        row = self.get_selected_item()
        if row is None:
            return
        subtypes = list(Subtype.select())
        subtype = subtypes[row]
        form = SubtypeForm(self, subtype)
        if form.exec():
            SuccessNotification.show_success(self, 'Subtype updated successfully.')
            self.load_data()

    def handle_delete(self):
        row = self.get_selected_item()
        if row is None:
            return
        subtypes = list(Subtype.select())
        subtype = subtypes[row]
        if delete_with_confirmation(self, 'Subtype', subtype.id, DatabaseOperations.delete_subtype):
            self.load_data()
