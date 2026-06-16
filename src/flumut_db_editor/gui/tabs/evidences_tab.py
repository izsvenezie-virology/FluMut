from PySide6.QtWidgets import QTableWidgetItem

from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.evidence_form import EvidenceForm
from flumut_db_editor.models import Evidence


class EvidencesTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.new_button.clicked.connect(self.handle_new)
        self.load_data()

    def load_data(self):
        header = ['Marker', 'Paper', 'Effect', 'Subtype', 'Host', 'Notes']
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        evidences: list[Evidence] = Evidence.select()
        self.table.setRowCount(len(evidences))
        for row, evidence in enumerate(evidences):
            self.table.setItem(row, 0, QTableWidgetItem(evidence.marker.name or ''))
            self.table.setItem(row, 1, QTableWidgetItem(evidence.paper.short_name))
            self.table.setItem(row, 2, QTableWidgetItem(evidence.effect.name))
            self.table.setItem(row, 3, QTableWidgetItem(evidence.subtype.name))
            self.table.setItem(row, 4, QTableWidgetItem(evidence.host.name if evidence.host else ''))
            self.table.setItem(row, 5, QTableWidgetItem(evidence.notes or ''))
        self.table.resizeColumnsToContents()

    def handle_new(self):
        form = EvidenceForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Evidence created successfully.')
            self.load_data()

    def handle_edit(self):
        row = self.get_selected_item()
        if row is None:
            return
        evidences = list(Evidence.select())
        evidence = evidences[row]
        form = EvidenceForm(self, evidence)
        if form.exec():
            SuccessNotification.show_success(self, 'Evidence updated successfully.')
            self.load_data()

    def handle_delete(self):
        row = self.get_selected_item()
        if row is None:
            return
        evidences = list(Evidence.select())
        evidence = evidences[row]
        if delete_with_confirmation(self, 'Evidence', evidence.id, DatabaseOperations.delete_evidence):
            self.load_data()
