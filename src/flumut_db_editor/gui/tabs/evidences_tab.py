from PySide6.QtWidgets import QAbstractItemView

from flumut.flumutdb import loader
from flumut.flumutdb.models import Evidence
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.evidence_form import EvidenceForm
from flumut_db_editor.gui.tabs.base import BaseTableTab


class EvidencesTab(BaseTableTab[Evidence]):
    """Lists every evidence. Rows are selected in bulk and edited together in a single form."""

    def __init__(self) -> None:
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        headers = ['Marker', 'Paper', 'Effect', 'Subtype', 'Host', 'Notes']
        self.table.setColumnCount(len(headers))
        self.table.setHorizontalHeaderLabels(headers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)

        self.new_btn.setText('New evidences')
        self.edit_btn.setText('Edit selected')

        self.refresh()

    def refresh(self, selected: Evidence | None = None) -> None:
        rows = {evidence: self.row_texts(evidence) for evidence in sorted(loader.get(Evidence), key=self.sort_key)}
        self.populate_table(rows, selected)
        self.table.resizeColumnsToContents()

    def row_texts(self, evidence: Evidence) -> list[str]:
        return [*self.sort_key(evidence), evidence.notes or '']

    def sort_key(self, evidence: Evidence) -> tuple[str, ...]:
        values = (evidence.marker, evidence.paper, evidence.effect, evidence.subtype, evidence.host)
        return tuple(str(value) if value else '' for value in values)

    def on_new_requested(self) -> None:
        self.edit(EvidenceForm(self))

    def on_edit_requested(self) -> None:
        if instances := self.get_selected_items():
            self.edit(EvidenceForm(self, instances))

    def on_delete_requested(self) -> None:
        instance = self.get_selected_item()
        if instance and DeleteForm.confirm_and_delete(instance, self):
            self.refresh()

    def edit(self, form: EvidenceForm) -> None:
        """Run `form` and, when it saves, show the table again around its first row."""
        if form.exec():
            self.refresh(next(iter(form.instances), None))
