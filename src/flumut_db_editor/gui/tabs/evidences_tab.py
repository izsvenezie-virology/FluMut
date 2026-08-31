from collections import defaultdict
from collections.abc import Callable

from PySide6.QtWidgets import QAbstractItemView, QComboBox, QLabel

from flumut.flumutdb import loader
from flumut.flumutdb.models import BaseModel, Evidence
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.evidence_form import EvidenceForm
from flumut_db_editor.gui.tabs.base import BaseTreeTab

COLUMNS = ('Marker', 'Paper', 'Effect', 'Subtype', 'Host', 'Notes')

GROUPS: dict[str, Callable[[Evidence], BaseModel | None]] = {
    'Marker': lambda evidence: evidence.marker,
    'Paper': lambda evidence: evidence.paper,
    'Effect': lambda evidence: evidence.effect,
    'Nothing': lambda evidence: None,
}


class EvidencesTab(BaseTreeTab[Evidence | BaseModel]):
    """Lists every evidence under the marker, paper or effect it belongs to.

    Evidences are selected in bulk and edited together in a single form;
    selecting a group is a shorthand for selecting all the evidences under it.
    """

    def __init__(self) -> None:
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        self.tree.setColumnCount(len(COLUMNS))
        self.tree.setHeaderLabels(COLUMNS)
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)

        self.new_btn.setText('New evidences')
        self.edit_btn.setText('Edit selected')

        self.group_combo = QComboBox()
        self.group_combo.addItems(GROUPS)  # pyright: ignore[reportArgumentType]
        self.header.addWidget(QLabel('Group by:'))
        self.header.addWidget(self.group_combo)
        self.group_combo.currentTextChanged.connect(self.on_group_changed)

        self.refresh()

    def populate(self, selected: Evidence | BaseModel | None = None) -> None:
        self.tree.clear()

        for key, evidences in self.grouped_evidences().items():
            parent = None
            if key is not None:
                parent = self.create_item(key, [f'{key} ({len(evidences)})'])
                if key == selected:
                    self.set_selected_item(parent)

            for evidence in evidences:
                item = self.create_item(evidence, self.row_texts(evidence), parent)
                if evidence == selected:
                    self.set_selected_item(item)

        for column in range(len(COLUMNS)):
            self.tree.resizeColumnToContents(column)

    def grouped_evidences(self) -> dict[BaseModel | None, list[Evidence]]:
        key_of = GROUPS[self.group_combo.currentText()]
        groups: dict[BaseModel | None, list[Evidence]] = defaultdict(list)
        for evidence in sorted(loader.get(Evidence), key=self.sort_key):
            groups[key_of(evidence)].append(evidence)
        return dict(sorted(groups.items(), key=lambda group: str(group[0])))

    def row_texts(self, evidence: Evidence) -> list[str]:
        return [*self.sort_key(evidence), evidence.notes or '']

    def sort_key(self, evidence: Evidence) -> tuple[str, ...]:
        values = (evidence.marker, evidence.paper, evidence.effect, evidence.subtype, evidence.host)
        return tuple(str(value) if value else '' for value in values)

    def get_selected_evidences(self) -> list[Evidence]:
        evidences = []
        for item in self.tree.selectedItems():
            instance = self.get_data(item)
            if isinstance(instance, Evidence):
                evidences.append(instance)
            else:
                evidences.extend(self.get_data(item.child(row)) for row in range(item.childCount()))
        return list(dict.fromkeys(evidences))

    def on_group_changed(self) -> None:
        self.refresh(self.get_selected_instance())

    def on_new_requested(self) -> None:
        self.edit(EvidenceForm(self))

    def on_edit_requested(self) -> None:
        if evidences := self.get_selected_evidences():
            self.edit(EvidenceForm(self, evidences))

    def on_delete_requested(self) -> None:
        refresh = False
        for instance in self.get_selected_evidences():
            if isinstance(instance, Evidence) and DeleteForm.confirm_and_delete(instance, self):
                refresh = True
        if refresh:
            self.refresh()

    def edit(self, form: EvidenceForm) -> None:
        if form.exec():
            self.refresh(next(iter(form.instances), None))
