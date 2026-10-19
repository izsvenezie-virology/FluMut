from collections.abc import Iterable
from functools import partial
from itertools import product

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QWidget,
)

from flumut.flumutdb import loader
from flumut.flumutdb.models import BaseModel, Effect, Evidence, Host, Marker, Paper, Subtype, Target
from flumut_db_editor.gui.forms.base import MultiInstanceForm, TransactionalForm
from flumut_db_editor.gui.forms.effect_form import EffectForm
from flumut_db_editor.gui.forms.host_form import HostForm
from flumut_db_editor.gui.forms.marker_form import MarkerForm
from flumut_db_editor.gui.forms.paper_form import PaperForm
from flumut_db_editor.gui.forms.subtype_form import SubtypeForm
from flumut_db_editor.gui.widgets import FilterableList

COLUMNS = ('Marker', 'Paper', 'Effect', 'Subtype', 'Host', 'Target', 'Notes')
NOTES_COLUMN = COLUMNS.index('Notes')


class EvidenceForm(MultiInstanceForm[Evidence]):
    model = Evidence

    def __init__(self, parent: QWidget | None = None, instances: Iterable[Evidence] | None = None) -> None:
        super().__init__(parent, instances or ())
        self.removed: list[Evidence] = []
        self.__init_ui()

    def __init_ui(self) -> None:
        self.resize(900, 600)

        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked | QAbstractItemView.EditTrigger.EditKeyPressed
        )
        self.table.horizontalHeader().setSectionResizeMode(NOTES_COLUMN, QHeaderView.ResizeMode.Stretch)

        self.remove_btn = QPushButton('Remove selected')
        remove_row = QHBoxLayout()
        remove_row.addWidget(self.remove_btn)
        remove_row.addStretch()

        self.marker_list = FilterableList[Marker]('Marker', self)
        self.paper_list = FilterableList[Paper]('Paper', self)
        self.effect_list = FilterableList[Effect]('Effect', self)
        self.subtype_list = FilterableList[Subtype]('Subtype', self)
        self.host_list = FilterableList[Host]('Host (optional)', self)
        self.target_list = FilterableList[Target]('Target (optional)', self)

        # Every list draws from a model and creates new items of it through its own form.
        self.sources: tuple[tuple[FilterableList, type[BaseModel], type[TransactionalForm]], ...] = (
            (self.marker_list, Marker, MarkerForm),
            (self.paper_list, Paper, PaperForm),
            (self.effect_list, Effect, EffectForm),
            (self.subtype_list, Subtype, SubtypeForm),
            (self.host_list, Host, HostForm),
            (self.target_list, Target, HostForm),
        )

        lists_row = QHBoxLayout()
        for widget, _, form_class in self.sources:
            lists_row.addWidget(widget)
            widget.new_btn.clicked.connect(partial(self.on_new_instance_requested, widget, form_class))

        self.add_btn = QPushButton('Add combinations')
        add_row = QHBoxLayout()
        add_row.addStretch()
        add_row.addWidget(self.add_btn)

        self.form_layout.addWidget(self.table)
        self.form_layout.addLayout(remove_row)
        self.form_layout.addLayout(lists_row)
        self.form_layout.addLayout(add_row)

        self.remove_btn.clicked.connect(self.on_remove_requested)
        self.add_btn.clicked.connect(self.on_add_requested)

        self.refresh_lists()
        self.refresh_table()

    def form_title(self) -> str:
        return 'Evidences'

    def refresh_lists(self) -> None:
        """Reload every list, keeping whatever the user had selected in it."""
        for widget, model, _ in self.sources:
            selected = widget.selected_instances()
            widget.set_instances(sorted(loader.get(model), key=str))
            widget.select_instances(selected)

    def on_new_instance_requested(self, widget: FilterableList, form_class: type[TransactionalForm]) -> None:
        """Create an item for `widget` with its own form, then select it so it joins the next combination."""
        form = form_class(self)
        if not form.exec():
            return
        widget.filter_field.clear()
        self.refresh_lists()
        widget.select_instances([form.instance])

    def refresh_table(self) -> None:
        self.table.clearContents()
        self.table.setRowCount(len(self.instances))
        for row, evidence in enumerate(self.instances):
            for col, text in enumerate(self.row_texts(evidence)):
                item = QTableWidgetItem(text)
                if col == NOTES_COLUMN:
                    item.setToolTip('Double-click to edit the notes of this evidence.')
                else:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if col == 0:
                    item.setData(Qt.ItemDataRole.UserRole, evidence)
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()

    def row_texts(self, evidence: Evidence) -> list[str]:
        return [str(value) if value else '' for value in self.combination(evidence)] + [evidence.notes or '']

    def combination(self, evidence: Evidence) -> tuple:
        return (evidence.marker, evidence.paper, evidence.effect, evidence.subtype, evidence.host, evidence.target)

    def on_add_requested(self) -> None:
        self.populate_instances()  # rebuilding the table redraws notes from the instances, so sync what was typed
        combinations = product(
            self.marker_list.selected_instances(),
            self.paper_list.selected_instances(),
            self.effect_list.selected_instances(),
            self.subtype_list.selected_instances(),
            self.host_list.selected_instances() or [None],
            self.target_list.selected_instances() or [None],
        )
        present = [self.combination(evidence) for evidence in self.instances]
        for marker, paper, effect, subtype, host, target in combinations:
            if (marker, paper, effect, subtype, host, target) in present:
                continue
            self.instances.append(Evidence(marker=marker, paper=paper, effect=effect, subtype=subtype, host=host, target=target))
        self.refresh_table()

    def on_remove_requested(self) -> None:
        self.populate_instances()  # rebuilding the table redraws notes from the instances, so sync what was typed
        for evidence in self.selected_evidences():
            self.instances.remove(evidence)
            if evidence.get_id() is not None:
                self.removed.append(evidence)
        self.refresh_table()

    def selected_evidences(self) -> list[Evidence]:
        rows = self.table.selectionModel().selectedRows()
        return [self.table.item(index.row(), 0).data(Qt.ItemDataRole.UserRole) for index in rows]  # pyright: ignore[reportOptionalMemberAccess]

    def populate_instances(self) -> None:
        """Copy the notes column back onto its evidence; every other column is read-only."""
        self.commit_open_editor()
        for row in range(self.table.rowCount()):
            evidence: Evidence = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)  # pyright: ignore[reportOptionalMemberAccess]
            notes = self.table.item(row, NOTES_COLUMN).text().strip()  # pyright: ignore[reportOptionalMemberAccess]
            evidence.notes = notes or None

    def commit_open_editor(self) -> None:
        """Close the cell editor, if one is open, so a note being typed is not lost."""
        if self.table.state() == QAbstractItemView.State.EditingState:
            self.table.setFocus()

    def instances_to_delete(self) -> Iterable[BaseModel]:
        return self.removed
