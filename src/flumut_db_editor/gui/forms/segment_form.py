from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb.models import Protein, Segment
from flumut_db_editor.gui.dialogs import ConfirmationDialog, ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm
from flumut_db_editor.gui.forms.protein_form import ProteinForm


class SegmentForm(TransactionalForm):
    model = Segment

    def __init__(self, parent: QWidget | None = None, instance: Segment | None = None) -> None:
        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: Segment | None

    def init_ui(self) -> None:
        super().init_ui()
        self.setMinimumSize(480, 420)

        self.name_field = QLineEdit()
        self.form_layout.addWidget(QLabel('Name:'))
        self.form_layout.addWidget(self.name_field)
        self.form_layout.addWidget(QLabel('Proteins:'))

        proteins_row = QHBoxLayout()
        self.proteins_list = QListWidget()
        proteins_row.addWidget(self.proteins_list, stretch=1)

        btn_col = QVBoxLayout()
        self.add_protein_btn = QPushButton('Add')
        self.edit_protein_btn = QPushButton('Edit')
        self.remove_protein_btn = QPushButton('Remove')
        btn_col.addWidget(self.add_protein_btn)
        btn_col.addWidget(self.edit_protein_btn)
        btn_col.addWidget(self.remove_protein_btn)
        btn_col.addStretch()
        proteins_row.addLayout(btn_col)

        self.form_layout.addLayout(proteins_row)

        self.add_protein_btn.clicked.connect(self._add_protein)
        self.edit_protein_btn.clicked.connect(self._edit_protein)
        self.remove_protein_btn.clicked.connect(self._remove_protein)

        if self.instance:
            self.name_field.setText(self.instance.name)

        self._set_proteins_enabled(self.instance is not None)
        self._refresh_list()
        self.name_field.setFocus()

    def _set_proteins_enabled(self, enabled: bool) -> None:
        self.proteins_list.setEnabled(enabled)
        self.edit_protein_btn.setEnabled(enabled)
        self.remove_protein_btn.setEnabled(enabled)

    def _refresh_list(self) -> None:
        self.proteins_list.clear()
        if self.instance:
            for p in self.instance.proteins:
                item = QListWidgetItem(p.name)
                item.setData(Qt.ItemDataRole.UserRole, p)
                self.proteins_list.addItem(item)

    def _add_protein(self) -> None:
        if self.instance is None:
            if not self.validate():
                return
            self.instance = Segment.create(name=self.name_field.text().strip())
            self._set_proteins_enabled(True)

        protein_form = ProteinForm(self, None, self.instance)
        if protein_form.exec():
            self._refresh_list()

    def _edit_protein(self) -> None:
        if protein := self._get_selected_protein():
            protein_form = ProteinForm(self, protein)
            if protein_form.exec():
                self._refresh_list()

    def _remove_protein(self) -> None:
        if protein := self._get_selected_protein():
            if ConfirmationDialog.ask(self, 'Delete Protein', f'Are you sure you want to delete {protein}?'):
                protein.delete_instance()
                self._refresh_list()

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'Segment name cannot be empty.')
            return False
        existing = Segment.get_or_none(Segment.name == name)
        if existing is not None:
            if self.instance is None or existing.get_id() != self.instance.get_id():
                ValidationErrorDialog.show_validation_error(self, 'Name', 'A segment with this name already exists.')
                return False
        return True

    def field_values(self) -> dict:
        return {'name': self.name_field.text().strip()}

    def _get_selected_protein(self) -> Protein | None:
        if item := self.proteins_list.currentItem():
            return item.data(Qt.ItemDataRole.UserRole)
