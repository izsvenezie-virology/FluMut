from typing import TYPE_CHECKING

import peewee
from PySide6.QtWidgets import (
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb.models import Protein, Segment
from flumut_db_editor.database_operations import DeletionValidator
from flumut_db_editor.gui.dialogs import ForeignKeyViolationDialog, ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class SegmentForm(TransactionalForm):
    model = Segment

    def __init__(self, parent: QWidget | None = None, instance: Segment | None = None) -> None:
        super().__init__(parent, instance)
        if self.instance:
            self._begin_transaction()
        if TYPE_CHECKING:
            self.instance: Segment

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

    def _set_proteins_enabled(self, enabled: bool) -> None:
        self.proteins_list.setEnabled(enabled)
        self.edit_protein_btn.setEnabled(enabled)
        self.remove_protein_btn.setEnabled(enabled)

    def _refresh_list(self) -> None:
        self.proteins_list.clear()
        if self.instance:
            for p in Protein.select().where(Protein.segment == self.instance):
                self.proteins_list.addItem(p.name)

    def _add_protein(self) -> None:
        if self.instance is None:
            segment_name = self.name_field.text().strip()
            if not segment_name:
                ValidationErrorDialog.show_validation_error(self, 'Segment Name', 'Please enter a segment name before adding proteins.')
                return
            if Segment.select().where(Segment.name == segment_name).exists():
                ValidationErrorDialog.show_validation_error(self, 'Segment Name', 'A segment with this name already exists.')
                return
            self._begin_transaction()
            self.instance = Segment.create(name=segment_name)
            self._set_proteins_enabled(True)

        name, ok = QInputDialog.getText(self, 'Add Protein', 'Protein name:')
        if not ok:
            return
        name = name.strip()
        if not name:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'Protein name cannot be empty.')
            return

        try:
            with self._db.savepoint():
                Protein.create(name=name, segment=self.instance)
        except peewee.IntegrityError:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'A protein with this name already exists.')
            return

        self._refresh_list()

    def _edit_protein(self) -> None:
        row = self.proteins_list.currentRow()
        if row < 0:
            return
        proteins: list[Protein] = list(Protein.select().where(Protein.segment == self.instance))
        protein = proteins[row]

        name, ok = QInputDialog.getText(self, 'Edit Protein', 'Protein name:', text=protein.name)
        if not ok:
            return
        name = name.strip()
        if not name:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'Protein name cannot be empty.')
            return

        try:
            with self._db.savepoint():
                protein.name = name
                protein.save()
        except peewee.IntegrityError:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'A protein with this name already exists.')
            return

        self._refresh_list()
        self.proteins_list.setCurrentRow(row)

    def _remove_protein(self) -> None:
        row = self.proteins_list.currentRow()
        if row < 0:
            return
        proteins: list[Protein] = list(Protein.select().where(Protein.segment == self.instance))
        protein = proteins[row]

        can_delete, violations = DeletionValidator.can_delete_protein(protein.get_id())
        if not can_delete:
            ForeignKeyViolationDialog.show_violation(self, 'Protein', protein.name, violations)
            return

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
