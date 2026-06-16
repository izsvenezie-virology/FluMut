from typing import TYPE_CHECKING

from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QLineEdit, QWidget

from flumut.flumutdb.models import Protein, Segment
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class ProteinForm(TransactionalForm):
    model = Protein

    def __init__(self, parent: QWidget | None = None, instance: Protein | None = None, segment: Segment | None = None) -> None:
        self.segment = segment
        if instance:
            self.segment = instance.segment

        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: Protein | None

    def init_ui(self) -> None:
        super().init_ui()
        self.name_field = QLineEdit()
        self.segment_combo = QComboBox()

        segments: list[Segment] = Segment.select()
        for segment in segments:
            self.segment_combo.addItem(segment.name, segment)

        if self.instance:
            self.name_field.setText(self.instance.name)

        if self.segment:
            self.segment_combo.setCurrentText(self.segment.name)
            self.segment_combo.setEnabled(False)

        name_row = QHBoxLayout()
        name_row.insertWidget(0, QLabel('Name:'))
        name_row.insertWidget(1, self.name_field)
        segment_row = QHBoxLayout()
        segment_row.insertWidget(0, QLabel('Segment:'))
        segment_row.insertWidget(1, self.segment_combo, 1)

        self.form_layout.addLayout(name_row)
        self.form_layout.addLayout(segment_row)
        self.form_layout.addStretch()

        self.name_field.setFocus()

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'Name cannot be empty.')
            self.name_field.setFocus()
            return False
        if existing := Protein.get_or_none(Protein.name == name):
            if self.instance is None or existing.get_id() != self.instance.get_id():
                ValidationErrorDialog.show_validation_error(self, 'Name', 'A protein with this name already exists.')
                self.name_field.setFocus()
                return False
        if self.segment_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Segment', 'Please select a segment.')
            self.name_field.setFocus()
            return False
        return True

    def field_values(self) -> dict:
        return {
            'name': self.name_field.text().strip(),
            'segment': self.segment_combo.currentData(),
        }
