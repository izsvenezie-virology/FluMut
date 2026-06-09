from PySide6.QtWidgets import QComboBox, QLabel, QLineEdit

from flumut.flumutdb.models import Reference, Segment
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import BaseForm


class ReferenceForm(BaseForm):
    def __init__(self, parent=None, reference=None):
        self.reference = reference
        super().__init__(parent, 'Reference')

    def init_ui(self):
        super().init_ui()
        self.name_field = QLineEdit()
        self.segment_combo = QComboBox()

        segments = Segment.select()
        for segment in segments:
            self.segment_combo.addItem(segment.name, segment.id)

        self.form_layout.insertWidget(0, QLabel('Name:'))
        self.form_layout.insertWidget(1, self.name_field)
        self.form_layout.insertWidget(2, QLabel('Segment:'))
        self.form_layout.insertWidget(3, self.segment_combo)

        if self.reference:
            self.name_field.setText(self.reference.name)
            self.segment_combo.setCurrentText(self.reference.segment.name)

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(
                self, 'Name', 'Name cannot be empty.'
            )
            return False
        if self.segment_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(
                self, 'Segment', 'Please select a segment.'
            )
            return False
        return True

    def save_to_db(self) -> None:
        name = self.name_field.text().strip()
        segment_id = self.segment_combo.currentData()
        if self.reference:
            self.reference.name = name
            self.reference.segment_id = segment_id
            self.reference.save()
        else:
            Reference.create(name=name, segment_id=segment_id)

