from PySide6.QtWidgets import QLabel, QLineEdit

from flumut.flumutdb.models import Segment
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import BaseForm


class SegmentForm(BaseForm):
    def __init__(self, parent=None, segment=None):
        self.segment = segment
        super().__init__(parent, 'Segment')

    def init_ui(self):
        super().init_ui()
        self.name_field = QLineEdit()
        self.form_layout.insertWidget(0, QLabel('Name:'))
        self.form_layout.insertWidget(1, self.name_field)

        if self.segment:
            self.name_field.setText(self.segment.name)

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(
                self, 'Name', 'Name cannot be empty.'
            )
            return False
        if not self.segment:
            if Segment.select().where(Segment.name == name).exists():
                ValidationErrorDialog.show_validation_error(
                    self, 'Name', 'A segment with this name already exists.'
                )
                return False
        return True

    def save_to_db(self) -> None:
        name = self.name_field.text().strip()
        if self.segment:
            self.segment.name = name
            self.segment.save()
        else:
            Segment.create(name=name)

