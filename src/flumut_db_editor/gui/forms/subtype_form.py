from PySide6.QtWidgets import QLabel, QLineEdit

from flumut.flumutdb.models import Subtype
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import BaseForm


class SubtypeForm(BaseForm):
    def __init__(self, parent=None, subtype=None):
        self.subtype = subtype
        super().__init__(parent, 'Subtype')

    def init_ui(self):
        super().init_ui()
        self.name_field = QLineEdit()
        self.form_layout.insertWidget(0, QLabel('Name:'))
        self.form_layout.insertWidget(1, self.name_field)

        if self.subtype:
            self.name_field.setText(self.subtype.name)

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'Name cannot be empty.')
            return False
        return True

    def save_to_db(self) -> None:
        name = self.name_field.text().strip()
        if self.subtype:
            self.subtype.name = name
            self.subtype.save()
        else:
            Subtype.create(name=name)
