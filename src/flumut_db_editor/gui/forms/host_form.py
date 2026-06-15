from PySide6.QtWidgets import QLabel, QLineEdit

from flumut.flumutdb.models import Host
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class HostForm(TransactionalForm):
    model = Host

    def __init__(self, parent=None, instance=None):
        super().__init__(parent, 'Host', instance)

    def init_ui(self):
        super().init_ui()
        self.name_field = QLineEdit()
        self.form_layout.insertWidget(0, QLabel('Name:'))
        self.form_layout.insertWidget(1, self.name_field)

        if self.instance:
            self.name_field.setText(self.instance.name)

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'Name cannot be empty.')
            return False
        return True

    def field_values(self) -> dict:
        return {'name': self.name_field.text().strip()}
