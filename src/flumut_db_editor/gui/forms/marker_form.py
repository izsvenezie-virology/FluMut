from PySide6.QtWidgets import QLabel, QLineEdit

from flumut_db_editor.gui.forms.base import BaseForm


class MarkerForm(BaseForm):
    def __init__(self, parent=None, marker=None):
        self.marker = marker
        super().__init__(parent, 'Marker')

    def init_ui(self):
        super().init_ui()
        self.name_field = QLineEdit()
        self.form_layout.insertWidget(0, QLabel('Name:'))
        self.form_layout.insertWidget(1, self.name_field)

        if self.marker:
            self.name_field.setText(self.marker.name)

    def get_data(self):
        return {'name': self.name_field.text()}
