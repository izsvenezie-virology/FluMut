from PySide6.QtWidgets import QLabel, QLineEdit

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

    def get_data(self):
        return {'name': self.name_field.text()}
