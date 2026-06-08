from PySide6.QtWidgets import QComboBox, QLabel, QLineEdit

from flumut_db_editor.gui.forms.base import BaseForm
from flumut_db_editor.models import Segment


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

    def get_data(self):
        return {
            'name': self.name_field.text(),
            'segment_id': self.segment_combo.currentData(),
        }
