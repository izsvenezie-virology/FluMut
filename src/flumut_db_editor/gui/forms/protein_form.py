from typing import TYPE_CHECKING

from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QLineEdit, QWidget

from flumut.flumutdb.models import Protein, Segment
from flumut_db_editor.gui.forms.base import TransactionalForm


class ProteinForm(TransactionalForm):
    model = Protein

    def __init__(self, parent: QWidget | None = None, instance: Protein | None = None, segment: Segment | None = None) -> None:
        self.force_segment = segment
        if instance:
            self.force_segment = instance.segment

        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: Protein | None

    def init_ui(self) -> None:
        super().init_ui()
        self.resize(480, 120)

        self.name_field = QLineEdit()
        self.segment_combo = QComboBox()

        segments: list[Segment] = sorted(Segment.select())
        for segment in segments:
            self.segment_combo.addItem(segment.name, segment)

        if self.instance:
            self.name_field.setText(self.instance.name)

        if self.force_segment:
            self.segment_combo.setCurrentText(self.force_segment.name)
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

    @property
    def name(self) -> str:
        return self.name_field.text().strip()

    @property
    def segment(self) -> Segment:
        return self.segment_combo.currentData()

    def validate(self) -> bool:
        if not self.check_unique_required('name', self.name, 'Name', self.name_field):
            return False
        return True

    def field_values(self) -> dict:
        return {
            'name': self.name,
            'segment': self.segment,
        }
