from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QLineEdit, QWidget

from flumut.flumutdb.models import Protein, Segment
from flumut_db_editor.gui.forms.base import TransactionalForm


class ProteinForm(TransactionalForm[Protein]):
    model = Protein

    def __init__(self, parent: QWidget | None = None, instance: Protein | None = None, segment: Segment | None = None) -> None:
        super().__init__(parent, instance)
        self.force_segment = segment
        if instance:
            self.force_segment = instance.segment
        self.__init_ui()

    def __init_ui(self) -> None:
        self.resize(480, 120)

        self.name_field = QLineEdit()
        self.segment_combo = QComboBox()

        segments: list[Segment] = sorted(Segment.select())
        for segment in segments:
            self.segment_combo.addItem(segment.name, segment)
        if self.force_segment:
            self.segment_combo.setCurrentText(self.force_segment.name)
            self.segment_combo.setEnabled(False)

        if self.instance:
            self.name_field.setText(self.instance.name)

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

    def field_values(self) -> dict:
        return {
            'name': self.name,
            'segment': self.segment,
        }
