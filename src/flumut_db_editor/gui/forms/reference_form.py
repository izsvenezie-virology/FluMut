from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QTextEdit,
    QWidget,
)

from flumut.core.io.input import sanitize_sequence
from flumut.flumutdb.models import Reference, Segment
from flumut_db_editor.gui.forms.base import TransactionalForm


class ReferenceForm(TransactionalForm[Reference]):
    model = Reference

    def __init__(self, parent: QWidget | None = None, instance: Reference | None = None, segment: Segment | None = None) -> None:
        self.force_segment = segment
        if instance:
            self.force_segment = instance.segment

        super().__init__(parent, instance)

    def _init_ui(self):
        super()._init_ui()
        self.resize(1000, 500)

        self.name_field = QLineEdit()
        self.segment_combo = QComboBox()
        self.source_field = QLineEdit()
        self.sequence_txt = QTextEdit()

        segments: list[Segment] = sorted(Segment.select())
        for segment in segments:
            self.segment_combo.addItem(segment.name, segment)
        if self.force_segment:
            self.segment_combo.setCurrentText(self.force_segment.name)
            self.segment_combo.setEnabled(False)

        if self.instance:
            self.name_field.setText(self.instance.name)
            self.source_field.setText(self.instance.source)
            self.sequence_txt.setPlainText(self.instance.sequence)

        name_row = QHBoxLayout()
        name_row.insertWidget(0, QLabel('Name:'))
        name_row.insertWidget(1, self.name_field)
        segment_row = QHBoxLayout()
        segment_row.insertWidget(0, QLabel('Segment:'))
        segment_row.insertWidget(1, self.segment_combo, 1)
        source_row = QHBoxLayout()
        source_row.insertWidget(0, QLabel('Source:'))
        source_row.insertWidget(1, self.source_field)

        self.form_layout.addLayout(name_row)
        self.form_layout.addLayout(segment_row)
        self.form_layout.addLayout(source_row)
        self.form_layout.addWidget(QLabel('Sequence:'))
        self.form_layout.addWidget(self.sequence_txt)

        self.name_field.setFocus()

    @property
    def name(self) -> str:
        return self.name_field.text().strip()

    @property
    def segment(self) -> str:
        return self.segment_combo.currentData()

    @property
    def source(self) -> str:
        return self.source_field.text().strip()

    @property
    def sequence(self) -> str:
        return sanitize_sequence(self.sequence_txt.toPlainText())

    def field_values(self) -> dict:
        return {
            'name': self.name,
            'segment': self.segment,
            'source': self.source,
            'sequence': self.sequence,
        }
