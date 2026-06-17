from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb.models import Annotation, Reference, Segment
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.annotation_form import AnnotationForm
from flumut_db_editor.gui.forms.base import TransactionalForm
from flumut_db_editor.gui.forms.delete_form import DeleteForm


class ReferenceForm(TransactionalForm):
    model = Reference

    def __init__(self, parent: QWidget | None = None, instance: Reference | None = None, segment: Segment | None = None) -> None:
        self.force_segment = segment
        if instance:
            self.force_segment = instance.segment

        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: Reference | None

    def init_ui(self):
        super().init_ui()

        self.name_field = QLineEdit()
        self.segment_combo = QComboBox()
        self.source_field = QLineEdit()
        self.sequence_txt = QTextEdit()

        segments: list[Segment] = Segment.select()
        for segment in segments:
            self.segment_combo.addItem(segment.name, segment)

        if self.instance:
            self.name_field.setText(self.instance.name)
            self.source_field.setText(self.instance.source)
            self.sequence_txt.setPlainText(self.instance.sequence)

        if self.force_segment:
            self.segment_combo.setCurrentText(self.force_segment.name)
            self.segment_combo.setEnabled(False)

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
        self.form_layout.addStretch()

        annotations_row = QHBoxLayout()
        self.annotations_list = QListWidget()
        annotations_row.addWidget(self.annotations_list, stretch=1)

        btn_col = QVBoxLayout()
        self.add_annotation_btn = QPushButton('Add')
        self.edit_annotation_btn = QPushButton('Edit')
        self.remove_annotation_btn = QPushButton('Remove')
        btn_col.addWidget(self.add_annotation_btn)
        btn_col.addWidget(self.edit_annotation_btn)
        btn_col.addWidget(self.remove_annotation_btn)
        btn_col.addStretch()
        annotations_row.addLayout(btn_col)

        self.form_layout.addLayout(annotations_row)

        self.add_annotation_btn.clicked.connect(self._add_annotation)
        self.edit_annotation_btn.clicked.connect(self._edit_annotation)
        self.remove_annotation_btn.clicked.connect(self._remove_annotation)

        self._set_annotations_enabled(self.instance is not None)
        self._refresh_list()
        self.name_field.setFocus()

    def _set_annotations_enabled(self, enable: bool) -> None:
        self.annotations_list.setEnabled(enable)
        self.edit_annotation_btn.setEnabled(enable)
        self.remove_annotation_btn.setEnabled(enable)

    def _add_annotation(self) -> None:
        if self.instance is None:
            self.save_to_db()
            if self.instance is None:
                return
        self._set_annotations_enabled(True)

        form = AnnotationForm(self.instance, self, None)
        if form.exec():
            self._refresh_list()

    def _edit_annotation(self) -> None:
        if annotation := self._get_selected_annotation():
            form = AnnotationForm(annotation.reference, self, annotation)
            if form.exec():
                self._refresh_list()

    def _remove_annotation(self) -> None:
        if annotation := self._get_selected_annotation():
            if DeleteForm.confirm_and_delete(annotation, self):
                self._refresh_list()

    def _refresh_list(self) -> None:
        self.annotations_list.clear()
        if self.instance:
            for a in self.instance.annotations:
                item = QListWidgetItem(f'{a.protein}: {a.start}-{a.end}')
                item.setData(Qt.ItemDataRole.UserRole, a)
                self.annotations_list.addItem(item)

    def _get_selected_annotation(self) -> Annotation | None:
        if item := self.annotations_list.currentItem():
            return item.data(Qt.ItemDataRole.UserRole)

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
        return ''.join(self.sequence_txt.toPlainText().split()).upper()

    def validate(self) -> bool:
        if not self.check_unique_required('name', self.name, 'Name', self.name_field):
            return False
        if not self.check_unique_required('source', self.source, 'Source', self.source_field):
            return False
        if not self.check_unique_required('sequence', self.sequence, 'Sequence', self.sequence_txt):
            return False
        if unknown_nucleotides := set(self.sequence).difference(set(['A', 'T', 'C', 'G'])):
            ValidationErrorDialog.show_validation_error(
                self, 'Sequence', f'Reference sequences must contain only ATCG. Found unknown nucleotides: {", ".join(unknown_nucleotides)}'
            )
            return False
        return True

    def field_values(self) -> dict:
        return {
            'name': self.name,
            'segment': self.segment,
            'source': self.source,
            'sequence': self.sequence,
        }
