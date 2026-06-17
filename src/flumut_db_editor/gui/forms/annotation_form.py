from typing import TYPE_CHECKING

from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QSpinBox, QWidget

from flumut.flumutdb.models import Annotation, Protein, Reference
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class AnnotationForm(TransactionalForm):
    model = Annotation

    def __init__(
        self,
        reference: Reference,
        parent: QWidget | None = None,
        instance: Annotation | None = None,
        protein: Protein | None = None,
    ) -> None:
        self.force_reference = reference
        self.force_protein = protein
        if instance:
            self.force_reference = instance.reference
            self.force_protein = instance.protein

        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: Annotation | None

    def init_ui(self):
        super().init_ui()

        self.reference_combo = QComboBox()
        self.protein_combo = QComboBox()
        self.start_field = QSpinBox()
        self.end_field = QSpinBox()

        references: list[Reference] = Reference.select()
        for reference in references:
            self.reference_combo.addItem(reference.name, reference)

        proteins = self.force_reference.segment.proteins
        for protein in proteins:
            self.protein_combo.addItem(protein.name, protein)

        self.start_field.setMinimum(1)
        self.start_field.setMaximum(1_000_000)
        self.end_field.setMinimum(1)
        self.end_field.setMaximum(1_000_000)

        if self.instance:
            self.start_field.setValue(self.instance.start)
            self.end_field.setValue(self.instance.end)

        if self.force_reference:
            self.reference_combo.setCurrentText(self.force_reference.name)
            self.reference_combo.setEnabled(False)
            self.start_field.setMaximum(len(self.force_reference.sequence))
            self.end_field.setMaximum(len(self.force_reference.sequence))

        if self.force_protein:
            self.protein_combo.setCurrentText(self.force_protein.name)
            self.protein_combo.setEnabled(False)

        reference_row = QHBoxLayout()
        reference_row.addWidget(QLabel('Reference:'))
        reference_row.addWidget(self.reference_combo)

        protein_row = QHBoxLayout()
        protein_row.addWidget(QLabel('Protein:'))
        protein_row.addWidget(self.protein_combo)

        start_row = QHBoxLayout()
        start_row.addWidget(QLabel('Start:'))
        start_row.addWidget(self.start_field)

        end_row = QHBoxLayout()
        end_row.addWidget(QLabel('End:'))
        end_row.addWidget(self.end_field)

        self.form_layout.addLayout(reference_row)
        self.form_layout.addLayout(protein_row)
        self.form_layout.addLayout(start_row)
        self.form_layout.addLayout(end_row)

    @property
    def reference(self) -> Reference:
        return self.reference_combo.currentData()

    @property
    def protein(self) -> Protein:
        return self.protein_combo.currentData()

    @property
    def start(self) -> int:
        return self.start_field.value()

    @property
    def end(self) -> int:
        return self.end_field.value()

    def validate(self) -> bool:
        if self.start > self.end:
            ValidationErrorDialog.show_validation_error(self, 'Range', 'Start position cannot be higher than end position.')
            return False
        if self.protein_combo.currentData() is None:
            ValidationErrorDialog.show_validation_error(self, 'Protein', 'The protein cannot be empty.')
            return False
        if self.reference.segment != self.protein.segment:
            ValidationErrorDialog.show_validation_error(self, 'Segment', 'Reference and protein belongs to different segments.')
            return False
        return True

    def field_values(self) -> dict:
        return {
            'reference': self.reference,
            'protein': self.protein,
            'start': self.start,
            'end': self.end,
        }
