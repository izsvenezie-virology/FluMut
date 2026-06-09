from PySide6.QtWidgets import QComboBox, QLabel, QSpinBox

from flumut.flumutdb.models import Annotation, Protein
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import BaseForm


class AnnotationForm(BaseForm):
    def __init__(self, parent=None, annotation=None, reference_id=None):
        self.annotation = annotation
        self.reference_id = reference_id
        super().__init__(parent, 'Annotation')

    def init_ui(self):
        super().init_ui()
        self.protein_combo = QComboBox()
        self.start_field = QSpinBox()
        self.end_field = QSpinBox()

        proteins = Protein.select()
        for protein in proteins:
            self.protein_combo.addItem(protein.name, protein.id)

        self.start_field.setMinimum(0)
        self.end_field.setMinimum(0)

        self.form_layout.insertWidget(0, QLabel('Protein:'))
        self.form_layout.insertWidget(1, self.protein_combo)
        self.form_layout.insertWidget(2, QLabel('Start:'))
        self.form_layout.insertWidget(3, self.start_field)
        self.form_layout.insertWidget(4, QLabel('End:'))
        self.form_layout.insertWidget(5, self.end_field)

        if self.annotation:
            self.protein_combo.setCurrentText(self.annotation.protein.name)
            self.start_field.setValue(self.annotation.start)
            self.end_field.setValue(self.annotation.end)

    def validate(self) -> bool:
        if self.protein_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Protein', 'Please select a protein.')
            return False
        if self.start_field.value() > self.end_field.value():
            ValidationErrorDialog.show_validation_error(self, 'Range', 'Start position must be <= End position.')
            return False
        return True

    def save_to_db(self) -> None:
        protein_id = self.protein_combo.currentData()
        start = self.start_field.value()
        end = self.end_field.value()

        if self.annotation:
            self.annotation.protein_id = protein_id
            self.annotation.start = start
            self.annotation.end = end
            self.annotation.save()
        else:
            Annotation.create(
                protein_id=protein_id,
                reference_id=self.reference_id,
                start=start,
                end=end,
            )
