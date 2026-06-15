from PySide6.QtWidgets import QComboBox, QLabel, QSpinBox

from flumut.flumutdb.models import Annotation, Protein
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class AnnotationForm(TransactionalForm):
    model = Annotation

    def __init__(self, parent=None, instance=None, reference_id=None):
        self.reference_id = reference_id
        super().__init__(parent, 'Annotation', instance)

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

        if self.instance:
            self.protein_combo.setCurrentText(self.instance.protein.name)
            self.start_field.setValue(self.instance.start)
            self.end_field.setValue(self.instance.end)

    def validate(self) -> bool:
        if self.protein_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Protein', 'Please select a protein.')
            return False
        if self.start_field.value() > self.end_field.value():
            ValidationErrorDialog.show_validation_error(self, 'Range', 'Start position must be <= End position.')
            return False
        return True

    def field_values(self) -> dict:
        return {
            'protein_id': self.protein_combo.currentData(),
            'start': self.start_field.value(),
            'end': self.end_field.value(),
        }

    def create_values(self) -> dict:
        return {**self.field_values(), 'reference_id': self.reference_id}
