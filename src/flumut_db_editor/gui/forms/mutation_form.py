from PySide6.QtWidgets import QComboBox, QLabel, QLineEdit

from flumut.flumutdb.models import Mutation, MutationType, Protein
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import BaseForm


class MutationForm(BaseForm):
    def __init__(self, parent=None, mutation=None):
        self.mutation = mutation
        super().__init__(parent, 'Mutation')

    def init_ui(self):
        super().init_ui()
        self.name_field = QLineEdit()
        self.type_combo = QComboBox()
        self.protein_combo = QComboBox()

        types = [MutationType.SNP]
        for t in types:
            self.type_combo.addItem(t.value if hasattr(t, 'value') else str(t), t)

        proteins = Protein.select()
        for protein in proteins:
            self.protein_combo.addItem(protein.name, protein.id)

        self.form_layout.insertWidget(0, QLabel('Name:'))
        self.form_layout.insertWidget(1, self.name_field)
        self.form_layout.insertWidget(2, QLabel('Type:'))
        self.form_layout.insertWidget(3, self.type_combo)
        self.form_layout.insertWidget(4, QLabel('Protein:'))
        self.form_layout.insertWidget(5, self.protein_combo)

        if self.mutation:
            self.name_field.setText(self.mutation.name)
            self.protein_combo.setCurrentText(self.mutation.protein.name)

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(
                self, 'Name', 'Name cannot be empty.'
            )
            return False
        if self.protein_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(
                self, 'Protein', 'Please select a protein.'
            )
            return False
        if not self.mutation:
            if Mutation.select().where(Mutation.name == name).exists():
                ValidationErrorDialog.show_validation_error(
                    self, 'Name', 'A mutation with this name already exists.'
                )
                return False
        return True

    def save_to_db(self) -> None:
        name = self.name_field.text().strip()
        type_val = self.type_combo.currentData()
        protein_id = self.protein_combo.currentData()
        if self.mutation:
            self.mutation.name = name
            self.mutation.type = type_val
            self.mutation.protein_id = protein_id
            self.mutation.save()
        else:
            Mutation.create(name=name, type=type_val, protein_id=protein_id)

