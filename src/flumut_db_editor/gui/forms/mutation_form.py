from PySide6.QtWidgets import QComboBox, QLabel, QLineEdit

from flumut.flumutdb.models import Mutation, MutationType, Protein
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class MutationForm(TransactionalForm):
    model = Mutation

    def __init__(self, parent=None, instance=None):
        super().__init__(parent, 'Mutation', instance)

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

        if self.instance:
            self.name_field.setText(self.instance.name)
            self.protein_combo.setCurrentText(self.instance.protein.name)

    def validate(self) -> bool:
        name = self.name_field.text().strip()
        if not name:
            ValidationErrorDialog.show_validation_error(self, 'Name', 'Name cannot be empty.')
            return False
        if self.protein_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Protein', 'Please select a protein.')
            return False
        if not self.instance:
            if Mutation.select().where(Mutation.name == name).exists():
                ValidationErrorDialog.show_validation_error(self, 'Name', 'A mutation with this name already exists.')
                return False
        return True

    def field_values(self) -> dict:
        return {
            'name': self.name_field.text().strip(),
            'type': self.type_combo.currentData(),
            'protein_id': self.protein_combo.currentData(),
        }
