from PySide6.QtWidgets import QComboBox, QLabel, QLineEdit

from flumut_db_editor.gui.forms.base import BaseForm
from flumut_db_editor.models import MutationType, Protein


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

    def get_data(self):
        return {
            'name': self.name_field.text(),
            'type': self.type_combo.currentData(),
            'protein_id': self.protein_combo.currentData(),
        }
