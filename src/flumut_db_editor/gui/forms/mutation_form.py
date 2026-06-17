from typing import TYPE_CHECKING

from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QLineEdit, QSpinBox, QWidget

from flumut.flumutdb.models import Mutation, MutationType, Protein
from flumut_db_editor.gui.forms.base import TransactionalForm


class MutationForm(TransactionalForm):
    model = Mutation

    def __init__(self, parent: QWidget | None = None, instance: Mutation | None = None, protein: Protein | None = None) -> None:
        self.force_protein = protein
        if instance:
            self.force_protein = instance.protein

        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: Mutation | None

    def init_ui(self):
        super().init_ui()

        self.name_field = QLineEdit()
        self.type_combo = QComboBox()
        self.protein_combo = QComboBox()
        self.default_pos_field = QSpinBox()

        for t in MutationType:
            self.type_combo.addItem(t.name, t)

        proteins: list[Protein] = Protein.select()
        for protein in proteins:
            self.protein_combo.addItem(protein.name, protein)

        self.default_pos_field.setMinimum(0)

        if self.instance:
            self.name_field.setText(self.instance.name)
            self.type_combo.setCurrentText(self.instance.type)
            if self.instance.default_position:
                self.default_pos_field.setValue(self.instance.default_position)

        if self.force_protein:
            self.protein_combo.setCurrentText(self.force_protein.name)
            self.protein_combo.setEnabled(False)

        name_row = QHBoxLayout()
        name_row.addWidget(QLabel('Name:'))
        name_row.addWidget(self.name_field)

        type_row = QHBoxLayout()
        type_row.addWidget(QLabel('Type:'))
        type_row.addWidget(self.type_combo, 1)

        protein_row = QHBoxLayout()
        protein_row.addWidget(QLabel('Protein:'))
        protein_row.addWidget(self.protein_combo, 1)

        default_pos_row = QHBoxLayout()
        default_pos_row.addWidget(QLabel('Default position:'))
        default_pos_row.addWidget(self.default_pos_field, 1)

        self.form_layout.addLayout(name_row)
        self.form_layout.addLayout(type_row)
        self.form_layout.addLayout(protein_row)
        self.form_layout.addLayout(default_pos_row)
        self.form_layout.addStretch()

    @property
    def name(self) -> str:
        return self.name_field.text().strip()

    @property
    def type(self) -> str:
        return self.type_combo.currentText().strip()

    @property
    def protein(self) -> Protein:
        return self.protein_combo.currentData()

    @property
    def default_position(self) -> int | None:
        if value := self.default_pos_field.value() > 0:
            return value
        return None

    def validate(self) -> bool:
        if not self.check_unique_required('name', self.name, 'Name', self.name_field):
            return False
        return True

    def field_values(self) -> dict:
        return {
            'name': self.name,
            'type': self.type,
            'protein': self.protein,
            'default_position': self.default_position,
        }
