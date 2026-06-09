from PySide6.QtWidgets import QComboBox, QLabel, QLineEdit, QSpinBox

from flumut.flumutdb.models import Mapping, Reference
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import BaseForm


class MappingForm(BaseForm):
    def __init__(self, parent=None, mapping=None, mutation_id=None):
        self.mapping = mapping
        self.mutation_id = mutation_id
        super().__init__(parent, 'Mapping')

    def init_ui(self):
        super().init_ui()
        self.reference_combo = QComboBox()
        self.position_field = QSpinBox()
        self.alteration_field = QLineEdit()

        references = Reference.select()
        for ref in references:
            self.reference_combo.addItem(ref.name, ref.id)

        self.position_field.setMinimum(0)

        self.form_layout.insertWidget(0, QLabel('Reference:'))
        self.form_layout.insertWidget(1, self.reference_combo)
        self.form_layout.insertWidget(2, QLabel('Position:'))
        self.form_layout.insertWidget(3, self.position_field)
        self.form_layout.insertWidget(4, QLabel('Alteration:'))
        self.form_layout.insertWidget(5, self.alteration_field)

        if self.mapping:
            self.reference_combo.setCurrentText(self.mapping.reference.name)
            self.position_field.setValue(self.mapping.position)
            self.alteration_field.setText(self.mapping.alteration or '')

    def validate(self) -> bool:
        if self.reference_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Reference', 'Please select a reference.')
            return False
        return True

    def save_to_db(self) -> None:
        reference_id = self.reference_combo.currentData()
        position = self.position_field.value()
        alteration = self.alteration_field.text().strip()

        if self.mapping:
            self.mapping.reference_id = reference_id
            self.mapping.position = position
            self.mapping.alteration = alteration or None
            self.mapping.save()
        else:
            Mapping.create(
                mutation_id=self.mutation_id,
                reference_id=reference_id,
                position=position,
                alteration=alteration or None,
            )
