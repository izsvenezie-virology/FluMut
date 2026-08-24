from collections.abc import Iterable

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QWidget,
)

from flumut.flumutdb.models import BaseModel, Mapping, Mutation, MutationType, Protein
from flumut_db_editor.gui.forms.base import MasterDetailForm


class MutationForm(MasterDetailForm[Mutation, Mapping]):
    model = Mutation

    def __init__(self, parent: QWidget | None = None, instance: Mutation | None = None, protein: Protein | None = None) -> None:
        super().__init__(parent, instance)
        self.force_protein = protein
        if instance:
            self.force_protein = instance.protein

        self._mapping_editors: list[tuple[Mapping, QSpinBox, QLineEdit, QPushButton]] = []
        self.__init_ui()

    def __init_ui(self):
        self.setMinimumSize(560, 560)

        self.name_field = QLineEdit()
        self.type_combo = QComboBox()
        self.protein_combo = QComboBox()

        for t in MutationType:
            self.type_combo.addItem(t.name, t)

        proteins: list[Protein] = sorted(Protein.select())
        for protein in proteins:
            self.protein_combo.addItem(protein.name, protein)

        if self.instance:
            self.name_field.setText(self.instance.name)
            self.type_combo.setCurrentText(self.instance.type)

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

        self.mappings_table = QTableWidget(0, 4)
        self.mappings_table.setHorizontalHeaderLabels(['Reference', 'Position', 'Alteration'])
        self.mappings_table.verticalHeader().setVisible(False)
        header = self.mappings_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)

        self.form_layout.addLayout(name_row)
        self.form_layout.addLayout(type_row)
        self.form_layout.addLayout(protein_row)
        self.form_layout.addWidget(QLabel('Mappings:'))
        self.form_layout.addWidget(self.mappings_table, 1)

        self.protein_combo.currentIndexChanged.connect(self._populate_mappings_table)

        self._populate_mappings_table()

    def load_related(self) -> Iterable[Mapping]:
        related: list[Mapping] = []
        existing = {mapping.reference: mapping for mapping in self.instance.mappings}
        if self.instance.get_id() is not None:
            for reference in self.instance.protein.segment.references:
                if reference not in existing:
                    related.append(Mapping(mutation=self.instance, reference=reference, position=0, alteration=''))
                else:
                    related.append(existing[reference])
        return related

    def instances_to_save(self) -> Iterable[BaseModel]:
        return [self.instance, *[m for m in self.related if m.position > 0]]

    def instances_to_delete(self) -> Iterable[BaseModel]:
        return [m for m in self.related if m.position == 0 and m.get_id() is not None]

    def _populate_mappings_table(self) -> None:
        self.related.clear()
        existing_mappings = {mapping.reference: mapping for mapping in self.instance.mappings}
        self.related.extend(existing_mappings.values())

        protein: Protein = self.protein_combo.currentData()
        for reference in protein.segment.references:
            if reference not in existing_mappings:
                self.related.append(Mapping(mutation=self.instance, reference=reference, position=0, alteration=''))

        self._mapping_editors = []
        mappings = sorted(self.related, key=lambda m: m.reference)
        self.mappings_table.setRowCount(len(mappings))
        for row, mapping in enumerate(mappings):
            name_item = QTableWidgetItem(mapping.reference.name)
            name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.mappings_table.setItem(row, 0, name_item)

            position_spin = QSpinBox()
            position_spin.setMinimum(0)
            position_spin.setMaximum(len(mapping.reference.sequence))
            position_spin.setValue(mapping.position)
            self.mappings_table.setCellWidget(row, 1, position_spin)

            alteration_edit = QLineEdit(mapping.alteration or '')
            self.mappings_table.setCellWidget(row, 2, alteration_edit)

            name_button = QPushButton('Set name')
            name_button.clicked.connect(lambda checked=False, mapping=mapping: self.set_mutation_name(mapping))
            self.mappings_table.setCellWidget(row, 3, name_button)

            self._mapping_editors.append((mapping, position_spin, alteration_edit, name_button))

    def populate_related(self) -> None:
        for mapping, position_spin, alteration_edit, _ in self._mapping_editors:
            mapping.mutation = self.instance
            mapping.position = position_spin.value()
            mapping.alteration = alteration_edit.text().strip()

    def set_mutation_name(self, mapping: Mapping) -> None:
        self.populate_instance()
        self.populate_related()
        self.name_field.setText(self.create_name(mapping))

    def create_name(self, mapping: Mapping) -> str:
        return f'{mapping.mutation.protein.name}:{mapping.position}{mapping.alteration}'

    @property
    def name(self) -> str:
        return self.name_field.text().strip()

    @property
    def type(self) -> str:
        return self.type_combo.currentText().strip()

    @property
    def protein(self) -> Protein:
        return self.protein_combo.currentData()

    def field_values(self) -> dict:
        return {
            'name': self.name,
            'type': self.type,
            'protein': self.protein,
        }
