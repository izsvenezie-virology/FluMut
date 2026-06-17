from PySide6.QtCore import Qt
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.mapping_form import MappingForm
from flumut_db_editor.gui.forms.mutation_form import MutationForm
from flumut_db_editor.gui.tabs.base import BaseTab
from flumut_db_editor.models import Mapping, Mutation


class MutationsTab(BaseTab):
    def __init__(self) -> None:
        super().__init__()
        self._init_ui()

    def _init_ui(self) -> None:
        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.tab_layout.addWidget(self.tree)
        self.refresh()

    def refresh(self):
        self.tree.clear()
        mutations: list[Mutation] = Mutation.select()
        for mutation in mutations:
            mut_item = QTreeWidgetItem(self.tree)
            mut_item.setText(0, mutation.name)
            mut_item.setText(1, f'Type: {mutation.type} | Protein: {mutation.protein.name} | Markers: {len(mutation.markers)}')
            mut_item.setData(0, Qt.ItemDataRole.UserRole, mutation)

            for mapping in mutation.mappings:
                mapping_item = QTreeWidgetItem(mut_item)
                mapping_item.setText(0, f'{mapping.reference.name}:{mapping.position}')
                mapping_item.setText(1, mapping.alteration or '')
                mapping_item.setData(0, Qt.ItemDataRole.UserRole, mapping)

    def on_new_requested(self):
        form = MutationForm(self, None)
        if form.exec():
            self.refresh()

    def on_edit_requested(self):
        instance = self.get_selected_instance()
        if isinstance(instance, Mutation):
            form = MutationForm(self, instance)
        elif isinstance(instance, Mapping):
            form = MappingForm(self, instance)
        else:
            return

        if form.exec():
            self.refresh()

    def on_delete_requested(self):
        instance = self.get_selected_instance()
        if instance:
            if DeleteForm.confirm_and_delete(instance, self):
                self.refresh()

    def get_selected_instance(self) -> Mutation | Mapping | None:
        if item := self.tree.currentItem():
            return item.data(0, Qt.ItemDataRole.UserRole)
        return None
