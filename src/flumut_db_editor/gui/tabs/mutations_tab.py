from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.gui.base import HierarchicalTab
from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.mapping_form import MappingForm
from flumut_db_editor.gui.forms.mutation_form import MutationForm
from flumut_db_editor.models import Mapping, Mutation


class MutationsTab(HierarchicalTab):
    def __init__(self):
        super().__init__()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.new_button.clicked.connect(self.handle_new)
        self.load_data()

    def load_data(self):
        self.tree.clear()
        mutations: list[Mutation] = Mutation.select()
        for mutation in mutations:
            mut_item = QTreeWidgetItem(self.tree)
            mut_item.setText(0, mutation.name)
            mut_item.setText(1, f'Type: {mutation.type} | Protein: {mutation.protein.name}')
            mut_item.setData(0, 0x0100, ('mutation', mutation.id))

            mappings: list[Mapping] = Mapping.select().where(Mapping.mutation == mutation)
            for mapping in mappings:
                mapping_item = QTreeWidgetItem(mut_item)
                mapping_item.setText(0, f'{mapping.reference.name}:{mapping.position}')
                mapping_item.setText(1, mapping.alteration or '')
                mapping_item.setData(0, 0x0100, ('mapping', mapping.id))

    def handle_new(self):
        form = MutationForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Mutation created successfully.')
            self.load_data()

    def handle_edit(self):
        item = self.get_selected_item()
        if item is None:
            return
        item_type, item_id = self.get_item_data(item)
        if item_type == 'mutation':
            mutation = Mutation.get_by_id(item_id)
            form = MutationForm(self, mutation)
            if form.exec():
                SuccessNotification.show_success(self, 'Mutation updated successfully.')
                self.load_data()
        elif item_type == 'mapping':
            mapping = Mapping.get_by_id(item_id)
            form = MappingForm(self, mapping)
            if form.exec():
                SuccessNotification.show_success(self, 'Mapping updated successfully.')
                self.load_data()

    def handle_delete(self):
        item = self.get_selected_item()
        if item is None:
            return
        item_type, item_id = self.get_item_data(item)
        if item_type == 'mutation':
            if delete_with_confirmation(self, 'Mutation', item_id, DatabaseOperations.delete_mutation):
                self.load_data()
        elif item_type == 'mapping':
            if delete_with_confirmation(self, 'Mapping', item_id, DatabaseOperations.delete_mapping):
                self.load_data()
