from PySide6.QtWidgets import QPushButton

from flumut.flumutdb.models import Protein
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.mapping_form import MappingForm
from flumut_db_editor.gui.forms.mutation_form import MutationForm
from flumut_db_editor.gui.tabs.base import BaseSortableTreeTab
from flumut_db_editor.models import Mapping, Mutation


class MutationsTab(BaseSortableTreeTab[Protein | Mutation | Mapping]):
    def __init__(self) -> None:
        super().__init__()
        self._init_ui()

    def _init_ui(self) -> None:
        self.add_mapping_btn = QPushButton('Add mapping')
        self.add_mapping_btn.clicked.connect(self.on_add_mapping_requested)

        self.new_btn.setText('New mutation')
        self.header.insertWidget(1, self.add_mapping_btn)

        self.refresh()

    def refresh(self, selected=None):
        self.tree.clear()
        proteins: list[Protein] = sorted(Protein.select())
        for protein in proteins:
            protein_item = self.create_item(protein, [f'Protein {protein}', f'{len(protein.mutations)} mutations'])
            if protein == selected:
                self.set_selected_item(protein_item)

            for mutation in sorted(protein.mutations):
                mutation_item = self.create_item(mutation, [mutation.name], protein_item)
                if mutation == selected:
                    self.set_selected_item(mutation_item)
                for mapping in mutation.mappings:
                    mapping_item = self.create_item(mapping, [str(mapping)], mutation_item)
                    if mapping == selected:
                        self.set_selected_item(mapping_item)

    def on_new_requested(self):
        form = MutationForm(self, None)
        if form.exec() and form.instance:
            lst = self.get_sorted_list(form.instance)

            if selected := self.get_selected_mutation():
                lst.remove(form.instance)
                lst.insert(lst.index(selected) + 1, form.instance)

            self.update_order(lst)
            self.refresh(form.instance)

    def on_add_mapping_requested(self):
        pass

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
        if instance and DeleteForm.confirm_and_delete(instance, self):
            self.refresh()

    def get_selected_mutation(self) -> Mutation | None:
        selected = self.get_selected_instance()
        match selected:
            case Mutation():
                return selected
            case Mapping():
                return selected.mutation
            case _:
                return None

    def get_sorted_list(self, instance: Protein | Mutation | Mapping) -> list[Mutation]:
        match instance:
            case Mutation():
                return sorted(instance.protein.mutations)
            case _:
                return []
