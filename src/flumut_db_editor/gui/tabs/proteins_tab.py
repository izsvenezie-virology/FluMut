from collections.abc import Sequence

from PySide6.QtWidgets import QPushButton

from flumut.flumutdb import loader
from flumut.flumutdb.models import Protein, Segment
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.protein_form import ProteinForm
from flumut_db_editor.gui.forms.segment_form import SegmentForm
from flumut_db_editor.gui.tabs.base import BaseSortableTreeTab


class ProteinsTab(BaseSortableTreeTab[Segment | Protein]):
    def __init__(self):
        super().__init__()
        self.__init_ui()

    def __init_ui(self):
        self.new_protein_btn = QPushButton('New protein')
        self.new_protein_btn.clicked.connect(self.on_new_protein_requested)

        self.new_btn.setText('New segment')  # type: ignore
        self.header.insertWidget(1, self.new_protein_btn)

        self.refresh()

    def populate(self, selected=None):
        self.tree.clear()
        segments = sorted(loader.get(Segment))
        for segment in segments:
            segment_item = self.create_item(
                segment, [f'Segment: {segment.name}', f'{len(segment.proteins)} proteins, {len(segment.references)} references']
            )
            if segment == selected:
                self.set_selected_item(segment_item)

            proteins: Sequence[Protein] = sorted(segment.proteins)
            for protein in proteins:
                protein_item = self.create_item(protein, [f'Protein: {protein.name}', f'{len(protein.mutations)} mutations'], segment_item)
                if protein == selected:
                    self.set_selected_item(protein_item)

    def on_new_requested(self):
        form = SegmentForm(self, None)
        if form.exec() and form.instance:
            lst = self.get_sorted_list(form.instance)

            if selected := self.get_selected_segment():
                lst.remove(form.instance)
                lst.insert(lst.index(selected) + 1, form.instance)

            self.update_order(lst)
            self.refresh(form.instance)

    def on_new_protein_requested(self):
        form = ProteinForm(self, None, self.get_selected_segment())
        if form.exec() and form.instance:
            lst = self.get_sorted_list(form.instance)

            selected_protein = self.get_selected_protein()
            if selected_protein and selected_protein.segment == form.instance.segment:
                lst.remove(form.instance)
                lst.insert(lst.index(selected_protein) + 1, form.instance)

            self.update_order(lst)
            self.refresh(form.instance)

    def on_edit_requested(self):
        if instance := self.get_selected_instance():
            form = self.get_form(instance)
            if form.exec() and form.instance:
                self.refresh(form.instance)

    def on_delete_requested(self):
        instance = self.get_selected_instance()
        if instance and DeleteForm.confirm_and_delete(instance, self):
            list_to_order = self.get_sorted_list(instance)
            self.update_order(list_to_order)
            self.refresh()

    def get_selected_segment(self) -> Segment | None:
        selected = self.get_selected_instance()
        match selected:
            case Segment():
                return selected
            case Protein():
                return selected.segment
            case _:
                return None

    def get_selected_protein(self) -> Protein | None:
        selected = self.get_selected_instance()
        match selected:
            case Protein():
                return selected
            case _:
                return None

    def get_sorted_list(self, instance: Segment | Protein) -> list[Segment | Protein]:
        match instance:
            case Segment():
                return sorted(Segment.select())
            case Protein():
                return sorted(instance.segment.proteins)
            case _:
                return []

    def get_form(self, instance: Segment | Protein) -> SegmentForm | ProteinForm:
        match instance:
            case Segment():
                return SegmentForm(self, instance)
            case Protein():
                return ProteinForm(self, instance)
