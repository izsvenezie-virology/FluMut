from collections.abc import Sequence

from PySide6.QtWidgets import QPushButton, QTreeWidgetItem

from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.protein_form import ProteinForm
from flumut_db_editor.gui.forms.segment_form import SegmentForm
from flumut_db_editor.gui.tabs.base import BaseSortableTreeTab
from flumut_db_editor.models import Protein, Segment


class ProteinsTab(BaseSortableTreeTab[Segment | Protein]):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        self.refresh()

        self.new_protein_btn = QPushButton('New protein')

        self.header.itemAt(0).widget().setText('New segment')  # type: ignore
        self.header.insertWidget(1, self.new_protein_btn)

        self.new_protein_btn.clicked.connect(self.on_new_protein_requested)

    def refresh(self, selected: Segment | Protein | None = None):
        self.tree.clear()
        segments: Sequence[Segment] = sorted(Segment.select())
        for segment in segments:
            segment_item = QTreeWidgetItem(self.tree)
            segment_item.setText(0, segment.name)
            segment_item.setText(1, f'{len(segment.proteins)} proteins, {len(segment.references)} references')
            self.set_data(segment_item, segment)
            segment_item.setExpanded(self.is_expanded(segment))
            if segment == selected:
                self.tree.setCurrentItem(segment_item)
                self.tree.scrollToItem(segment_item)

            proteins: Sequence[Protein] = sorted(segment.proteins)
            for protein in proteins:
                protein_item = QTreeWidgetItem(segment_item)
                protein_item.setText(0, protein.name)
                protein_item.setText(1, f'{len(protein.mutations)} mutations')
                self.set_data(protein_item, protein)
                if protein == selected:
                    segment_item.setExpanded(True)
                    self.tree.setCurrentItem(protein_item)
                    self.tree.scrollToItem(protein_item)

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
        if instance:
            if DeleteForm.confirm_and_delete(instance, self):
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
