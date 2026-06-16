from PySide6.QtCore import Qt
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

from flumut_db_editor.gui.forms.protein_form import ProteinForm
from flumut_db_editor.gui.forms.segment_form import SegmentForm
from flumut_db_editor.gui.tabs.base import BaseTab
from flumut_db_editor.models import Protein, Segment


class ProteinsTab(BaseTab):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.tab_layout.addWidget(self.tree)
        self.refresh()

    def refresh(self):
        self.tree.clear()
        segments: list[Segment] = Segment.select()
        for segment in segments:
            segment_item = QTreeWidgetItem(self.tree)
            segment_item.setText(0, segment.name)
            segment_item.setText(1, f'{len(segment.proteins)} proteins, {len(segment.references)} references')
            segment_item.setData(0, Qt.ItemDataRole.UserRole, segment)
            proteins: list[Protein] = Protein.select().where(Protein.segment == segment)
            for protein in proteins:
                protein_item = QTreeWidgetItem(segment_item)
                protein_item.setText(0, protein.name)
                protein_item.setText(1, f'{len(protein.mutations)} mutations')
                protein_item.setData(0, Qt.ItemDataRole.UserRole, protein)

    def on_new_requested(self):
        form = SegmentForm(self, None)
        if form.exec():
            self.refresh()

    def on_edit_requested(self):
        instance = self.get_selected_instance()
        if instance is None:
            return

        if isinstance(instance, Segment):
            form = SegmentForm(self, instance)
        elif isinstance(instance, Protein):
            form = ProteinForm(self, instance)
        else:
            return

        if form.exec():
            self.refresh()

    def on_delete_requested(self):
        instance = self.get_selected_instance()
        if instance is None:
            return

        if isinstance(instance, Segment):
            form = SegmentForm(self, instance)
        elif isinstance(instance, Protein):
            form = ProteinForm(self, instance)
        else:
            return

        if form.exec():
            self.refresh()

    def get_selected_instance(self) -> Segment | Protein | None:
        if item := self.tree.currentItem():
            return item.data(0, Qt.ItemDataRole.UserRole)
        return None
