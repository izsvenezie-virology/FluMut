from PySide6.QtCore import Qt
from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.gui.base import HierarchicalTab
from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.segment_form import SegmentForm
from flumut_db_editor.models import Protein, Segment


class ProteinsTab(HierarchicalTab):
    def __init__(self):
        super().__init__()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.new_button.clicked.connect(self.handle_new)
        self.load_data()

    def load_data(self):
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
                protein_item.setText(1, '')
                protein_item.setData(0, Qt.ItemDataRole.UserRole, protein)

    def handle_new(self):
        form = SegmentForm(self, None)
        if form.exec():
            self.load_data()

    def handle_edit(self):
        item = self.get_selected_item()
        if item is None:
            return
        instance = item.data(0, Qt.ItemDataRole.UserRole)

        if isinstance(instance, Segment):
            segment = instance
        elif isinstance(instance, Protein):
            segment = instance.segment
        else:
            return

        form = SegmentForm(self, segment)
        if form.exec():
            SuccessNotification.show_success(self, 'Segment updated successfully.')
            self.load_data()

    def handle_delete(self):
        item = self.get_selected_item()
        if item is None:
            return
        instance = item.data(0, Qt.ItemDataRole.UserRole)

        if isinstance(instance, Segment):
            if delete_with_confirmation(self, 'Segment', instance.get_id(), lambda _: DatabaseOperations.delete_segment(instance.name)):
                self.load_data()
        elif isinstance(instance, Protein):
            if delete_with_confirmation(self, 'Protein', instance.get_id(), DatabaseOperations.delete_protein):
                self.load_data()
