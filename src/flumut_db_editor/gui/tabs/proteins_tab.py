from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.database_operations import DatabaseOperations
from flumut_db_editor.gui.base import HierarchicalTab
from flumut_db_editor.gui.crud_mixin import delete_with_confirmation
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
        segments = Segment.select()
        for segment in segments:
            segment_item = QTreeWidgetItem(self.tree)
            segment_item.setText(0, segment.name)
            segment_item.setText(1, '')
            segment_item.setData(0, 0x0100, ('segment', segment.id))
            proteins = Protein.select().where(Protein.segment == segment)
            for protein in proteins:
                protein_item = QTreeWidgetItem(segment_item)
                protein_item.setText(0, protein.name)
                protein_item.setText(1, '')
                protein_item.setData(0, 0x0100, ('protein', protein.id))

    def handle_new(self):
        form = SegmentForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Segment created successfully.')
            self.load_data()

    def handle_edit(self):
        item = self.get_selected_item()
        if item is None:
            return
        data = item.data(0, 0x0100)
        if not data:
            return

        if data[0] == 'segment':
            segment = Segment.get_by_id(data[1])
        elif data[0] == 'protein':
            protein = Protein.get_by_id(data[1])
            segment = protein.segment
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
        data = item.data(0, 0x0100)
        if not data:
            return

        if data[0] == 'segment':
            segment = Segment.get_by_id(data[1])
            if delete_with_confirmation(self, 'Segment', segment.id, lambda _: DatabaseOperations.delete_segment(segment.name)):
                self.load_data()
        elif data[0] == 'protein':
            if delete_with_confirmation(self, 'Protein', data[1], DatabaseOperations.delete_protein):
                self.load_data()
