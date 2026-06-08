from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.gui.base import HierarchicalTab
from flumut_db_editor.models import Protein, Segment


class ProteinsTab(HierarchicalTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        self.tree.clear()
        segments = Segment.select()
        for segment in segments:
            segment_item = QTreeWidgetItem(self.tree)
            segment_item.setText(0, segment.name)
            proteins = Protein.select().where(Protein.segment == segment)
            for protein in proteins:
                protein_item = QTreeWidgetItem(segment_item)
                protein_item.setText(0, protein.name)
