from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.gui.base import HierarchicalTab
from flumut_db_editor.models import Marker, Mutation


class MarkersTab(HierarchicalTab):
    def __init__(self):
        super().__init__()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.load_data()

    def load_data(self):
        self.tree.clear()
        markers: list[Marker] = Marker.select()
        for marker in markers:
            marker_item = QTreeWidgetItem(self.tree)
            marker_item.setText(0, marker.name)
            marker_item.setText(1, '')
            marker_item.setData(0, 0x0100, ('marker', marker.id))

            mutations = marker.mutations
            for mutation in mutations:
                mutation_item = QTreeWidgetItem(marker_item)
                mutation_item.setText(
                    0, f'{mutation.name} ({mutation.type})'
                )
                mutation_item.setText(1, f'Protein: {mutation.protein.name}')
                mutation_item.setData(0, 0x0100, ('mutation', mutation.id))
