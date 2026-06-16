from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.marker_form import MarkerForm
from flumut_db_editor.gui.tabs.base import HierarchicalTab
from flumut_db_editor.models import Marker


class MarkersTab(HierarchicalTab):
    def __init__(self):
        super().__init__()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.new_button.clicked.connect(self.handle_new)
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
                mutation_item.setText(0, f'{mutation.name} ({mutation.type})')
                mutation_item.setText(1, f'Protein: {mutation.protein.name}')
                mutation_item.setData(0, 0x0100, ('marker_mutation', marker.id, mutation.id))

    def handle_new(self):
        form = MarkerForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Marker created successfully.')
            self.load_data()

    def handle_edit(self):
        item = self.get_selected_item()
        if item is None:
            return
        data = item.data(0, 0x0100)
        if not data:
            return

        if data[0] == 'marker':
            marker_id = data[1]
            marker = Marker.get_by_id(marker_id)
            form = MarkerForm(self, marker)
            if form.exec():
                SuccessNotification.show_success(self, 'Marker updated successfully.')
                self.load_data()

    def handle_delete(self):
        item = self.get_selected_item()
        if item is None:
            return
        data = item.data(0, 0x0100)
        if not data:
            return

        if data[0] == 'marker':
            marker_id = data[1]
            if delete_with_confirmation(self, 'Marker', marker_id, DatabaseOperations.delete_marker):
                self.load_data()
        elif data[0] == 'marker_mutation':
            marker_id = data[1]
            mutation_id = data[2]
            if delete_with_confirmation(
                self,
                'Mutation Association',
                mutation_id,
                lambda mid=mutation_id: DatabaseOperations.remove_marker_mutation_association(marker_id, mid),
            ):
                self.load_data()
