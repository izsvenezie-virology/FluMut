from collections.abc import Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QPushButton, QSpacerItem, QTreeWidget, QTreeWidgetItem

from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.protein_form import ProteinForm
from flumut_db_editor.gui.forms.segment_form import SegmentForm
from flumut_db_editor.gui.tabs.base import BaseTab
from flumut_db_editor.models import Protein, Segment


class ProteinsTab(BaseTab):
    def __init__(self):
        super().__init__()

        self._expansion_status: dict[Segment, bool] = {}

        self._init_ui()

    def _init_ui(self):
        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.tab_layout.addWidget(self.tree)
        self.refresh_tree()

        self.new_protein_btn = QPushButton('New protein')
        self.up_btn = QPushButton('Move up')
        self.down_btn = QPushButton('Move down')

        self.header.itemAt(0).widget().setText('New segment')  # type: ignore
        self.header.insertWidget(1, self.new_protein_btn)
        self.header.addSpacerItem(QSpacerItem(10, 10))
        self.header.addWidget(self.up_btn)
        self.header.addWidget(self.down_btn)

        self.new_protein_btn.clicked.connect(self.on_new_protein_requested)
        self.up_btn.clicked.connect(self.on_move_up_requested)
        self.down_btn.clicked.connect(self.on_move_down_requested)
        self.tree.itemExpanded.connect(self.on_item_expanded)
        self.tree.itemCollapsed.connect(self.on_item_collapsed)

    def refresh_tree(self, selected: Segment | Protein | None = None):
        self.tree.clear()
        segments: Sequence[Segment] = sorted(Segment.select())
        for segment in segments:
            segment_item = QTreeWidgetItem(self.tree)
            segment_item.setText(0, segment.name)
            segment_item.setText(1, f'{len(segment.proteins)} proteins, {len(segment.references)} references')
            segment_item.setData(0, Qt.ItemDataRole.UserRole, segment)
            segment_item.setExpanded(self._expansion_status.get(segment, False))
            if segment == selected:
                self.tree.setCurrentItem(segment_item)
                self.tree.scrollToItem(segment_item)

            proteins: Sequence[Protein] = sorted(segment.proteins)
            for protein in proteins:
                protein_item = QTreeWidgetItem(segment_item)
                protein_item.setText(0, protein.name)
                protein_item.setText(1, f'{len(protein.mutations)} mutations')
                protein_item.setData(0, Qt.ItemDataRole.UserRole, protein)
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
            self.refresh_tree(form.instance)

    def on_new_protein_requested(self):
        form = ProteinForm(self, None, self.get_selected_segment())
        if form.exec() and form.instance:
            lst = self.get_sorted_list(form.instance)

            selected_protein = self.get_selected_protein()
            if selected_protein and selected_protein.segment == form.instance.segment:
                lst.remove(form.instance)
                lst.insert(lst.index(selected_protein) + 1, form.instance)

            self.update_order(lst)
            self.refresh_tree(form.instance)

    def on_edit_requested(self):
        match instance := self.get_selected_instance():
            case Segment():
                form = SegmentForm(self, instance)
            case Protein():
                form = ProteinForm(self, instance)
            case _:
                return

        if form.exec():
            self.refresh_tree()

    def on_delete_requested(self):
        instance = self.get_selected_instance()
        if instance:
            if DeleteForm.confirm_and_delete(instance, self):
                list_to_order = self.get_sorted_list(instance)
                self.update_order(list_to_order)
                self.refresh_tree()

    def on_move_up_requested(self) -> None:
        self.move_selected_instance(True)

    def on_move_down_requested(self) -> None:
        self.move_selected_instance(False)

    def on_item_expanded(self, item: QTreeWidgetItem):
        if segment := item.data(0, Qt.ItemDataRole.UserRole):
            self._expansion_status[segment] = True

    def on_item_collapsed(self, item: QTreeWidgetItem):
        if segment := item.data(0, Qt.ItemDataRole.UserRole):
            self._expansion_status[segment] = False

    def get_selected_instance(self) -> Segment | Protein | None:
        if item := self.tree.currentItem():
            return item.data(0, Qt.ItemDataRole.UserRole)
        return None

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

    def move_selected_instance(self, up: bool) -> None:
        instance = self.get_selected_instance()
        if not instance:
            return
        list_to_move = self.get_sorted_list(instance)
        if not list_to_move:
            return

        idx = list_to_move.index(instance)
        new_idx = idx + (-1 if up else 1)
        if not 0 <= new_idx < len(list_to_move):
            return

        list_to_move.insert(new_idx, list_to_move.pop(idx))
        self.update_order(list_to_move)
        self.refresh_tree(instance)

    def update_order(self, values: Sequence[Segment | Protein]) -> None:
        for order, instance in enumerate(values):
            instance.order = order
            instance.save()
