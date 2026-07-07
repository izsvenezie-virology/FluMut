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
        self.refresh()

        self.up_btn = QPushButton('Move up')
        self.down_btn = QPushButton('Move down')

        self.header.addSpacerItem(QSpacerItem(10, 10))
        self.header.addWidget(self.up_btn)
        self.header.addWidget(self.down_btn)

        self.up_btn.clicked.connect(self.move_up)
        self.down_btn.clicked.connect(self.move_down)
        self.tree.itemExpanded.connect(self.on_item_expanded)
        self.tree.itemCollapsed.connect(self.on_item_collapsed)

    def on_item_expanded(self, item: QTreeWidgetItem):
        if segment := item.data(0, Qt.ItemDataRole.UserRole):
            self._expansion_status[segment] = True

    def on_item_collapsed(self, item: QTreeWidgetItem):
        if segment := item.data(0, Qt.ItemDataRole.UserRole):
            self._expansion_status[segment] = False

    def refresh(self, selected: Segment | Protein | None = None):
        self.tree.clear()
        segments: list[Segment] = sorted(Segment.select())
        for segment in segments:
            segment_item = QTreeWidgetItem(self.tree)
            segment_item.setText(0, segment.name)
            segment_item.setText(1, f'{len(segment.proteins)} proteins, {len(segment.references)} references')
            segment_item.setData(0, Qt.ItemDataRole.UserRole, segment)
            segment_item.setExpanded(self._expansion_status.get(segment, False))
            if segment == selected:
                self.tree.setCurrentItem(segment_item)
                self.tree.scrollToItem(segment_item)

            proteins: list[Protein] = sorted(segment.proteins)
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
        if form.exec():
            self.refresh()

    def on_edit_requested(self):
        instance = self.get_selected_instance()
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
        if instance:
            if DeleteForm.confirm_and_delete(instance, self):
                if isinstance(instance, Segment):
                    list_to_order: list[Segment | Protein] = sorted(Segment.select())
                elif isinstance(instance, Protein):
                    list_to_order = sorted(instance.segment.proteins)
                else:
                    return
                self.update_order(list_to_order)
                self.refresh()

    def get_selected_instance(self) -> Segment | Protein | None:
        if item := self.tree.currentItem():
            return item.data(0, Qt.ItemDataRole.UserRole)
        return None

    def move_instance(self, instance: Segment | Protein, up: bool) -> None:
        if isinstance(instance, Segment):
            list_to_move: list[Segment | Protein] = sorted(Segment.select())
        elif isinstance(instance, Protein):
            list_to_move = sorted(instance.segment.proteins)
        else:
            return

        idx = list_to_move.index(instance)
        new_idx = idx + (-1 if up else 1)
        if not 0 <= new_idx < len(list_to_move):
            return

        list_to_move.insert(new_idx, list_to_move.pop(idx))
        self.update_order(list_to_move)
        self.refresh(instance)

    def update_order(self, values: list[Segment | Protein]) -> None:
        for order, instance in enumerate(values):
            instance.order = order
            instance.save()

    def move_up(self) -> None:
        instance = self.get_selected_instance()
        if isinstance(instance, Segment | Protein):
            self.move_instance(instance, True)

    def move_down(self) -> None:
        instance = self.get_selected_instance()
        if isinstance(instance, Segment | Protein):
            self.move_instance(instance, False)
