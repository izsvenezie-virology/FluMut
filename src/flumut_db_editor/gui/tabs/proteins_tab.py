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

    def refresh(self):
        self.tree.clear()
        segments: list[Segment] = Segment.select()
        for segment in sorted(segments, key=lambda x: x.order):
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
            segment_item.setExpanded(self._expansion_status.get(segment, False))

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
                self.refresh()

    def get_selected_instance(self) -> Segment | Protein | None:
        if item := self.tree.currentItem():
            return item.data(0, Qt.ItemDataRole.UserRole)
        return None

    def select_by_instance(self, instance: Segment | Protein) -> None:
        """Select the tree item backed by ``instance`` (a segment or one of its
        proteins), expanding and scrolling to it. No-op if it is not found."""
        for i in range(self.tree.topLevelItemCount()):
            segment_item = self.tree.topLevelItem(i)
            if segment_item is None:
                continue
            if segment_item.data(0, Qt.ItemDataRole.UserRole) == instance:
                self.tree.setCurrentItem(segment_item)
                self.tree.scrollToItem(segment_item)
                return
            for j in range(segment_item.childCount()):
                protein_item = segment_item.child(j)
                if protein_item is None:
                    continue
                if protein_item.data(0, Qt.ItemDataRole.UserRole) == instance:
                    segment_item.setExpanded(True)
                    self.tree.setCurrentItem(protein_item)
                    self.tree.scrollToItem(protein_item)
                    return

    def move_segment(self, segment: Segment, up: bool) -> None:
        old_pos = segment.order
        new_pos = old_pos + (-1 if up else 1)
        switch = Segment.get_or_none(Segment.order == new_pos)
        if switch is None:
            return
        segment.order = new_pos
        switch.order = old_pos
        segment.save()
        switch.save()

        self.refresh()
        self.select_by_instance(segment)

    def move_up(self) -> None:
        instance = self.get_selected_instance()
        if isinstance(instance, Segment):
            self.move_segment(instance, True)

    def move_down(self) -> None:
        instance = self.get_selected_instance()
        if isinstance(instance, Segment):
            self.move_segment(instance, False)
