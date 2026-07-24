from typing import Generic, Sequence, TypeVar

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QSpacerItem,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb.models import BaseModel

ModelT = TypeVar('ModelT', bound=BaseModel)


class BaseTab(QWidget):
    """Common tab chrome: a New/Refresh header above an item view.

    Subclasses implement :meth:`create_view` to supply the table or tree; the
    CRUD behaviour (context menu, refresh, delete) comes from the mixin.
    """

    def __init__(self):
        super().__init__()
        self._init_header()

    def _init_header(self) -> None:
        self.tab_layout = QVBoxLayout(self)

        self.header = QHBoxLayout()
        self.new_button = QPushButton('New')
        self.edit_button = QPushButton('Edit')
        self.delete_button = QPushButton('Delete')
        self.header.addWidget(self.new_button)
        self.header.addWidget(self.edit_button)
        self.header.addWidget(self.delete_button)

        self.new_button.clicked.connect(self.on_new_requested)
        self.edit_button.clicked.connect(self.on_edit_requested)
        self.delete_button.clicked.connect(self.on_delete_requested)

        self.header.addStretch()
        self.tab_layout.addLayout(self.header)

    def on_new_requested(self) -> None:
        raise NotImplementedError('New requested action must be implemented in child classes.')

    def on_edit_requested(self) -> None:
        raise NotImplementedError('Edit requested action must be implemented in child classes.')

    def on_delete_requested(self) -> None:
        raise NotImplementedError('Delete requested action must be implemented in child classes.')


class BaseTreeTab(BaseTab, Generic[ModelT]):
    def __init__(self):
        super().__init__()
        self._expansion_status: dict[ModelT, bool] = {}
        self._init_tree()

    def _init_tree(self) -> None:
        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.tab_layout.addWidget(self.tree)

        self.tree.itemExpanded.connect(self.on_item_expanded)
        self.tree.itemCollapsed.connect(self.on_item_collapsed)

    def refresh(self, selected: ModelT | None = None) -> None:
        raise NotImplementedError('Refresh action must be implemented in child classes.')

    def on_item_expanded(self, item: QTreeWidgetItem):
        if instance := item.data(0, Qt.ItemDataRole.UserRole):
            self._expansion_status[instance] = True

    def on_item_collapsed(self, item: QTreeWidgetItem):
        if instance := item.data(0, Qt.ItemDataRole.UserRole):
            self._expansion_status[instance] = False

    def is_expanded(self, instance: ModelT) -> bool:
        return self._expansion_status.get(instance, False)

    def get_selected_instance(self) -> ModelT | None:
        if item := self.tree.currentItem():
            return self.get_data(item)
        return None

    def create_item(self, instance: ModelT, texts: Sequence[str], parent: QTreeWidgetItem | None = None) -> QTreeWidgetItem:
        item = QTreeWidgetItem(parent or self.tree)
        for i, text in enumerate(texts):
            item.setText(i, text)
        self.set_data(item, instance)
        item.setExpanded(self.is_expanded(instance))
        return item

    def set_selected_item(self, item: QTreeWidgetItem) -> None:
        parent = item.parent()
        while isinstance(parent, QTreeWidgetItem):
            parent.setExpanded(True)
            parent = parent.parent()
        self.tree.setCurrentItem(item)
        self.tree.scrollToItem(item)

    def set_data(self, item: QTreeWidgetItem, instance: ModelT) -> None:
        item.setData(0, Qt.ItemDataRole.UserRole, instance)

    def get_data(self, item: QTreeWidgetItem) -> ModelT | None:
        return item.data(0, Qt.ItemDataRole.UserRole)


class BaseSortableTreeTab(BaseTreeTab[ModelT]):
    def __init__(self):
        super().__init__()
        self._init_sort_buttons()

    def _init_sort_buttons(self) -> None:
        self.up_btn = QPushButton('Move up')
        self.down_btn = QPushButton('Move down')

        self.header.addSpacerItem(QSpacerItem(10, 10))
        self.header.addWidget(self.up_btn)
        self.header.addWidget(self.down_btn)

        self.up_btn.clicked.connect(self.on_move_up_requested)
        self.down_btn.clicked.connect(self.on_move_down_requested)

    def get_sorted_list(self, instance: ModelT) -> list[ModelT]:
        raise NotImplementedError('Get sorted list action must be implemented in child classes.')

    def on_move_up_requested(self) -> None:
        self.move_selected_instance(True)

    def on_move_down_requested(self) -> None:
        self.move_selected_instance(False)

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
        self.refresh(instance)

    def update_order(self, sorted_list: Sequence[ModelT]) -> None:
        for index, instance in enumerate(sorted_list):
            instance.order = index  # type: ignore
            instance.save()
