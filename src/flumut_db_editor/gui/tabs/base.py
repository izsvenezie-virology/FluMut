from typing import Generic, TypeVar

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb.models import BaseModel
from flumut_db_editor.gui.crud_mixin import TableCrudMixin, TreeCrudMixin

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
        self.header.addWidget(self.new_button)
        self.header.addWidget(self.edit_button)

        self.new_button.clicked.connect(self.on_new_requested)
        self.edit_button.clicked.connect(self.on_edit_requested)

        self.header.addStretch()
        self.tab_layout.addLayout(self.header)

    def on_new_requested(self) -> None:
        raise NotImplementedError('New requested action must be implemented in child classes.')

    def on_edit_requested(self) -> None:
        raise NotImplementedError('Edit requested action must be implemented in child classes.')

    def on_delete_requested(self) -> None:
        raise NotImplementedError('Delete requested action must be implemented in child classes.')


class BaseTableTab(BaseTab, TableCrudMixin):
    def __init__(self):
        super().__init__()


class HierarchicalTab(BaseTab, TreeCrudMixin):
    def create_view(self) -> QAbstractItemView:
        self.tree = QTreeWidget()
        self.tree.setColumnCount(1)
        self.tree.setHeaderLabel('Items')
        return self.tree


class BaseTreeTab(BaseTab, Generic[ModelT]):
    def __init__(self):
        super().__init__()

        self._expansion_status: dict[ModelT, bool] = {}

        self._init_ui()

    def _init_ui(self) -> None:
        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.tab_layout.addWidget(self.tree)

        self.tree.itemExpanded.connect(self.on_item_expanded)
        self.tree.itemCollapsed.connect(self.on_item_collapsed)

        self.refresh()

    def refresh(self) -> None:
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

    def set_data(self, item: QTreeWidgetItem, instance: ModelT) -> None:
        item.setData(0, Qt.ItemDataRole.UserRole, instance)

    def get_data(self, item: QTreeWidgetItem) -> ModelT | None:
        return item.data(0, Qt.ItemDataRole.UserRole)
