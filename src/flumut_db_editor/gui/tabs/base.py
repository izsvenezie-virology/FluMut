from collections.abc import Iterable, Mapping, Sequence
from typing import Generic, TypeVar

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QSpacerItem,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb import loader
from flumut.flumutdb.models import BaseModel, SortableModel
from flumut_db_editor.gui.forms.base import EvidenceTermsForm
from flumut_db_editor.gui.forms.delete_form import DeleteForm

ModelT = TypeVar('ModelT', bound=BaseModel)
ModelS = TypeVar('ModelS', bound=SortableModel)


class BaseTab(QWidget):
    """Common tab chrome: a New/Refresh header above an item view.

    Subclasses implement :meth:`create_view` to supply the table or tree; the
    CRUD behaviour (context menu, refresh, delete) comes from the mixin.
    """

    def __init__(self):
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        self.tab_layout = QVBoxLayout(self)

        self.header = QHBoxLayout()
        self.new_btn = QPushButton('New')
        self.edit_btn = QPushButton('Edit')
        self.delete_btn = QPushButton('Delete')
        self.filter_field = QLineEdit()
        self.filter_field.setPlaceholderText('Filter...')
        self.filter_field.setClearButtonEnabled(True)
        self.filter_field.setMinimumWidth(200)

        self.header.addWidget(self.new_btn)
        self.header.addWidget(self.edit_btn)
        self.header.addWidget(self.delete_btn)
        self.header.addWidget(self.filter_field)

        self.new_btn.clicked.connect(self.on_new_requested)
        self.edit_btn.clicked.connect(self.on_edit_requested)
        self.delete_btn.clicked.connect(self.on_delete_requested)
        self.filter_field.textChanged.connect(self.apply_filter)

        self.header.addStretch()
        self.tab_layout.addLayout(self.header)

    def refresh(self, selected: BaseModel | None = None) -> None:
        """Rebuild the view, then hide again whatever the filter leaves out."""
        self.populate(selected)
        self.apply_filter(self.filter_field.text())

    def populate(self, selected: BaseModel | None = None) -> None:
        raise NotImplementedError('populate must be implemented in child classes.')

    def apply_filter(self, text: str) -> None:
        raise NotImplementedError('apply_filter must be implemented in child classes.')

    def on_new_requested(self) -> None:
        raise NotImplementedError('on_new_requested must be implemented in child classes.')

    def on_edit_requested(self) -> None:
        raise NotImplementedError('on_edit_requested must be implemented in child classes.')

    def on_delete_requested(self) -> None:
        raise NotImplementedError('on_delete_requested must be implemented in child classes.')


class BaseTableTab(BaseTab, Generic[ModelT]):
    def __init__(self):
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        self.table = QTableWidget(self)
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(['Name', 'Notes'])
        self.table.verticalHeader().setVisible(False)
        self.table.setColumnWidth(0, 200)
        self.tab_layout.addWidget(self.table)

    def apply_filter(self, text: str) -> None:
        needle = text.strip().lower()
        for row in range(self.table.rowCount()):
            items = (self.table.item(row, column) for column in range(self.table.columnCount()))
            self.table.setRowHidden(row, needle not in ' '.join(item.text() for item in items if item).lower())

    def get_selected_instance(self) -> ModelT | None:
        if item := self.table.currentItem():
            return self.get_data(item)
        return None

    def populate_table(self, rows: Mapping[ModelT, Iterable[str]], selected: ModelT | None = None) -> None:
        self.table.clearContents()
        self.table.setRowCount(len(rows))
        for row, (instance, texts) in enumerate(rows.items()):
            items = self.create_row(instance, texts)
            for col, item in enumerate(items):
                self.table.setItem(row, col, item)
            if instance == selected:
                self.table.selectRow(row)
                self.table.scrollTo(self.table.selectedIndexes()[0])

    def create_row(self, instance: ModelT, texts: Iterable[str]) -> list[QTableWidgetItem]:
        row = []
        for text in texts:
            item = QTableWidgetItem(text)
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            row.append(item)
        self.set_data(row[0], instance)
        return row

    def set_selected_item(self, item: QTableWidgetItem) -> None:
        self.table.setCurrentItem(item)
        self.table.scrollToItem(item)

    def get_selected_item(self) -> ModelT | None:
        row = self.table.currentRow()
        return self.get_data(self.table.item(row, 0))

    def get_selected_items(self) -> list[ModelT]:
        """Every selected instance. Only meaningful while the table selects whole rows."""
        rows = self.table.selectionModel().selectedRows()
        return [instance for index in rows if (instance := self.get_data(self.table.item(index.row(), 0)))]

    def set_data(self, item: QTableWidgetItem, instance: ModelT) -> None:
        item.setData(Qt.ItemDataRole.UserRole, instance)

    def get_data(self, item: QTableWidgetItem | None) -> ModelT | None:
        if item is None:
            return None
        return item.data(Qt.ItemDataRole.UserRole)


class EvidenceTermsTab(BaseTableTab[ModelT]):
    model: type[ModelT]
    form: type[EvidenceTermsForm]

    def __init__(self):
        super().__init__()
        self.refresh()

    def populate(self, selected=None):
        terms: list[ModelT] = list(self.model.select())
        rows = {}
        for term in terms:
            rows[term] = [term.name, term.notes or '']  # pyright: ignore[reportAttributeAccessIssue]
        self.populate_table(rows, selected)  # pyright: ignore[reportArgumentType]
        self.table.resizeColumnsToContents()

    def on_new_requested(self):
        form = self.form(self, None)
        if form.exec():
            self.refresh(form.instance)

    def on_edit_requested(self):
        instance = self.get_selected_item()
        if instance is None:
            return
        form = self.form(self, instance)
        if form.exec():
            self.refresh(form.instance)

    def on_delete_requested(self):
        instance = self.get_selected_item()
        if instance and DeleteForm.confirm_and_delete(instance, self):
            self.refresh()


class BaseTreeTab(BaseTab, Generic[ModelT]):
    def __init__(self):
        super().__init__()
        self._expansion_status: dict[ModelT, bool] = {}
        self.__init_ui()

    def __init_ui(self) -> None:
        self.tree = QTreeWidget(self)
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.tree.setColumnWidth(0, 200)
        self.tab_layout.addWidget(self.tree)

        self.tree.itemExpanded.connect(self.on_item_expanded)
        self.tree.itemCollapsed.connect(self.on_item_collapsed)

    def apply_filter(self, text: str) -> None:
        needle = text.strip().lower()
        for row in range(self.tree.topLevelItemCount()):
            self.filter_item(self.tree.topLevelItem(row), needle)  # pyright: ignore[reportArgumentType]

    def filter_item(self, item: QTreeWidgetItem, needle: str) -> bool:
        """Hide `item` unless it or one of its children matches `needle`. Returns whether it stayed visible."""
        children = [self.filter_item(item.child(row), needle) for row in range(item.childCount())]
        texts = ' '.join(item.text(column) for column in range(self.tree.columnCount())).lower()
        matches = any(children) or needle in texts
        item.setHidden(not matches)
        if needle and any(children):
            item.setExpanded(True)
        return matches

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

    def get_selected_instances(self) -> list[ModelT]:
        return [instance for item in self.tree.selectedItems() if (instance := self.get_data(item))]

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


class BaseSortableTreeTab(BaseTreeTab[ModelS]):
    def __init__(self):
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        self.up_btn = QPushButton('Move up')
        self.down_btn = QPushButton('Move down')

        self.header.addSpacerItem(QSpacerItem(10, 10))
        self.header.addWidget(self.up_btn)
        self.header.addWidget(self.down_btn)

        self.up_btn.clicked.connect(self.on_move_up_requested)
        self.down_btn.clicked.connect(self.on_move_down_requested)

    def get_sorted_list(self, instance: ModelS) -> Sequence[ModelS]:
        raise NotImplementedError('Get sorted list action must be implemented in child classes.')

    def on_move_up_requested(self) -> None:
        self.move_selected_instance(True)

    def on_move_down_requested(self) -> None:
        self.move_selected_instance(False)

    def move_selected_instance(self, up: bool) -> None:
        instance = self.get_selected_instance()
        if not instance:
            return
        list_to_move = list(self.get_sorted_list(instance))
        if not list_to_move:
            return

        idx = list_to_move.index(instance)
        new_idx = idx + (-1 if up else 1)
        if not 0 <= new_idx <= len(list_to_move):
            return

        list_to_move.insert(new_idx, list_to_move.pop(idx))
        self.update_order(list_to_move)
        self.refresh(instance)

    def update_order(self, sorted_list: Sequence[ModelS]) -> None:
        for index, instance in enumerate(sorted_list):
            instance.order = index + 1
            instance.save()
        loader.load()
