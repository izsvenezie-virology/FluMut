from collections.abc import Iterable
from typing import Generic, TypeVar

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb.models import BaseModel

ModelT = TypeVar('ModelT', bound=BaseModel)


class FilterableList(QWidget, Generic[ModelT]):
    """A titled list of instances narrowed down by a filter box.

    Items not matching the filter are hidden and deselected, so
    :meth:`selected_instances` only ever returns items the user can see.
    """

    def __init__(self, title: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.__init_ui(title)

    def __init_ui(self, title: str) -> None:
        self.filter_field = QLineEdit()
        self.filter_field.setPlaceholderText('Filter...')
        self.list = QListWidget()
        self.list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.new_btn = QPushButton('New')

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        title_row = QHBoxLayout()
        title_row.addWidget(QLabel(title))
        title_row.addStretch()
        title_row.addWidget(self.new_btn)

        layout.addLayout(title_row)
        layout.addWidget(self.filter_field)
        layout.addWidget(self.list)

        self.filter_field.textChanged.connect(self.apply_filter)

    def set_instances(self, instances: Iterable[ModelT]) -> None:
        self.list.clear()
        for instance in instances:
            item = QListWidgetItem(str(instance))
            item.setData(Qt.ItemDataRole.UserRole, instance)
            self.list.addItem(item)
        self.apply_filter(self.filter_field.text())

    def selected_instances(self) -> list[ModelT]:
        return [item.data(Qt.ItemDataRole.UserRole) for item in self.list.selectedItems()]

    def select_instances(self, instances: Iterable[ModelT]) -> None:
        """Add `instances` to the selection, skipping the ones the filter is hiding."""
        wanted = list(instances)
        for row in range(self.list.count()):
            item = self.list.item(row)
            if item.isHidden() or item.data(Qt.ItemDataRole.UserRole) not in wanted:
                continue
            item.setSelected(True)
            self.list.scrollToItem(item)

    def apply_filter(self, text: str) -> None:
        needle = text.strip().lower()
        for row in range(self.list.count()):
            item = self.list.item(row)
            hidden = needle not in item.text().lower()
            item.setHidden(hidden)
            if hidden:
                item.setSelected(False)
