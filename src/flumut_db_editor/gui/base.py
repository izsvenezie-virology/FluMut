from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QPushButton,
    QTableWidget,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)

from flumut_db_editor.gui.crud_mixin import CrudMixin, TableCrudMixin, TreeCrudMixin


class BaseTab(QWidget, CrudMixin):
    """Common tab chrome: a New/Refresh header above an item view.

    Subclasses implement :meth:`create_view` to supply the table or tree; the
    CRUD behaviour (context menu, refresh, delete) comes from the mixin.
    """

    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.addLayout(self._build_header())

        view = self.create_view()
        view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        view.customContextMenuRequested.connect(self.show_context_menu)
        layout.addWidget(view)

    def _build_header(self) -> QHBoxLayout:
        header = QHBoxLayout()
        self.new_button = QPushButton('New')
        self.refresh_button = QPushButton('Refresh')
        self.refresh_button.clicked.connect(self.on_refresh)
        header.addWidget(self.new_button)
        header.addWidget(self.refresh_button)
        header.addStretch()
        return header

    def create_view(self) -> QAbstractItemView:
        """Create and return the tab's item view - to be overridden."""
        raise NotImplementedError


class BaseTableTab(BaseTab, TableCrudMixin):
    def create_view(self) -> QAbstractItemView:
        self.table = QTableWidget()
        self.table.setColumnCount(0)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        return self.table


class HierarchicalTab(BaseTab, TreeCrudMixin):
    def create_view(self) -> QAbstractItemView:
        self.tree = QTreeWidget()
        self.tree.setColumnCount(1)
        self.tree.setHeaderLabel('Items')
        return self.tree
