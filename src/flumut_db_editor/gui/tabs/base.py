from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QPushButton,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)

from flumut_db_editor.gui.crud_mixin import TableCrudMixin, TreeCrudMixin


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
