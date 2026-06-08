from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QTableWidget,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)


class BaseTableTab(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        header_layout = QHBoxLayout()
        self.new_button = QPushButton('New')
        self.refresh_button = QPushButton('Refresh')
        header_layout.addWidget(self.new_button)
        header_layout.addWidget(self.refresh_button)
        header_layout.addStretch()
        layout.addLayout(header_layout)

        self.table = QTableWidget()
        self.table.setColumnCount(0)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.show_context_menu)
        layout.addWidget(self.table)

    def show_context_menu(self, pos):
        pass


class HierarchicalTab(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)

        header_layout = QHBoxLayout()
        self.new_button = QPushButton('New')
        self.refresh_button = QPushButton('Refresh')
        header_layout.addWidget(self.new_button)
        header_layout.addWidget(self.refresh_button)
        header_layout.addStretch()
        layout.addLayout(header_layout)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(1)
        self.tree.setHeaderLabel('Items')
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self.show_context_menu)
        layout.addWidget(self.tree)

    def show_context_menu(self, pos):
        pass
