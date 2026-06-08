from PySide6.QtWidgets import QHBoxLayout, QHeaderView, QPushButton, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget

from flumut_db_editor.models import Effect, Host, Subtype


class EffectsTab(QWidget):
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

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.effects_table = self._create_table(['ID', 'Name'])
        self.subtypes_table = self._create_table(['ID', 'Name'])
        self.hosts_table = self._create_table(['ID', 'Name'])

        self.tabs.addTab(self.effects_table, 'Effects')
        self.tabs.addTab(self.subtypes_table, 'Subtypes')
        self.tabs.addTab(self.hosts_table, 'Hosts')

        self.load_data()

    def _create_table(self, headers):
        table = QTableWidget()
        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        return table

    def load_data(self):
        self._load_effects()
        self._load_subtypes()
        self._load_hosts()

    def _load_effects(self):
        effects = Effect.select()
        self.effects_table.setRowCount(len(effects))
        for row, effect in enumerate(effects):
            self.effects_table.setItem(row, 0, QTableWidgetItem(str(effect.id)))
            self.effects_table.setItem(row, 1, QTableWidgetItem(effect.name))

    def _load_subtypes(self):
        subtypes = Subtype.select()
        self.subtypes_table.setRowCount(len(subtypes))
        for row, subtype in enumerate(subtypes):
            self.subtypes_table.setItem(row, 0, QTableWidgetItem(str(subtype.id)))
            self.subtypes_table.setItem(row, 1, QTableWidgetItem(subtype.name))

    def _load_hosts(self):
        hosts = Host.select()
        self.hosts_table.setRowCount(len(hosts))
        for row, host in enumerate(hosts):
            self.hosts_table.setItem(row, 0, QTableWidgetItem(str(host.id)))
            self.hosts_table.setItem(row, 1, QTableWidgetItem(host.name))
