from PySide6.QtWidgets import QTableWidgetItem

from flumut.flumutdb.models import Host
from flumut_db_editor.gui.base import BaseTableTab


class HostsTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        header = ['Name', 'Notes']
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        hosts: list[Host] = Host.select()
        self.table.setRowCount(len(hosts))
        for row, host in enumerate(hosts):
            self.table.setItem(row, 0, QTableWidgetItem(host.name))
            self.table.setItem(row, 1, QTableWidgetItem(host.notes or ''))
        self.table.resizeColumnsToContents()
