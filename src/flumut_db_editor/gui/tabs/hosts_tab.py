from PySide6.QtWidgets import QTableWidgetItem

from flumut.flumutdb.models import Host
from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.host_form import HostForm


class HostsTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.new_button.clicked.connect(self.handle_new)
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

    def handle_new(self):
        form = HostForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Host created successfully.')
            self.load_data()

    def handle_edit(self):
        row = self.get_selected_item()
        if row is None:
            return
        hosts = list(Host.select())
        host = hosts[row]
        form = HostForm(self, host)
        if form.exec():
            SuccessNotification.show_success(self, 'Host updated successfully.')
            self.load_data()

    def handle_delete(self):
        row = self.get_selected_item()
        if row is None:
            return
        hosts = list(Host.select())
        host = hosts[row]
        if delete_with_confirmation(self, 'Host', host.id, DatabaseOperations.delete_host):
            self.load_data()
