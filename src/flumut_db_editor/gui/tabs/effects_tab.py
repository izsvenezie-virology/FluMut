from PySide6.QtWidgets import QTableWidgetItem

from flumut.flumutdb.models import Effect
from flumut_db_editor.gui.base import BaseTableTab
from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.effect_form import EffectForm


class EffectsTab(BaseTableTab):
    def __init__(self):
        super().__init__()
        self.new_button.clicked.connect(self.handle_new)
        self.load_data()

    def load_data(self):
        header = ['Name', 'Notes']
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        effects: list[Effect] = Effect.select()
        self.table.setRowCount(len(effects))
        for row, effect in enumerate(effects):
            self.table.setItem(row, 0, QTableWidgetItem(effect.name))
            self.table.setItem(row, 1, QTableWidgetItem(effect.notes or ''))
        self.table.resizeColumnsToContents()

    def handle_new(self):
        form = EffectForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Effect created successfully.')
            self.load_data()

    def handle_edit(self):
        row = self.get_selected_item()
        if row is None:
            return
        effects = list(Effect.select())
        effect = effects[row]
        form = EffectForm(self, effect)
        if form.exec():
            SuccessNotification.show_success(self, 'Effect updated successfully.')
            self.load_data()

    def handle_delete(self):
        row = self.get_selected_item()
        if row is None:
            return
        effects = list(Effect.select())
        effect = effects[row]
        if delete_with_confirmation(self, 'Effect', effect.id, DatabaseOperations.delete_effect):
            self.load_data()
