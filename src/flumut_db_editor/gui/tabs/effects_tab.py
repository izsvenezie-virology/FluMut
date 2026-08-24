from flumut.flumutdb.models import Effect
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.effect_form import EffectForm
from flumut_db_editor.gui.tabs.base import BaseTableTab


class EffectsTab(BaseTableTab[Effect]):
    def __init__(self):
        super().__init__()
        self.load_data()

    def load_data(self):
        effects: list[Effect] = list(Effect.select())
        rows = {}
        for effect in effects:
            rows[effect] = [effect.name, effect.notes or '']
        self.populate_table(rows)
        self.table.resizeColumnsToContents()

    def on_new_requested(self):
        form = EffectForm(self, None)
        if form.exec():
            self.load_data()

    def on_edit_requested(self):
        instance = self.get_selected_item()
        if instance is None:
            return
        form = EffectForm(self, instance)
        if form.exec():
            self.load_data()

    def on_delete_requested(self):
        instance = self.get_selected_item()
        if instance and DeleteForm.confirm_and_delete(instance, self):
            self.load_data()
