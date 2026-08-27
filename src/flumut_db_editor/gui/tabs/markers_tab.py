from flumut.flumutdb import loader
from flumut.flumutdb.models import Marker, Mutation
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.marker_form import MarkerForm
from flumut_db_editor.gui.tabs.base import BaseTreeTab


class MarkersTab(BaseTreeTab[Marker | Mutation]):
    def __init__(self):
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        self.new_btn.setText('New marker')
        self.refresh()

    def refresh(self, selected=None):
        self.tree.clear()
        markers = loader.get(Marker)
        for marker in markers:
            marker_item = self.create_item(marker, [f'Marker: {marker.name}', f'{len(marker.mutations)} mutations'])
            if marker == selected:
                self.set_selected_item(marker_item)

            for mutation in marker.mutations:
                mutation_item = self.create_item(mutation, [f'Mutation: {mutation.name}'], marker_item)
                if mutation == selected:
                    self.set_selected_item(mutation_item)

    def on_new_requested(self) -> None:
        form = MarkerForm(self, None)
        if form.exec():
            self.refresh(form.instance)

    def on_edit_requested(self) -> None:
        form = MarkerForm(self, self.get_selected_marker())
        if form.exec():
            self.refresh(form.instance)

    def on_delete_requested(self) -> None:
        instance = self.get_selected_instance()
        if instance and DeleteForm.confirm_and_delete(instance, self):
            self.refresh()

    def get_selected_marker(self) -> Marker | None:
        selected = self.get_selected_instance()
        match selected:
            case Marker():
                return selected
            case Mutation():
                return self.get_data(self.tree.currentItem().parent())  # pyright: ignore[reportReturnType]
            case _:
                return None
