from collections.abc import Iterable, Iterator

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem, QPushButton, QVBoxLayout

from flumut.flumutdb import loader
from flumut.flumutdb.models import BaseModel, Marker, MarkerMutation, Mutation
from flumut_db_editor.gui.forms.base import TransactionalForm


class MarkerForm(TransactionalForm[Marker]):
    model = Marker

    def __init__(self, parent=None, instance=None):
        super().__init__(parent, instance)
        self.initial_links: list[MarkerMutation] = list(self.load_links())
        self.selected_mutations: list[Mutation] = [link.mutation for link in self.initial_links]
        self.__init_ui()

    def __init_ui(self):
        name_row = QHBoxLayout()
        self.name_field = QLineEdit()
        self.generate_name_btn = QPushButton('Generate')

        name_row.addWidget(QLabel('Name:'))
        name_row.addWidget(self.name_field)
        name_row.addWidget(self.generate_name_btn)
        self.form_layout.addLayout(name_row)

        self.selectable_mutations_list = QListWidget()
        self.selected_mutations_list = QListWidget()

        self.add_mutation_btn = QPushButton('>')
        self.remove_mutation_btn = QPushButton('<')
        mutation_btn_column = QVBoxLayout()
        mutation_btn_column.addStretch()
        mutation_btn_column.addWidget(self.add_mutation_btn)
        mutation_btn_column.addWidget(self.remove_mutation_btn)
        mutation_btn_column.addStretch()

        self.form_layout.addWidget(QLabel('Mutations'))
        mutations_row = QHBoxLayout()
        mutations_row.addWidget(self.selectable_mutations_list)
        mutations_row.addLayout(mutation_btn_column)
        mutations_row.addWidget(self.selected_mutations_list)
        self.form_layout.addLayout(mutations_row)

        self.add_mutation_btn.clicked.connect(self.add_mutation)
        self.remove_mutation_btn.clicked.connect(self.remove_mutation)
        self.selectable_mutations_list.doubleClicked.connect(self.add_mutation)
        self.selected_mutations_list.doubleClicked.connect(self.remove_mutation)
        self.generate_name_btn.clicked.connect(self.on_generate_name_requested)

        if self.instance:
            self.name_field.setText(self.instance.name)

        self.refresh_lists()

    def add_mutation(self) -> None:
        item = self.selectable_mutations_list.currentItem()
        if item:
            self.selected_mutations.append(self.get_data(item))
        self.refresh_lists()

    def remove_mutation(self) -> None:
        item = self.selected_mutations_list.currentItem()
        if item:
            self.selected_mutations.remove(self.get_data(item))
        self.refresh_lists()

    def refresh_lists(self) -> None:
        selection = []
        if item := self.selectable_mutations_list.currentItem():
            selection.append(self.get_data(item))
        if item := self.selected_mutations_list.currentItem():
            selection.append(self.get_data(item))

        self.selected_mutations_list.clear()
        self.selectable_mutations_list.clear()

        for mut in self.selected_mutations:
            item = self.create_item(mut)
            self.selected_mutations_list.addItem(item)
            if mut in selection:
                self.selected_mutations_list.setCurrentItem(item)

        mutations = sorted(loader.get(Mutation))
        for mut in mutations:
            if mut not in self.selected_mutations:
                item = self.create_item(mut)
                self.selectable_mutations_list.addItem(item)
                if mut in selection:
                    self.selected_mutations_list.setCurrentItem(item)

    def on_generate_name_requested(self) -> None:
        self.name_field.setText(self.generate_name())

    def generate_name(self) -> str:
        return '; '.join(sorted(mutation.name for mutation in self.selected_mutations))

    def create_item(self, mutation: Mutation) -> QListWidgetItem:
        item = QListWidgetItem(mutation.name)
        item.setData(Qt.ItemDataRole.UserRole, mutation)
        return item

    def get_data(self, item: QListWidgetItem) -> Mutation:
        return item.data(Qt.ItemDataRole.UserRole)

    def load_links(self) -> Iterable[MarkerMutation]:
        if self.instance.get_id() is None:
            return ()
        return MarkerMutation.select(MarkerMutation, Mutation).join(Mutation).where(MarkerMutation.marker == self.instance)

    def field_values(self) -> dict:
        return {'name': self.name_field.text().strip()}

    def instances_to_save(self) -> Iterator[BaseModel]:
        yield self.instance
        linked_mutations = [link.mutation for link in self.initial_links]
        for mutation in self.selected_mutations:
            if mutation not in linked_mutations:
                yield MarkerMutation(marker=self.instance, mutation=mutation)

    def instances_to_delete(self) -> Iterable[BaseModel]:
        return [link for link in self.initial_links if link.mutation not in self.selected_mutations]
