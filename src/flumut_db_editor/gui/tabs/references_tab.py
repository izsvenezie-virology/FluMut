from PySide6.QtCore import Qt
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

from flumut_db_editor.gui.forms.annotation_form import AnnotationForm
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.reference_form import ReferenceForm
from flumut_db_editor.gui.tabs.base import BaseTab
from flumut_db_editor.models import Annotation, Reference


class ReferencesTab(BaseTab):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self) -> None:
        self.tree = QTreeWidget()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.tab_layout.addWidget(self.tree)
        self.refresh()

    def refresh(self):
        self.tree.clear()
        references: list[Reference] = Reference.select()
        for ref in sorted(references, key=(lambda x: x.name)):
            ref_item = QTreeWidgetItem(self.tree)
            ref_item.setText(0, ref.name)
            proteins_count = len(set([mapping.mutation.protein for mapping in ref.mappings]))
            ref_item.setText(
                1,
                f'Segment: {ref.segment.name} | Source: {ref.source} | Sequence {len(ref.sequence)} | Mappings: {len(ref.mappings)} over {proteins_count} protein{"" if proteins_count == 1 else "s"}',
            )
            ref_item.setData(0, Qt.ItemDataRole.UserRole, ref)

            for annot in ref.annotations:
                annot_item = QTreeWidgetItem(ref_item)
                annot_item.setText(0, f'{annot.protein.name}')
                annot_item.setText(1, f'{annot.start}-{annot.end}')
                annot_item.setData(0, Qt.ItemDataRole.UserRole, annot)

    def on_new_requested(self):
        form = ReferenceForm(self, None)
        if form.exec():
            self.refresh()

    def on_edit_requested(self):
        instance = self.get_selected_instance()
        if isinstance(instance, Reference):
            form = ReferenceForm(self, instance)
        elif isinstance(instance, Annotation):
            form = AnnotationForm(self, instance)
        else:
            return

        if form.exec():
            self.refresh()

    def on_delete_requested(self):
        instance = self.get_selected_instance()
        if instance:
            if DeleteForm.confirm_and_delete(instance, self):
                self.refresh()

    def get_selected_instance(self) -> Reference | Annotation | None:
        if item := self.tree.currentItem():
            return item.data(0, Qt.ItemDataRole.UserRole)
        return None
