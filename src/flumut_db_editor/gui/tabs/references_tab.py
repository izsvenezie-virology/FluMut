from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.gui.dialogs import SuccessNotification
from flumut_db_editor.gui.forms.annotation_form import AnnotationForm
from flumut_db_editor.gui.forms.reference_form import ReferenceForm
from flumut_db_editor.gui.tabs.base import HierarchicalTab
from flumut_db_editor.models import Annotation, Reference


class ReferencesTab(HierarchicalTab):
    def __init__(self):
        super().__init__()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.new_button.clicked.connect(self.handle_new)
        self.load_data()

    def load_data(self):
        self.tree.clear()
        references: list[Reference] = Reference.select()
        for ref in references:
            ref_item = QTreeWidgetItem(self.tree)
            ref_item.setText(0, ref.name)
            ref_item.setText(1, f'Segment: {ref.segment.name} | Source: {ref.source}')
            ref_item.setData(0, 0x0100, ('reference', ref.id))

            annotations: list[Annotation] = Annotation.select().where(Annotation.reference == ref)
            for annot in annotations:
                annot_item = QTreeWidgetItem(ref_item)
                annot_item.setText(
                    0,
                    f'{annot.protein.name} ({annot.start}-{annot.end})',
                )
                annot_item.setText(1, annot.notes or '')
                annot_item.setData(0, 0x0100, ('annotation', annot.id))

    def handle_new(self):
        form = ReferenceForm(self, None)
        if form.exec():
            SuccessNotification.show_success(self, 'Reference created successfully.')
            self.load_data()

    def handle_edit(self):
        item = self.get_selected_item()
        if item is None:
            return
        item_type, item_id = self.get_item_data(item)
        if item_type == 'reference':
            reference = Reference.get_by_id(item_id)
            form = ReferenceForm(self, reference)
            if form.exec():
                SuccessNotification.show_success(self, 'Reference updated successfully.')
                self.load_data()
        elif item_type == 'annotation':
            annotation = Annotation.get_by_id(item_id)
            form = AnnotationForm(self, annotation)
            if form.exec():
                SuccessNotification.show_success(self, 'Annotation updated successfully.')
                self.load_data()

    def handle_delete(self):
        item = self.get_selected_item()
        if item is None:
            return
        item_type, item_id = self.get_item_data(item)
        if item_type == 'reference':
            if delete_with_confirmation(self, 'Reference', item_id, DatabaseOperations.delete_reference):
                self.load_data()
        elif item_type == 'annotation':
            if delete_with_confirmation(self, 'Annotation', item_id, DatabaseOperations.delete_annotation):
                self.load_data()
