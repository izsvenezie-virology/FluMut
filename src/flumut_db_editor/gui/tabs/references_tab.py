from PySide6.QtWidgets import QTreeWidgetItem

from flumut_db_editor.gui.base import HierarchicalTab
from flumut_db_editor.models import Annotation, Reference


class ReferencesTab(HierarchicalTab):
    def __init__(self):
        super().__init__()
        self.tree.setColumnCount(2)
        self.tree.setHeaderLabels(['Name', 'Details'])
        self.load_data()

    def load_data(self):
        self.tree.clear()
        references: list[Reference] = Reference.select()
        for ref in references:
            ref_item = QTreeWidgetItem(self.tree)
            ref_item.setText(0, ref.name)
            ref_item.setText(1, f'Segment: {ref.segment.name} | Source: {ref.source}')
            ref_item.setData(0, 0x0100, ('reference', ref.id))

            annotations: list[Annotation] = Annotation.select().where(
                Annotation.reference == ref
            )
            for annot in annotations:
                annot_item = QTreeWidgetItem(ref_item)
                annot_item.setText(
                    0,
                    f'{annot.protein.name} ({annot.start}-{annot.end})',
                )
                annot_item.setText(1, annot.notes or '')
                annot_item.setData(0, 0x0100, ('annotation', annot.id))
