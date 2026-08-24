from collections.abc import Sequence

from PySide6.QtWidgets import QPushButton

from flumut.flumutdb.models import Protein, Segment
from flumut_db_editor.gui.forms.annotation_form import AnnotationForm
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.reference_form import ReferenceForm
from flumut_db_editor.gui.tabs.base import BaseSortableTreeTab
from flumut_db_editor.models import Annotation, Reference


class ReferencesTab(BaseSortableTreeTab[Segment | Reference | Protein | Annotation]):
    def __init__(self):
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        self.new_annotation_btn = QPushButton('New annotation')
        self.new_annotation_btn.clicked.connect(self.on_new_annotation_requested)

        self.new_btn.setText('New reference')  # type: ignore
        self.header.insertWidget(1, self.new_annotation_btn)

        self.refresh()

    def refresh(self, selected=None) -> None:
        self.tree.clear()
        segments: Sequence[Segment] = sorted(Segment.select())
        for segment in segments:
            segment_item = self.create_item(
                segment, [f'Segment: {segment.name}', f'{len(segment.proteins)} proteins, {len(segment.references)} references']
            )
            if segment not in self._expansion_status:
                self._expansion_status[segment] = True
            if segment == selected:
                self.set_selected_item(segment_item)

            references: Sequence[Reference] = sorted(segment.references)
            for reference in references:
                message = self.get_reference_description_text(reference)
                reference_item = self.create_item(reference, [f'Reference: {reference.name}', message], segment_item)
                if reference == selected:
                    self.set_selected_item(reference_item)

                for protein in sorted({a.protein for a in reference.annotations}):
                    annotations = self.get_annotations(reference, protein)
                    if protein not in self._expansion_status:
                        self._expansion_status[protein] = True
                    protein_item = self.create_item(
                        protein, [f'Protein: {protein.name}', f'{len(annotations)} annotations'], reference_item
                    )
                    if protein == selected:
                        self.set_selected_item(protein_item)

                    for annotation in annotations:
                        annotation_item = self.create_item(
                            annotation, [f'{annotation.order} - {annotation.start}-{annotation.end}'], protein_item
                        )
                        if annotation == selected:
                            self.set_selected_item(annotation_item)

    def on_new_requested(self):
        segment = self.get_selected_segment()
        form = ReferenceForm(self, None, segment)
        if form.exec() and form.instance:
            lst = self.get_sorted_list(form.instance)

            if selected := self.get_selected_reference():
                lst.remove(form.instance)
                lst.insert(lst.index(selected) + 1, form.instance)

            self.update_order(lst)
            self.refresh(form.instance)

    def on_new_annotation_requested(self):
        reference = self.get_selected_reference()
        if not reference:
            return
        protein = self.get_selected_protein()
        form = AnnotationForm(reference, self, None, protein)
        if form.exec() and form.instance:
            lst = self.get_sorted_list(form.instance)

            if selected := self.get_selected_annotation():
                lst.remove(form.instance)
                lst.insert(lst.index(selected), form.instance)

            self.update_order(lst)
            self.refresh(form.instance)

    def on_edit_requested(self):
        instance = self.get_selected_instance()
        if (form := self.get_form(instance)) and form.exec() and form.instance:
            self.refresh()

    def on_delete_requested(self):
        instance = self.get_selected_instance()
        if instance and DeleteForm.confirm_and_delete(instance, self):
            list_to_order = self.get_sorted_list(instance)
            self.update_order(list_to_order)
            self.refresh()

    def get_reference_description_text(self, reference: Reference) -> str:
        missing_annotated_proteins = {a.protein for a in reference.annotations} - set(reference.segment.proteins)
        message = f'{len(reference.sequence)} bp, {len(reference.mappings)} mapped mutations'
        if missing_annotated_proteins:
            message = f'⚠ Missing annotation for protein(s) {missing_annotated_proteins} - {message}'
        return message

    def get_selected_segment(self) -> Segment | None:
        selected = self.get_selected_instance()
        match selected:
            case Segment():
                return selected
            case Reference():
                return selected.segment
            case Protein():
                return selected.segment
            case Annotation():
                return selected.reference.segment
            case _:
                return None

    def get_selected_reference(self) -> Reference | None:
        selected = self.get_selected_instance()
        match selected:
            case Reference():
                return selected
            case Annotation():
                return selected.reference
            case Protein():
                reference_item = self.tree.currentItem().parent()
                return self.get_data(reference_item)  # type: ignore
            case _:
                return None

    def get_selected_protein(self) -> Protein | None:
        selected = self.get_selected_instance()
        match selected:
            case Protein():
                return selected
            case Annotation():
                return selected.protein
            case _:
                return None

    def get_selected_annotation(self) -> Annotation | None:
        selected = self.get_selected_instance()
        match selected:
            case Annotation():
                return selected
            case _:
                return None

    def get_sorted_list(self, instance) -> list[Reference | Annotation]:
        match instance:
            case Reference():
                return sorted(instance.segment.references)
            case Annotation():
                return sorted(self.get_annotations(instance.reference, instance.protein))
            case _:
                return []

    def get_form(self, instance: Segment | Reference | Protein | Annotation | None) -> ReferenceForm | AnnotationForm | None:
        match instance:
            case Reference():
                return ReferenceForm(self, instance)
            case Annotation():
                return AnnotationForm(instance.reference, self, instance, instance.protein)
            case _:
                return None

    def get_annotations(self, reference: Reference, protein: Protein) -> Sequence[Annotation]:
        return [ann for ann in reference.annotations if ann.protein == protein]
