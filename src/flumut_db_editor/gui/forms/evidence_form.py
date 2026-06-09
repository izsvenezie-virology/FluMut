from PySide6.QtWidgets import QComboBox, QLabel, QTextEdit

from flumut.flumutdb.models import Effect, Evidence, Host, Marker, Paper, Subtype
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import BaseForm


class EvidenceForm(BaseForm):
    def __init__(self, parent=None, evidence=None):
        self.evidence = evidence
        super().__init__(parent, 'Evidence')

    def init_ui(self):
        super().init_ui()
        self.marker_combo = QComboBox()
        self.paper_combo = QComboBox()
        self.effect_combo = QComboBox()
        self.subtype_combo = QComboBox()
        self.host_combo = QComboBox()
        self.notes_field = QTextEdit()

        markers = Marker.select()
        for marker in markers:
            self.marker_combo.addItem(marker.name or str(marker.id), marker.id)

        papers = Paper.select()
        for paper in papers:
            self.paper_combo.addItem(paper.short_name, paper.id)

        effects = Effect.select()
        for effect in effects:
            self.effect_combo.addItem(effect.name, effect.id)

        subtypes = Subtype.select()
        for subtype in subtypes:
            self.subtype_combo.addItem(subtype.name, subtype.id)

        self.host_combo.addItem('(None)', None)
        hosts = Host.select()
        for host in hosts:
            self.host_combo.addItem(host.name, host.id)

        self.form_layout.insertWidget(0, QLabel('Marker:'))
        self.form_layout.insertWidget(1, self.marker_combo)
        self.form_layout.insertWidget(2, QLabel('Paper:'))
        self.form_layout.insertWidget(3, self.paper_combo)
        self.form_layout.insertWidget(4, QLabel('Effect:'))
        self.form_layout.insertWidget(5, self.effect_combo)
        self.form_layout.insertWidget(6, QLabel('Subtype:'))
        self.form_layout.insertWidget(7, self.subtype_combo)
        self.form_layout.insertWidget(8, QLabel('Host:'))
        self.form_layout.insertWidget(9, self.host_combo)
        self.form_layout.insertWidget(10, QLabel('Notes:'))
        self.form_layout.insertWidget(11, self.notes_field)

        if self.evidence:
            self.marker_combo.setCurrentIndex(self.marker_combo.findData(self.evidence.marker_id))
            self.paper_combo.setCurrentIndex(self.paper_combo.findData(self.evidence.paper_id))
            self.effect_combo.setCurrentIndex(self.effect_combo.findData(self.evidence.effect_id))
            self.subtype_combo.setCurrentIndex(self.subtype_combo.findData(self.evidence.subtype_id))
            self.host_combo.setCurrentIndex(self.host_combo.findData(self.evidence.host_id))
            self.notes_field.setText(self.evidence.notes or '')

    def validate(self) -> bool:
        if self.marker_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Marker', 'Please select a marker.')
            return False
        if self.paper_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Paper', 'Please select a paper.')
            return False
        if self.effect_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Effect', 'Please select an effect.')
            return False
        if self.subtype_combo.count() == 0:
            ValidationErrorDialog.show_validation_error(self, 'Subtype', 'Please select a subtype.')
            return False
        return True

    def save_to_db(self) -> None:
        marker_id = self.marker_combo.currentData()
        paper_id = self.paper_combo.currentData()
        effect_id = self.effect_combo.currentData()
        subtype_id = self.subtype_combo.currentData()
        host_id = self.host_combo.currentData()
        notes = self.notes_field.toPlainText().strip() or None

        if self.evidence:
            self.evidence.marker_id = marker_id
            self.evidence.paper_id = paper_id
            self.evidence.effect_id = effect_id
            self.evidence.subtype_id = subtype_id
            self.evidence.host_id = host_id
            self.evidence.notes = notes
            self.evidence.save()
        else:
            Evidence.create(
                marker_id=marker_id,
                paper_id=paper_id,
                effect_id=effect_id,
                subtype_id=subtype_id,
                host_id=host_id,
                notes=notes,
            )
