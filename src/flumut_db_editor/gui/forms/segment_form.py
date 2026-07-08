from typing import TYPE_CHECKING

from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QWidget,
)

from flumut.flumutdb.models import Segment
from flumut_db_editor.gui.forms.base import TransactionalForm


class SegmentForm(TransactionalForm):
    model = Segment

    def __init__(self, parent: QWidget | None = None, instance: Segment | None = None) -> None:
        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: Segment | None

    def init_ui(self) -> None:
        super().init_ui()
        self.resize(480, 85)

        name_layout = QHBoxLayout()

        self.name_field = QLineEdit()
        name_layout.addWidget(QLabel('Name:'))
        name_layout.addWidget(self.name_field)

        self.form_layout.addLayout(name_layout)
        self.form_layout.addStretch()

        if self.instance:
            self.name_field.setText(self.instance.name)
        self.name_field.setFocus()

    @property
    def name(self) -> str:
        return self.name_field.text().strip()

    def validate(self) -> bool:
        return self.check_unique_required('name', self.name, 'Name', self.name_field)

    def field_values(self) -> dict:
        return {'name': self.name}
