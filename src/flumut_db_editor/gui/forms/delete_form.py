from peewee import Model
from PySide6.QtWidgets import (
    QLabel,
    QWidget,
)

from flumut_db_editor.gui.forms.base import TransactionalForm


class DeleteForm(TransactionalForm):
    def __init__(self, instance: Model, parent: QWidget | None = None) -> None:
        super().__init__(parent, instance)

        instance.delete_instance()

    def init_ui(self) -> None:
        super().init_ui()
        self.setMinimumSize(480, 420)

        self.form_layout.addWidget(QLabel(f'Are you sure you want to delete {self.instance}?'))
