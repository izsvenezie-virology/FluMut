from typing import TYPE_CHECKING

from PySide6.QtWidgets import (
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from flumut.flumutdb.models import BaseModel
from flumut_db_editor.delete_validator import DeleteValidator
from flumut_db_editor.gui.dialogs import ForeignKeyViolationDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class DeleteForm(TransactionalForm):
    def __init__(self, instance: BaseModel, validator: DeleteValidator, parent: QWidget | None = None) -> None:
        self.validator = validator  # set before super() so init_ui() can use it
        super().__init__(parent, instance)
        if TYPE_CHECKING:
            self.instance: BaseModel

    def init_ui(self) -> None:
        super().init_ui()
        self.setMinimumSize(480, 420)

        self.form_layout.addWidget(QLabel(f'Are you sure you want to delete {self.instance}?'))

        if self.validator.cascade_items:
            message = QLabel(
                f'Delete {type(self.validator.instance)} "{self.validator.instance}" will delete the following items that depend on it:'
            )
            self.form_layout.addWidget(message)

            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            content_widget = QWidget()
            content_layout = QVBoxLayout(content_widget)

            for dep_type, items in self.validator.cascade_items.items():
                content_layout.addWidget(QLabel(f'{dep_type} ({len(items)}):'))
                for item in items:
                    content_layout.addWidget(QLabel(f'  - {item}'))

            content_layout.addStretch()
            scroll.setWidget(content_widget)
            self.form_layout.addWidget(scroll)

    def save_to_db(self) -> None:
        """Persist the form by deleting the instance (committed by ``on_accept``)."""
        self.instance.delete_instance()

    @staticmethod
    def confirm_and_delete(instance: BaseModel, parent: QWidget | None = None) -> bool:
        """Validate, confirm, and delete ``instance``. Returns True if deleted."""
        validator = DeleteValidator(instance)
        if not validator.can_delete():
            ForeignKeyViolationDialog.show_violation(parent, validator)
            return False
        return bool(DeleteForm(instance, validator, parent).exec())
