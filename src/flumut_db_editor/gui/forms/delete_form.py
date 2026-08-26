from peewee import DatabaseError
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from flumut.core.globals import DATABASE_PROXY
from flumut.flumutdb import loader
from flumut.flumutdb.models import BaseModel
from flumut_db_editor.gui.dialogs import ErrorDialog, ForeignKeyViolationDialog
from flumut_db_editor.validator import DeleteValidator


class DeleteForm(QDialog):
    def __init__(self, validator: DeleteValidator, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.validator = validator
        self.instance = validator.instance

        self._init_ui()

    def _init_ui(self) -> None:
        self.resize(400, 300)
        self.setWindowTitle(f'Deleting {type(self.instance).__name__}')

        self.form_layout = QVBoxLayout()
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)

        main_layout = QVBoxLayout(self)
        main_layout.addLayout(self.form_layout)
        main_layout.addWidget(self.buttons)

        self.form_layout.addWidget(QLabel(f'Are you sure you want to delete {self.instance}?'))

        if self.validator.cascade_items:
            message = QLabel(
                f'Delete {type(self.instance).__name__} "{self.validator.instance}" will delete the following items that depend on it:'
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

        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

    def accept(self) -> None:
        if self.save_to_db():
            super().accept()

    def save_to_db(self) -> bool:
        try:
            with DATABASE_PROXY.atomic():
                self.instance.delete_instance()
        except DatabaseError as error:
            ErrorDialog.show_error(self, 'Deletion failed', f'Could not delete {type(self.instance).__name__}.', str(error))
            return False
        loader.load()
        return True

    @staticmethod
    def confirm_and_delete(instance: BaseModel, parent: QWidget | None = None) -> bool:
        """Validate, confirm, and delete ``instance``. Returns True if deleted."""
        validator = DeleteValidator(instance)
        if not validator.can_delete():
            ForeignKeyViolationDialog.show_violation(parent, validator)
            return False
        return bool(DeleteForm(validator, parent).exec())
