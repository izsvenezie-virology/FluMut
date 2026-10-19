from collections.abc import Iterable

from PySide6.QtWidgets import QLabel, QScrollArea, QVBoxLayout, QWidget

from flumut.flumutdb.models import BaseModel
from flumut_db_editor.gui.dialogs import ForeignKeyViolationDialog
from flumut_db_editor.gui.forms.base import BaseForm
from flumut_db_editor.validator import DeleteValidator


class DeleteForm(BaseForm):
    """Confirms the deletion of one instance, then lets the base form carry it out.

    The dialog itself is the confirmation: nothing is saved, so the instance to
    remove is the whole of `instances_to_delete`. Blocking dependants are ruled
    out by `DeleteValidator` before the dialog is ever built.
    """

    def __init__(self, validator: DeleteValidator, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.validator = validator
        self.instance = validator.instance
        self.__init_ui()

    def __init_ui(self) -> None:
        self.setWindowTitle(self.form_title())
        self.form_layout.addWidget(QLabel(f'Are you sure you want to delete {self.instance}?'))
        self.form_layout.addStretch()

        if not self.validator.cascade_items:
            return

        message = QLabel(f'Delete {type(self.instance).__name__} "{self.instance}" will delete the following items that depend on it:')
        self.form_layout.addWidget(message)

        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        for dep_type, items in self.validator.cascade_items.items():
            content_layout.addWidget(QLabel(f'{dep_type} ({len(items)}):'))
            for item in items:
                content_layout.addWidget(QLabel(f'  - {item}'))
        content_layout.addStretch()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(content_widget)
        self.form_layout.addWidget(scroll)

    def form_title(self) -> str:
        return f'Deleting {type(self.instance).__name__}'

    def instances_to_save(self) -> Iterable[BaseModel]:
        return ()

    def instances_to_delete(self) -> Iterable[BaseModel]:
        return (self.instance,)

    def failure_message(self, target: BaseModel | None) -> tuple[str, str]:
        return 'Deletion failed', f'Could not delete {type(self.instance).__name__}.'

    @staticmethod
    def confirm_and_delete(instance: BaseModel, parent: QWidget | None = None) -> bool:
        """Validate, confirm, and delete ``instance``. Returns True if deleted."""
        validator = DeleteValidator(instance)
        if not validator.can_delete():
            ForeignKeyViolationDialog.show_violation(parent, validator)
            return False
        return bool(DeleteForm(validator, parent).exec())
