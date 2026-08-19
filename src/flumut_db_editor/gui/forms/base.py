from typing import Generic, TypeVar

from peewee import DatabaseError
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QWidget

from flumut.core.globals import DATABASE_PROXY
from flumut.flumutdb.models import BaseModel
from flumut_db_editor.gui.dialogs import DataErrorDialog, ErrorDialog
from flumut_db_editor.validator import DataValidator

ModelT = TypeVar('ModelT', bound=BaseModel)
RelatedT = TypeVar('RelatedT', bound=BaseModel | None)


class TransactionalForm(QDialog, Generic[ModelT]):
    model: type[ModelT]

    def __init__(self, parent: QWidget | None = None, instance: ModelT | None = None) -> None:
        super().__init__(parent)
        self.instance: ModelT = self.model.get_by_id(instance.get_id()) if instance is not None else self.model()
        self.validator = DataValidator(self.instance)
        self._init_ui()

    def _init_ui(self) -> None:
        self.resize(400, 300)
        self.setWindowTitle(f'New {self.model.__name__}' if self.instance.get_id() is None else f'Edit {self.instance}')

        self.form_layout = QVBoxLayout()
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)

        main_layout = QVBoxLayout(self)
        main_layout.addLayout(self.form_layout)
        main_layout.addWidget(self.buttons)

        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

    def accept(self) -> None:
        if not self.save_to_db():
            return
        super().accept()

    def save_to_db(self) -> bool:
        self.populate_instance()
        self.validator.validate()
        if self.validator.errors:
            DataErrorDialog(self, self.validator).exec()
            return False
        try:
            with DATABASE_PROXY.atomic():
                self.instance.save()
        except DatabaseError as error:
            ErrorDialog.show_error(self, 'Save failed', f'Could not save {self.model.__name__}.', str(error))
            return False
        return True

    def field_values(self) -> dict:
        """Override in subclass. Map model field names to current widget values."""
        raise NotImplementedError('field_values must be implemented in child classes.')

    def populate_instance(self) -> None:
        for field, value in self.field_values().items():
            setattr(self.instance, field, value)
