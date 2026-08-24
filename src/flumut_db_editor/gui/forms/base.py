from collections.abc import Iterable
from typing import Generic, TypeVar

from peewee import DatabaseError
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QWidget

from flumut.core.globals import DATABASE_PROXY
from flumut.flumutdb.models import BaseModel
from flumut_db_editor.gui.dialogs import DataErrorDialog, ErrorDialog
from flumut_db_editor.validator import DataValidator

ModelT = TypeVar('ModelT', bound=BaseModel)
RelatedT = TypeVar('RelatedT', bound=BaseModel)


class BaseForm(QDialog):
    """Dialog shell that validates and saves any collection of model instances in one transaction.

    Subclasses decide which instances take part through the hooks below:
    `form_title`, `populate_instances`, `instances_to_save` and, when rows can
    be removed, `instances_to_delete`. The instances may belong to different
    model classes.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

    def _init_ui(self) -> None:
        self.resize(400, 300)
        self.setWindowTitle(self.form_title())

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
        self.populate_instances()
        for validator in self.validators():
            validator.validate()
            if validator.errors:
                DataErrorDialog(self, validator).exec()
                return False
        target: BaseModel | None = None
        try:
            with DATABASE_PROXY.atomic():
                for target in self.instances_to_save():
                    target.save()
                for target in self.instances_to_delete():
                    target.delete_instance()
        except DatabaseError as error:
            name = type(target).__name__ if target is not None else 'data'
            ErrorDialog.show_error(self, 'Save failed', f'Could not save {name}.', str(error))
            return False
        return True

    def validators(self) -> Iterable[DataValidator]:
        return (DataValidator(instance) for instance in self.instances_to_save())

    def form_title(self) -> str:
        raise NotImplementedError('form_title must be implemented in child classes.')

    def instances_to_save(self) -> Iterable[BaseModel]:
        raise NotImplementedError('instances_to_save must be implemented in child classes.')

    def instances_to_delete(self) -> Iterable[BaseModel]:
        return ()

    def populate_instances(self) -> None:
        """Sync widget values into the instances before validation. Override when widgets do not edit the instances directly."""

    def field_values(self) -> dict:
        """Override in subclass. Map model field names to current widget values."""
        raise NotImplementedError('field_values must be implemented in child classes.')


class TransactionalForm(BaseForm, Generic[ModelT]):
    """Edits a single instance of `model`."""

    model: type[ModelT]

    def __init__(self, parent: QWidget | None = None, instance: ModelT | None = None) -> None:
        super().__init__(parent)
        self._init_state(instance)
        self._init_ui()

    def _init_state(self, instance: ModelT | None) -> None:
        self.instance: ModelT = self.model.get_by_id(instance.get_id()) if instance is not None else self.model()
        self.validator = DataValidator(self.instance)

    def form_title(self) -> str:
        return f'New {self.model.__name__}' if self.instance.get_id() is None else f'Edit {self.instance}'

    def validators(self) -> Iterable[DataValidator]:
        return (self.validator,)

    def instances_to_save(self) -> Iterable[BaseModel]:
        return (self.instance,)

    def populate_instances(self) -> None:
        self.populate_instance()

    def populate_instance(self) -> None:
        for field, value in self.field_values().items():
            setattr(self.instance, field, value)


class MasterDetailForm(TransactionalForm[ModelT], Generic[ModelT, RelatedT]):
    """Edits a single `model` instance together with related instances of another class.

    `self.related` holds the related instances currently in the form; move an
    instance into `self.removed_related` to delete it on save. The master is
    saved first, so related rows can reference it even when it is new.
    """

    def _init_state(self, instance: ModelT | None) -> None:
        super()._init_state(instance)
        self.related: list[RelatedT] = list(self.load_related())
        self.removed_related: list[RelatedT] = []

    def load_related(self) -> Iterable[RelatedT]:
        raise NotImplementedError('load_related must be implemented in child classes.')

    def instances_to_save(self) -> Iterable[BaseModel]:
        return (self.instance, *self.related)

    def populate_instances(self) -> None:
        self.populate_instance()
        self.populate_related()

    def populate_related(self) -> None:
        raise NotImplementedError('populate_related must be implemented in child classes.')


class MultiInstanceForm(BaseForm, Generic[ModelT]):
    """Edits several instances of the same `model` at once, applying `field_values` to each."""

    model: type[ModelT]

    def __init__(self, parent: QWidget | None = None, instances: Iterable[ModelT] = ()) -> None:
        super().__init__(parent)
        self.instances: list[ModelT] = [
            self.model.get_by_id(instance.get_id()) if instance.get_id() is not None else instance for instance in instances
        ]
        self._init_ui()

    def form_title(self) -> str:
        return f'Edit {len(self.instances)} {self.model.__name__} instances'

    def instances_to_save(self) -> Iterable[BaseModel]:
        return self.instances

    def populate_instances(self) -> None:
        values = self.field_values()
        for instance in self.instances:
            for field, value in values.items():
                setattr(instance, field, value)
