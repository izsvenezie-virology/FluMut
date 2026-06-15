from contextlib import ExitStack
from typing import ClassVar

from peewee import Model
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QWidget

from flumut.core.globals import DATABASE_PROXY


class BaseForm(QDialog):
    """Base class for modal form dialogs without persistence.

    Provides the dialog chrome (title, a :attr:`form_layout` for input widgets,
    and an OK/Cancel button box) and gates closing on :meth:`validate`.
    Suitable for non-database dialogs such as a settings form. Subclasses
    override :meth:`init_ui` (calling ``super().init_ui()`` first) to add their
    widgets, override :meth:`validate` to reject invalid input, and override
    :meth:`on_accept` to act on a valid submission.
    """

    def __init__(self, parent: QWidget | None = None, title: str = 'Form') -> None:
        super().__init__(parent)
        self.form_layout = QVBoxLayout()
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.setWindowTitle(title)
        self.resize(400, 300)
        self.init_ui()

    def init_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.addLayout(self.form_layout)
        main_layout.addWidget(self.buttons)

        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

    def accept(self) -> None:
        if self.validate():
            self.on_accept()
            super().accept()

    def validate(self) -> bool:
        """Override in subclass. Return True if the form input is valid."""
        return True

    def on_accept(self) -> None:
        """Override in subclass to act on a valid submission before closing."""


class TransactionalForm(BaseForm):
    """A :class:`BaseForm` that persists a Peewee record on accept.

    Writes go to the DB as the user works but are only committed on a valid
    submission and rolled back on Cancel / close. Wrap individual writes in
    ``with self._db.savepoint():`` so that a single failed write does not
    invalidate the whole transaction.

    Single-model forms only need to set :attr:`model` and implement
    :meth:`field_values`; :meth:`save_to_db` then creates or updates
    :attr:`instance` automatically. Forms that touch several models override
    :meth:`save_to_db` directly.
    """

    model: ClassVar[type[Model]]  # subclasses persisting a single model set this
    submitted = Signal()

    def __init__(self, parent: QWidget | None = None, instance: Model | None = None) -> None:
        self.instance = instance
        self._db = DATABASE_PROXY
        self._exit_stack = ExitStack()

        title = f'Edit {instance}' if instance else f'New {self.model.__name__}'
        super().__init__(parent, title)

    def on_accept(self) -> None:
        self.save_to_db()
        self._commit()
        self.submitted.emit()

    def field_values(self) -> dict:
        """Override in subclass. Map model field names to current widget values."""
        return {}

    def create_values(self) -> dict:
        """Values used to create a new record. Defaults to :meth:`field_values`;
        override to add create-only fields such as a parent foreign key."""
        return self.field_values()

    def save_to_db(self) -> None:
        """Create or update :attr:`instance` from the form values.

        Override directly in forms that persist more than one model.
        """
        if self.instance is None:
            self.instance = self.model.create(**self.create_values())
        else:
            for field, value in self.field_values().items():
                setattr(self.instance, field, value)
            self.instance.save()

    def _begin_transaction(self) -> None:
        self._exit_stack.enter_context(self._db.manual_commit())
        self._db.begin()

    def _commit(self) -> None:
        if self._db.in_transaction():
            self._db.commit()
        self._exit_stack.close()

    def _rollback(self) -> None:
        if self._db.in_transaction():
            self._db.rollback()
        self._exit_stack.close()

    def reject(self) -> None:
        self._rollback()
        super().reject()
