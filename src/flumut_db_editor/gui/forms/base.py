from contextlib import ExitStack

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QDialog, QHBoxLayout, QPushButton, QVBoxLayout

from flumut.core.globals import DATABASE_PROXY


class BaseForm(QDialog):
    submitted = Signal()

    def __init__(self, parent=None, title='Form'):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setGeometry(100, 100, 400, 300)
        self.init_ui()

    def init_ui(self):
        self.setLayout(QVBoxLayout(self))
        self.form_layout = QVBoxLayout(self.layout())
        self.addWidget(self.form_layout)

        button_layout = QHBoxLayout()
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')
        button_layout.addWidget(self.cancel_button)
        button_layout.addWidget(self.ok_button)
        self.layout.addLayout(button_layout)

        self.ok_button.clicked.connect(self.submit)
        self.cancel_button.clicked.connect(self.reject)

    def submit(self):
        if self.validate():
            self.save_to_db()
            self.submitted.emit()
            self.accept()

    def validate(self) -> bool:
        """Override in subclass. Return True if form is valid."""
        return True

    def save_to_db(self) -> None:
        """Override in subclass. Save form data to database."""
        pass


class TransactionalForm(BaseForm):
    """BaseForm with an explicit SQLite transaction: all writes go to the DB
    immediately but are only committed when the user confirms (OK), or fully
    rolled back on Cancel / close.  Individual operations inside the form
    should use self._db.savepoint() so that a single failed write does not
    invalidate the whole transaction.
    """

    def __init__(self, parent=None, title='Form'):
        self._db = DATABASE_PROXY
        self._exit_stack = ExitStack()
        super().__init__(parent, title)

    def _begin_transaction(self):
        self._exit_stack.enter_context(self._db.manual_commit())
        self._db.begin()

    def _commit(self):
        if self._db.in_transaction():
            self._db.commit()
        self._exit_stack.close()

    def _rollback(self):
        if self._db.in_transaction():
            self._db.rollback()
        self._exit_stack.close()

    def reject(self):
        self._rollback()
        super().reject()
