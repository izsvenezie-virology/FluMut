from PySide6.QtWidgets import (
    QDialog,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from flumut_db_editor.delete_validator import DeleteValidator


class ConfirmationDialog(QDialog):
    """Confirmation dialog for delete operations."""

    def __init__(self, parent, title: str, message: str, details: str = ''):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setGeometry(100, 100, 400, 200)
        self.result_value = False

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(message))
        if details:
            layout.addWidget(QLabel(details))

        button_layout = QVBoxLayout()
        ok_button = QPushButton('Yes, Delete')
        cancel_button = QPushButton('Cancel')
        ok_button.clicked.connect(self.confirm)
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(ok_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)

    def confirm(self):
        self.result_value = True
        self.accept()

    @staticmethod
    def ask(parent, title: str, message: str, details: str = '') -> bool:
        dialog = ConfirmationDialog(parent, title, message, details)
        return dialog.exec() == QDialog.Accepted and dialog.result_value


class ForeignKeyViolationDialog(QDialog):
    """Dialog showing FK violations preventing deletion."""

    def __init__(self, parent, validator: DeleteValidator):
        super().__init__(parent)
        self.setWindowTitle('Cannot Delete')
        self.setGeometry(100, 100, 500, 400)

        layout = QVBoxLayout(self)

        message = QLabel(f'Cannot delete {type(validator.instance)} "{validator.instance}" because the following items depend on it:')
        layout.addWidget(message)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)

        for dep_type, items in validator.blocking_items.items():
            content_layout.addWidget(QLabel(f'{dep_type} ({len(items)}):'))
            for item in items:
                content_layout.addWidget(QLabel(f'  - {item}'))

        content_layout.addStretch()
        scroll.setWidget(content_widget)
        layout.addWidget(scroll)

        instruction = QLabel('Please delete these items first, then delete this item again.')
        layout.addWidget(instruction)

        ok_button = QPushButton('OK')
        ok_button.clicked.connect(self.accept)
        layout.addWidget(ok_button)

    @staticmethod
    def show_violation(parent, validator: DeleteValidator):
        dialog = ForeignKeyViolationDialog(parent, validator)
        dialog.exec()


class ErrorDialog:
    """Simple error message display."""

    @staticmethod
    def show_error(parent, title: str, message: str, details: str = ''):
        msg_box = QMessageBox(parent)
        msg_box.setWindowTitle(title)
        msg_box.setText(message)
        if details:
            msg_box.setDetailedText(details)
        msg_box.setIcon(QMessageBox.Icon.Critical)
        msg_box.exec()


class SuccessNotification:
    """Simple success notification."""

    @staticmethod
    def show_success(parent, message: str):
        msg_box = QMessageBox(parent)
        msg_box.setWindowTitle('Success')
        msg_box.setText(message)
        msg_box.setIcon(QMessageBox.Icon.Information)
        msg_box.exec()


class ValidationErrorDialog:
    """Dialog for form validation errors."""

    @staticmethod
    def show_validation_error(parent, field: str, message: str):
        msg_box = QMessageBox(parent)
        msg_box.setWindowTitle(f'Validation Error: {field}')
        msg_box.setText(f'{message}')
        msg_box.setIcon(QMessageBox.Icon.Warning)
        msg_box.exec()
