from PySide6.QtWidgets import (
    QDialog,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)


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

    def __init__(self, parent, item_type: str, item_name: str, violations: dict[str, list]):
        super().__init__(parent)
        self.setWindowTitle('Cannot Delete')
        self.setGeometry(100, 100, 500, 400)

        layout = QVBoxLayout(self)

        message = QLabel(f'Cannot delete {item_type} "{item_name}" because the following items depend on it:')
        layout.addWidget(message)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)

        for dep_type, items in violations.items():
            content_layout.addWidget(QLabel(f'{dep_type}:'))  # Bold header
            for item in items:
                item_name_str = getattr(item, 'name', str(item.id))
                content_layout.addWidget(QLabel(f'  - {item_name_str}'))

        content_layout.addStretch()
        scroll.setWidget(content_widget)
        layout.addWidget(scroll)

        instruction = QLabel('Please delete these items first, then try deleting this item again.')
        layout.addWidget(instruction)

        ok_button = QPushButton('OK')
        ok_button.clicked.connect(self.accept)
        layout.addWidget(ok_button)

    @staticmethod
    def show_violation(parent, item_type: str, item_name: str, violations: dict[str, list]):
        dialog = ForeignKeyViolationDialog(parent, item_type, item_name, violations)
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
        msg_box.setIcon(QMessageBox.Critical)
        msg_box.exec()


class SuccessNotification:
    """Simple success notification."""

    @staticmethod
    def show_success(parent, message: str):
        msg_box = QMessageBox(parent)
        msg_box.setWindowTitle('Success')
        msg_box.setText(message)
        msg_box.setIcon(QMessageBox.Information)
        msg_box.exec()


class ValidationErrorDialog:
    """Dialog for form validation errors."""

    @staticmethod
    def show_validation_error(parent, field: str, message: str):
        msg_box = QMessageBox(parent)
        msg_box.setWindowTitle('Validation Error')
        msg_box.setText(f'{field}: {message}')
        msg_box.setIcon(QMessageBox.Warning)
        msg_box.exec()
