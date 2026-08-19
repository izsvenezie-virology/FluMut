from PySide6.QtWidgets import (
    QDialog,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from flumut_db_editor.validator import DataValidator, DeleteValidator


class DataErrorDialog(QDialog):
    def __init__(self, parent: QWidget | None, validator: DataValidator) -> None:
        super().__init__(parent)
        self.setWindowTitle('Data validation errors')
        self.setGeometry(100, 100, 500, 400)

        layout = QVBoxLayout(self)

        message = QLabel(f'Cannot save {type(validator.instance)} "{validator.instance}" because of the following errors:')
        layout.addWidget(message)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)

        for attr, items in validator.errors.items():
            for msg in items:
                content_layout.addWidget(QLabel(f'{attr}: {msg}'))

        content_layout.addStretch()
        scroll.setWidget(content_widget)
        layout.addWidget(scroll)

        instruction = QLabel('Please, solve these errors, then try to save again.')
        layout.addWidget(instruction)

        ok_button = QPushButton('OK')
        ok_button.clicked.connect(self.accept)
        layout.addWidget(ok_button)


class ForeignKeyViolationDialog(QDialog):
    """Dialog showing FK violations preventing deletion."""

    def __init__(self, parent: QWidget | None, validator: DeleteValidator):
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
