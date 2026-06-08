from PySide6.QtCore import Signal
from PySide6.QtWidgets import QDialog, QHBoxLayout, QPushButton, QVBoxLayout


class BaseForm(QDialog):
    submitted = Signal()

    def __init__(self, parent=None, title='Form'):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setGeometry(100, 100, 400, 300)
        self.init_ui()

    def init_ui(self):
        self.layout = QVBoxLayout(self)
        self.form_layout = QVBoxLayout()
        self.layout.addLayout(self.form_layout)

        button_layout = QHBoxLayout()
        self.ok_button = QPushButton('OK')
        self.cancel_button = QPushButton('Cancel')
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        self.layout.addLayout(button_layout)

        self.ok_button.clicked.connect(self.submit)
        self.cancel_button.clicked.connect(self.reject)

    def submit(self):
        self.submitted.emit()
        self.accept()
