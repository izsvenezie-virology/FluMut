from PySide6.QtWidgets import QLineEdit


class FilterLineEdit(QLineEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setPlaceholderText('Filter...')
