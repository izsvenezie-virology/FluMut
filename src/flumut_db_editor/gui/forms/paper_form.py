from PySide6.QtWidgets import QLabel, QLineEdit, QSpinBox

from flumut.flumutdb.models import Paper
from flumut_db_editor.gui.dialogs import ValidationErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm


class PaperForm(TransactionalForm):
    model = Paper

    def __init__(self, parent=None, instance=None):
        super().__init__(parent, 'Paper', instance)

    def init_ui(self):
        super().init_ui()
        self.short_name_field = QLineEdit()
        self.title_field = QLineEdit()
        self.authors_field = QLineEdit()
        self.year_field = QSpinBox()
        self.year_field.setMinimum(1900)
        self.year_field.setMaximum(2100)

        self.form_layout.insertWidget(0, QLabel('Short Name:'))
        self.form_layout.insertWidget(1, self.short_name_field)
        self.form_layout.insertWidget(2, QLabel('Title:'))
        self.form_layout.insertWidget(3, self.title_field)
        self.form_layout.insertWidget(4, QLabel('Authors:'))
        self.form_layout.insertWidget(5, self.authors_field)
        self.form_layout.insertWidget(6, QLabel('Year:'))
        self.form_layout.insertWidget(7, self.year_field)

        if self.instance:
            self.short_name_field.setText(self.instance.short_name)
            self.title_field.setText(self.instance.title)
            self.authors_field.setText(self.instance.authors)
            self.year_field.setValue(self.instance.year)

    def validate(self) -> bool:
        short_name = self.short_name_field.text().strip()
        title = self.title_field.text().strip()

        if not short_name:
            ValidationErrorDialog.show_validation_error(self, 'Short Name', 'Short name cannot be empty.')
            return False
        if not title:
            ValidationErrorDialog.show_validation_error(self, 'Title', 'Title cannot be empty.')
            return False
        if not self.instance:
            if Paper.select().where(Paper.short_name == short_name).exists():
                ValidationErrorDialog.show_validation_error(self, 'Short Name', 'A paper with this short name already exists.')
                return False
        return True

    def field_values(self) -> dict:
        return {
            'short_name': self.short_name_field.text().strip(),
            'title': self.title_field.text().strip(),
            'authors': self.authors_field.text().strip(),
            'year': self.year_field.value(),
        }
