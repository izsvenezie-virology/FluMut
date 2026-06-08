from PySide6.QtWidgets import QLabel, QLineEdit, QSpinBox

from flumut_db_editor.gui.forms.base import BaseForm


class PaperForm(BaseForm):
    def __init__(self, parent=None, paper=None):
        self.paper = paper
        super().__init__(parent, 'Paper')

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

        if self.paper:
            self.short_name_field.setText(self.paper.short_name)
            self.title_field.setText(self.paper.title)
            self.authors_field.setText(self.paper.authors)
            self.year_field.setValue(self.paper.year)

    def get_data(self):
        return {
            'short_name': self.short_name_field.text(),
            'title': self.title_field.text(),
            'authors': self.authors_field.text(),
            'year': self.year_field.value(),
        }
