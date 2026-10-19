from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QPushButton, QSpinBox, QWidget

from flumut.flumutdb.models import Paper
from flumut_db_editor.gui.dialogs import ErrorDialog
from flumut_db_editor.gui.forms.base import TransactionalForm
from flumut_db_editor.papers import CrossrefError, build_short_name, disambiguate_short_name, fetch_metadata


class PaperForm(TransactionalForm[Paper]):
    model = Paper

    def __init__(self, parent: QWidget | None = None, instance: Paper | None = None) -> None:
        super().__init__(parent, instance)
        self.__init_ui()

    def __init_ui(self):
        self.short_name_field = QLineEdit()
        self.title_field = QLineEdit()
        self.authors_field = QLineEdit()
        self.year_field = QSpinBox()
        self.journal_field = QLineEdit()
        self.url_field = QLineEdit()
        self.doi_field = QLineEdit()

        self.year_field.setMinimum(1900)
        self.year_field.setMaximum(2100)

        self.generate_button = QPushButton('Generate')
        self.generate_button.setToolTip('Build the short name from the first author and the year.')
        self.generate_button.clicked.connect(self.on_generate_requested)
        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.setToolTip('Fill the form with the metadata Crossref holds for this DOI.')
        self.fetch_button.clicked.connect(self.on_fetch_requested)

        short_name_row = QHBoxLayout()
        short_name_row.addWidget(QLabel('Short name'))
        short_name_row.addWidget(self.short_name_field)
        short_name_row.addWidget(self.generate_button)
        title_row = QHBoxLayout()
        title_row.addWidget(QLabel('Title'))
        title_row.addWidget(self.title_field)
        authors_row = QHBoxLayout()
        authors_row.addWidget(QLabel('Authors'))
        authors_row.addWidget(self.authors_field)
        year_row = QHBoxLayout()
        year_row.addWidget(QLabel('Year'))
        year_row.addWidget(self.year_field)
        journal_row = QHBoxLayout()
        journal_row.addWidget(QLabel('Juornal'))
        journal_row.addWidget(self.journal_field)
        url_row = QHBoxLayout()
        url_row.addWidget(QLabel('URL'))
        url_row.addWidget(self.url_field)
        doi_row = QHBoxLayout()
        doi_row.addWidget(QLabel('DOI'))
        doi_row.addWidget(self.doi_field)
        doi_row.addWidget(self.fetch_button)

        self.form_layout.addLayout(short_name_row)
        self.form_layout.addLayout(title_row)
        self.form_layout.addLayout(authors_row)
        self.form_layout.addLayout(year_row)
        self.form_layout.addLayout(journal_row)
        self.form_layout.addLayout(url_row)
        self.form_layout.addLayout(doi_row)

        if self.instance:
            self.short_name_field.setText(self.instance.short_name)
            self.title_field.setText(self.instance.title)
            self.authors_field.setText(self.instance.authors)
            self.year_field.setValue(self.instance.year or 1900)
            self.journal_field.setText(self.instance.journal)
            self.doi_field.setText(self.instance.doi)
            self.url_field.setText(self.instance.url)

    def on_fetch_requested(self) -> None:
        self.fetch_button.setEnabled(False)
        QGuiApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            metadata = fetch_metadata(self.doi_field.text())
        except CrossrefError as error:
            ErrorDialog.show_error(self, 'Lookup failed', 'Could not fetch the metadata of this DOI.', str(error))
            return
        finally:
            QGuiApplication.restoreOverrideCursor()
            self.fetch_button.setEnabled(True)

        self.doi_field.setText(metadata.doi)
        self.title_field.setText(metadata.title)
        self.authors_field.setText(metadata.authors)
        if metadata.year:
            self.year_field.setValue(metadata.year)
        self.journal_field.setText(metadata.journal)
        if not self.url_field.text().strip():
            self.url_field.setText(metadata.url)
        if not self.short_name_field.text().strip():
            self.on_generate_requested()

    def on_generate_requested(self) -> None:
        short_name = build_short_name(self.authors_field.text().strip(), self.year_field.value())
        if not short_name:
            ErrorDialog.show_error(self, 'Cannot generate', 'Fill in the authors and the year first.')
            return
        self.short_name_field.setText(disambiguate_short_name(short_name, self.instance))

    @property
    def short_name(self) -> str:
        return self.short_name_field.text().strip()

    @property
    def title(self) -> str:
        return self.title_field.text().strip()

    @property
    def authors(self) -> str:
        return self.authors_field.text().strip()

    @property
    def year(self) -> int:
        return self.year_field.value()

    @property
    def journal(self) -> str | None:
        return self.journal_field.text().strip()

    @property
    def url(self) -> str | None:
        return self.url_field.text().strip()

    @property
    def doi(self) -> str | None:
        return self.doi_field.text().strip()

    def field_values(self) -> dict:
        return {
            'short_name': self.short_name,
            'title': self.title,
            'authors': self.authors,
            'year': self.year,
            'journal': self.journal or None,
            'url': self.url or None,
            'doi': self.doi or None,
        }
