from flumut.flumutdb.models import Paper
from flumut_db_editor.gui.forms.delete_form import DeleteForm
from flumut_db_editor.gui.forms.paper_form import PaperForm
from flumut_db_editor.gui.tabs.base import BaseTableTab


class PapersTab(BaseTableTab[Paper]):
    def __init__(self):
        super().__init__()
        self.__init_ui()

    def __init_ui(self) -> None:
        headers = ['Short Name', 'Authors', 'Year', 'Title', 'Journal', 'DOI', 'URL', 'Notes']
        self.table.setColumnCount(len(headers))
        self.table.setHorizontalHeaderLabels(headers)

        self.refresh()

    def populate(self, selected=None):
        papers: list[Paper] = sorted(Paper.select(), key=lambda p: p.short_name)
        self.table.clearContents()
        self.table.setRowCount(len(papers))

        rows = {}
        for paper in papers:
            rows[paper] = [
                paper.short_name,
                paper.authors,
                str(paper.year),
                paper.title,
                paper.journal or '',
                paper.doi or '',
                paper.url or '',
                paper.notes or '',
            ]
        self.populate_table(rows, selected)  # pyright: ignore[reportArgumentType]

    def on_new_requested(self):
        form = PaperForm(self, None)
        if form.exec():
            self.refresh(form.instance)

    def on_edit_requested(self):
        instance = self.get_selected_item()
        if instance is None:
            return
        form = PaperForm(self, instance)
        if form.exec():
            self.refresh(form.instance)

    def on_delete_requested(self):
        instance = self.get_selected_item()
        if instance and DeleteForm.confirm_and_delete(instance, self):
            self.refresh()
