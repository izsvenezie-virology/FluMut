"""The editor's front end for proposing a FluMutDB release.

All the work is done by :mod:`flumut_db_editor.release`; this only shows the
curator what it worked out, takes their summary, and reports what happened.
"""

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QApplication,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QStackedWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from flumut_db_editor.gui.dialogs import ErrorDialog, SuccessNotification
from flumut_db_editor.release import Proposal, ReleaseError, prepare, pull_request_url, release

LEVEL_EXPLANATIONS = {
    'patch': 'Analysis data changed: a cached analysis has to be re-analysed.',
    'minor': 'Preprocessing data changed: a cached analysis has to be preprocessed again.',
    'major': 'The schema changed: this is a breaking database release. Update DB_MAJOR_VERSION in src/flumut/core/globals.py by hand.',
}


class ReleaseTab(QWidget):
    """Shows the release the edited database calls for, and carries it out once the curator accepts."""

    def __init__(self) -> None:
        super().__init__()
        self.proposal: Proposal | None = None
        # The summary last filled in for the curator, so a refresh can tell it apart from their own.
        self.offered_summary = ''

        self.pages = QStackedWidget()
        self.problem_page = self._problem_page()
        self.proposal_page = self._proposal_page()
        self.pages.addWidget(self.problem_page)
        self.pages.addWidget(self.proposal_page)

        layout = QVBoxLayout(self)
        layout.addWidget(self.pages)

    def _problem_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        self.problem_label = QLabel()
        self.problem_label.setWordWrap(True)
        self.problem_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        check_again = QPushButton('Check again')
        check_again.clicked.connect(self.refresh)
        buttons = QHBoxLayout()
        buttons.addWidget(check_again)
        buttons.addStretch()

        layout.addWidget(QLabel('<b>Cannot propose a release</b>'))
        layout.addWidget(self.problem_label)
        layout.addLayout(buttons)
        layout.addStretch()
        return page

    def _proposal_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        versions = QGroupBox('Release')
        form = QFormLayout(versions)
        self.level_label = QLabel()
        self.level_label.setWordWrap(True)
        self.db_version_label = QLabel()
        self.flumut_version_label = QLabel()
        self.branch_label = QLabel()
        form.addRow('Bump:', self.level_label)
        form.addRow('FluMutDB:', self.db_version_label)
        form.addRow('FluMut:', self.flumut_version_label)
        form.addRow('Branch:', self.branch_label)

        changes = QGroupBox('What changed')
        self.changes_tree = QTreeWidget()
        self.changes_tree.setHeaderLabels(['Table', 'Rows changed'])
        self.changes_tree.setRootIsDecorated(False)
        QVBoxLayout(changes).addWidget(self.changes_tree)

        summary = QGroupBox('Summary')
        self.summary_field = QPlainTextEdit()
        self.summary_field.setPlaceholderText('What changed, and why. This becomes the commit message and the pull request.')
        self.summary_field.setMaximumHeight(90)
        QVBoxLayout(summary).addWidget(self.summary_field)

        note = QLabel(
            'Proposing commits the release on its own branch, re-records the end-to-end snapshot, pushes the branch '
            'and opens the pull request page in your browser. The release is never tagged from here: merging the pull '
            'request is what publishes it.'
        )
        note.setWordWrap(True)

        self.propose_btn = QPushButton('Propose release')
        self.propose_btn.clicked.connect(self.start)
        buttons = QHBoxLayout()
        buttons.addStretch()
        buttons.addWidget(self.propose_btn)

        layout.addWidget(versions)
        layout.addWidget(changes, stretch=1)
        layout.addWidget(summary)
        layout.addWidget(note)
        layout.addLayout(buttons)
        return page

    def refresh(self) -> None:
        """Work out the release the database calls for now, or why there cannot be one."""
        try:
            proposal = prepare()
        except ReleaseError as error:
            self.show_problem(str(error))
        else:
            self.show_proposal(proposal)

    def show_problem(self, message: str) -> None:
        self.proposal = None
        self.problem_label.setText(message)
        self.pages.setCurrentWidget(self.problem_page)

    def show_proposal(self, proposal: Proposal) -> None:
        self.proposal = proposal

        self.level_label.setText(f'<b>{proposal.level}</b> - {LEVEL_EXPLANATIONS[proposal.level]}')
        self.db_version_label.setText(f'{proposal.db_version} &rarr; <b>{proposal.new_db_version}</b>')
        self.flumut_version_label.setText(f'{proposal.flumut_version} &rarr; <b>{proposal.new_flumut_version}</b>')
        self.branch_label.setText(f'<code>{proposal.branch}</code> into <code>main</code>')

        self.changes_tree.clear()
        for table, rows in proposal.changes.items():
            QTreeWidgetItem(self.changes_tree, [table, str(rows)])
        self.changes_tree.resizeColumnToContents(0)

        # Keep whatever the curator wrote; only the offered sentence follows the data.
        current = self.summary_field.toPlainText().strip()
        if not current or current == self.offered_summary:
            self.offered_summary = proposal.default_summary
            self.summary_field.setPlainText(self.offered_summary)

        self.pages.setCurrentWidget(self.proposal_page)

    def start(self) -> None:
        """Carry out the release. The window is unresponsive until the snapshot and the push are done."""
        if self.proposal is None:
            return
        summary = self.summary_field.toPlainText().strip()

        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        QApplication.processEvents()
        try:
            release(self.proposal, summary)
            url = pull_request_url(self.proposal, summary)
        except ReleaseError as error:
            QApplication.restoreOverrideCursor()
            ErrorDialog.show_error(
                self,
                'The release could not be completed',
                f'The checkout is on "{self.proposal.branch}" with whatever was written, so no work is lost.',
                str(error),
            )
        else:
            QApplication.restoreOverrideCursor()
            QDesktopServices.openUrl(QUrl(url))
            SuccessNotification.show_success(
                self,
                f'{self.proposal.title} was pushed on "{self.proposal.branch}", and its pull request page is open in '
                f'your browser.\n\nThe checkout is still on "{self.proposal.branch}": switch back to main when done.',
            )
            self.offered_summary = ''
            self.summary_field.clear()

        self.refresh()
