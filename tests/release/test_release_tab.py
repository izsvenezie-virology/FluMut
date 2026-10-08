"""Smoke tests for the editor's release tab.

The tab carries none of the logic, so these only check that it builds and shows
the proposal it was given. They are skipped wherever PySide6 is not installed,
which includes the packaging CI.
"""

import os
from pathlib import Path

import pytest

pytest.importorskip('PySide6', reason='the database editor GUI is an optional extra')

from PySide6.QtWidgets import QApplication

from flumut_db_editor.release import Proposal

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture(scope='module')
def application() -> QApplication:
    """One Qt application for the module: a second instance would abort."""
    return QApplication.instance() or QApplication([])


@pytest.fixture
def proposal() -> Proposal:
    return Proposal(Path('.'), 'patch', {'evidence': 4, 'paper': 1}, '7.0.0', '1.0.0')


@pytest.fixture
def tab(application: QApplication, proposal: Proposal):
    from flumut_db_editor.gui.tabs.release_tab import ReleaseTab

    built = ReleaseTab()
    built.show_proposal(proposal)
    yield built
    built.deleteLater()


def test_the_tab_names_the_new_versions(tab) -> None:
    assert '7.0.1' in tab.db_version_label.text()
    assert '1.0.1' in tab.flumut_version_label.text()


def test_the_summary_starts_from_the_computed_change(tab) -> None:
    """The curator edits a sentence rather than writing one from nothing."""
    assert tab.summary_field.toPlainText() == 'Patch database update: evidence (4), paper (1).'


def test_a_refresh_keeps_the_summary_the_curator_wrote(tab, proposal: Proposal) -> None:
    tab.summary_field.setPlainText('Add the 2024 H5N1 evidence')
    tab.show_proposal(proposal)

    assert tab.summary_field.toPlainText() == 'Add the 2024 H5N1 evidence'


def test_every_changed_table_is_listed(tab) -> None:
    rows = [tab.changes_tree.topLevelItem(row).text(0) for row in range(tab.changes_tree.topLevelItemCount())]

    assert rows == ['evidence', 'paper']


def test_a_problem_replaces_the_proposal(tab) -> None:
    tab.show_problem('The checkout is on "dbeditor", not "main".')

    assert tab.pages.currentWidget() is tab.problem_page
    assert tab.proposal is None
