from PySide6.QtWidgets import QMainWindow, QTabWidget, QVBoxLayout, QWidget

from flumut_db_editor.gui.tabs.base import BaseTreeTab
from flumut_db_editor.gui.tabs.effects_tab import EffectsTab
from flumut_db_editor.gui.tabs.hosts_tab import HostsTab
from flumut_db_editor.gui.tabs.mutations_tab import MutationsTab
from flumut_db_editor.gui.tabs.papers_tab import PapersTab
from flumut_db_editor.gui.tabs.proteins_tab import ProteinsTab
from flumut_db_editor.gui.tabs.references_tab import ReferencesTab
from flumut_db_editor.gui.tabs.subtypes_tab import SubtypesTab


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('FluMut Database Editor')
        self.setGeometry(100, 100, 1200, 800)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout(central_widget)
        layout.setContentsMargins(0, 0, 0, 0)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        evidence_terms_tab = QTabWidget()
        evidence_terms_tab.addTab(EffectsTab(), 'Effects')
        evidence_terms_tab.addTab(SubtypesTab(), 'Subtypes')
        evidence_terms_tab.addTab(HostsTab(), 'Hosts')

        self.tabs.addTab(ProteinsTab(), 'Proteins')
        self.tabs.addTab(ReferencesTab(), 'References')
        self.tabs.addTab(MutationsTab(), 'Mutations')
        self.tabs.addTab(evidence_terms_tab, 'Evidence Terms')
        # self.tabs.addTab(MarkersTab(), 'Markers')
        self.tabs.addTab(PapersTab(), 'Papers')
        # self.tabs.addTab(EvidencesTab(), 'Evidences')

        self.tabs.currentChanged.connect(self.on_tab_changed)

    def on_tab_changed(self):
        match widget := self.tabs.currentWidget():
            case BaseTreeTab():
                widget.refresh()
            case _:
                return
