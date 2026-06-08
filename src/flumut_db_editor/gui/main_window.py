from PySide6.QtWidgets import QMainWindow, QTabWidget, QVBoxLayout, QWidget

from flumut_db_editor.gui.tabs.effects_tab import EffectsTab
from flumut_db_editor.gui.tabs.evidences_tab import EvidencesTab
from flumut_db_editor.gui.tabs.markers_tab import MarkersTab
from flumut_db_editor.gui.tabs.mutations_tab import MutationsTab
from flumut_db_editor.gui.tabs.papers_tab import PapersTab
from flumut_db_editor.gui.tabs.proteins_tab import ProteinsTab
from flumut_db_editor.gui.tabs.references_tab import ReferencesTab


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

        self.tabs.addTab(ProteinsTab(), 'Proteins')
        self.tabs.addTab(ReferencesTab(), 'References')
        self.tabs.addTab(PapersTab(), 'Papers')
        self.tabs.addTab(MutationsTab(), 'Mutations')
        self.tabs.addTab(MarkersTab(), 'Markers')
        self.tabs.addTab(EffectsTab(), 'Effects')
        self.tabs.addTab(EvidencesTab(), 'Evidences')
