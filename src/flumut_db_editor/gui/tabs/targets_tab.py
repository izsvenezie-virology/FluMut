from flumut.flumutdb.models import Target
from flumut_db_editor.gui.forms.target_form import TargetForm
from flumut_db_editor.gui.tabs.base import EvidenceTermsTab


class TargetsTab(EvidenceTermsTab[Target]):
    model = Target
    form = TargetForm
