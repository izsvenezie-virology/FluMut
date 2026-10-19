from flumut.flumutdb.models import Subtype
from flumut_db_editor.gui.forms.subtype_form import SubtypeForm
from flumut_db_editor.gui.tabs.base import EvidenceTermsTab


class SubtypesTab(EvidenceTermsTab[Subtype]):
    model = Subtype
    form = SubtypeForm
