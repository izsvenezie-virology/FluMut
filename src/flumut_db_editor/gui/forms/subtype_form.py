from flumut.flumutdb.models import Subtype
from flumut_db_editor.gui.forms.base import EvidenceTermsForm


class SubtypeForm(EvidenceTermsForm[Subtype]):
    model = Subtype
