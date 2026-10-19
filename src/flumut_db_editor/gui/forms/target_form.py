from flumut.flumutdb.models import Target
from flumut_db_editor.gui.forms.base import EvidenceTermsForm


class TargetForm(EvidenceTermsForm[Target]):
    model = Target
