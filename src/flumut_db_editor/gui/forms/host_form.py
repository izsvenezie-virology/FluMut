from flumut.flumutdb.models import Host
from flumut_db_editor.gui.forms.base import EvidenceTermsForm


class HostForm(EvidenceTermsForm[Host]):
    model = Host
