from flumut.flumutdb.models import Host
from flumut_db_editor.gui.forms.host_form import HostForm
from flumut_db_editor.gui.tabs.base import EvidenceTermsTab


class HostsTab(EvidenceTermsTab[Host]):
    model = Host
    form = HostForm
