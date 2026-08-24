from flumut.flumutdb.models import Host
from flumut_db_editor.gui.forms.base import NameForm


class HostForm(NameForm[Host]):
    model = Host
