from flumut.flumutdb.models import Subtype
from flumut_db_editor.gui.forms.base import NameForm


class SubtypeForm(NameForm[Subtype]):
    model = Subtype
