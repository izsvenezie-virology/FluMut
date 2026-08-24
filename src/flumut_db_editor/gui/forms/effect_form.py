from flumut.flumutdb.models import Effect
from flumut_db_editor.gui.forms.base import NameForm


class EffectForm(NameForm[Effect]):
    model = Effect
