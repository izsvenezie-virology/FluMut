from flumut.flumutdb.models import Effect
from flumut_db_editor.gui.forms.base import EvidenceTermsForm


class EffectForm(EvidenceTermsForm[Effect]):
    model = Effect
