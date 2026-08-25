from flumut.flumutdb.models import Effect
from flumut_db_editor.gui.forms.effect_form import EffectForm
from flumut_db_editor.gui.tabs.base import EvidenceTermsTab


class EffectsTab(EvidenceTermsTab[Effect]):
    model = Effect
    form = EffectForm
