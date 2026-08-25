from flumut.flumutdb.models import Segment
from flumut_db_editor.gui.forms.base import EvidenceTermsForm


class SegmentForm(EvidenceTermsForm[Segment]):
    model = Segment
