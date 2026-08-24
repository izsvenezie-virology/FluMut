from flumut.flumutdb.models import Segment
from flumut_db_editor.gui.forms.base import NameForm


class SegmentForm(NameForm[Segment]):
    model = Segment
