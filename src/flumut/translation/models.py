from dataclasses import dataclass

from flumutdb.models import Protein

from flumut.alignment.models import NucleotideAlignment
from flumut.core.models import Alignment


@dataclass
class AlignedProtein(Alignment):
    alignment: NucleotideAlignment
    protein: Protein
