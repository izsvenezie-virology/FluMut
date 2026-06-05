from dataclasses import dataclass, field

from flumut.core.globals import GAP_SYMBOL
from flumut.core.io.input import FastaSequence
from flumut.flumutdb import Protein, Reference


@dataclass
class Alignment:
    """Pairwise alignment between a reference and a query sequence as character lists.

    Both lists are gap-padded to equal length; gap positions are represented by
    the GAP_SYMBOL constant.

    Attributes:
        reference: Aligned reference sequence as a list of single characters.
        query: Aligned query sequence as a list of single characters.
    """

    reference: list[str] = field(default_factory=list)
    query: list[str] = field(default_factory=list)

    def get_positions(self) -> list[int]:
        """Map each alignment column to its 1-based reference position.

        Gap columns in the reference (insertions in the query) are mapped to
        None. Non-gap columns are numbered sequentially starting at 1.

        Returns:
            A list of the same length as ``reference``, where each element is
            the 1-based position number or None for gap columns.
        """
        positions = []
        last_position = 0
        for r in self.reference:
            if r == GAP_SYMBOL:
                positions.append(None)
            else:
                last_position += 1
                positions.append(last_position)
        return positions


@dataclass
class Nucleotide:
    """Bundles a query FASTA sequence with its best-matching reference and alignment.

    Attributes:
        query: The original query FASTA sequence.
        reference: The reference sequence chosen as the best alignment target.
        alignment: The pairwise nucleotide alignment between reference and query.
    """

    query: FastaSequence
    """Query sequence."""
    reference: Reference
    """Reference sequence."""
    alignment: Alignment


@dataclass
class CDS:
    """Coding sequence alignment slice for a specific protein.

    Contains the nucleotide alignment trimmed to the annotated CDS boundaries
    for one protein, along with any frameshifts detected during frame adjustment.

    Attributes:
        nucleotides: The full-length nucleotide alignment this CDS was extracted from.
        protein: The protein whose CDS boundaries define this slice.
        alignment: The nucleotide alignment restricted to the CDS region.
        frameshifts: Ranges of 1-based nucleotide positions affected by frameshifts,
            as ``(start, end)`` tuples where ``end`` is None for open-ended frameshifts.
    """

    nucleotides: Nucleotide
    protein: Protein
    alignment: Alignment

    # From nucleotide alignment properties
    frameshifts: list[tuple[int, int | None]] = field(default_factory=list)
