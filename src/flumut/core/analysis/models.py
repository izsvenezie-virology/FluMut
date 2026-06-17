from dataclasses import dataclass, field

from flumut.core.logger import LOGGER
from flumut.core.nucleotides.models import CDS, Alignment, Protein, Reference
from flumut.flumutdb import Mapping, Marker, Mutation, MutationType, Paper


@dataclass
class PositionScan:
    """Scan result for a single mutation position in a protein alignment.

    Attributes:
        mapping: Database mapping that links a reference position to the expected alteration.
        ammino_acid: Observed amino acid(s) in the query at this position.
        is_detected: True when the expected alteration is present in ``ammino_acid``.
    """

    mapping: Mapping
    ammino_acid: str
    is_detected: bool = field(init=False)

    def __post_init__(self) -> None:
        """Compute ``is_detected`` based on the mutation type of the mapping."""
        match self.mapping.mutation.type:
            case MutationType.SNP.value:
                self.is_detected = self.mapping.alteration in self.ammino_acid
            case _:
                raise NotImplementedError(f'Mutation type {self.mapping.mutation.type} not supported')

    @property
    def mutation(self) -> Mutation:
        """The Mutation associated with this position scan."""
        return self.mapping.mutation


@dataclass
class MarkerScan:
    """Scan result for a complete marker across all its positions.

    Attributes:
        marker: The Marker definition being evaluated.
        positions: Scanned position results for each mutation in the marker.
        detected_mutations: Subset of positions where the mutation is detected.
        is_detected: True when at least one mutation in the marker is detected.
        is_complete: True when all mutations in the marker are detected.
    """

    marker: Marker
    positions: list[PositionScan]

    detected_mutations: list[PositionScan] = field(default_factory=list, init=False)
    is_detected: bool = field(default=False, init=False)
    is_complete: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        """Derive ``detected_mutations``, ``is_detected``, and ``is_complete`` from positions."""
        self.detected_mutations = [position for position in self.positions if position.is_detected]
        self.is_detected = len(self.detected_mutations) > 0
        self.is_complete = len(self.detected_mutations) == len(self.marker.mutations)


@dataclass
class ProteinAlignment:
    """Holds the amino acid alignment for a single protein within a sample.

    Attributes:
        protein: The Protein whose CDS was translated.
        reference: The nucleotide Reference used for alignment.
        alignment: The resulting amino acid-level Alignment.
        cds: The underlying CDS used to produce the alignment, if available.
    """

    protein: Protein
    reference: Reference
    alignment: Alignment

    cds: CDS | None = None


@dataclass
class Sample:
    """Aggregates all per-protein alignments and scan results for one sequenced sample.

    Attributes:
        id: Sample identifier, derived from the FASTA header.
        alignments: Per-protein amino acid alignments computed from the nucleotide sequence.
        positions: Flat list of all PositionScan results across all alignments.
        marker_scans: Marker detection results for this sample.
        checks: Quality-check issues detected for this sample.
    """

    id: str
    alignments: list[ProteinAlignment] = field(default_factory=list)

    positions: list[PositionScan] = field(default_factory=list, init=False)
    marker_scans: list[MarkerScan] = field(default_factory=list, init=False)
    checks: list['Check'] = field(default_factory=list, init=False)


@dataclass
class Analysis:
    """Top-level container for a complete FluMut analysis run.

    Attributes:
        samples: Map from sample ID to its Sample data.
        mutations: Set of all detected Mutation objects across all samples.
        markers: Set of all detected Marker objects across all samples.
        literature: Set of Paper references cited by detected markers.
    """

    samples: dict[str, Sample] = field(default_factory=dict)
    mutations: set[Mutation] = field(default_factory=set)
    markers: set[Marker] = field(default_factory=set)
    literature: set[Paper] = field(default_factory=set)


class Check:
    """Base class for a quality-check warning raised during sample validation.

    Subclasses construct a human-readable message and pass it to this
    constructor, which stores the message and emits it as a WARNING log entry.
    """

    def __init__(self, message: str):
        """Store the check message and emit it as a warning.

        Args:
            message: Human-readable description of the issue detected.
        """
        self.message = message
        LOGGER.warning(self.message)


class TruncationCheck(Check):
    """Quality check raised when a premature stop codon is detected."""

    def __init__(self, sample_id: str, protein_name: str, position: int):
        """Create a truncation check for the given sample, protein, and position.

        Args:
            sample_id: ID of the sample containing the premature stop codon.
            protein_name: Name of the protein in which truncation was detected.
            position: 1-based amino acid position of the premature stop codon.
        """
        message = f'Premature stop codon detected in sample "{sample_id}" for protein "{protein_name}" at {position}'
        super().__init__(message)


class EnlongationCheck(Check):
    """Quality check raised when a protein extends beyond the expected stop codon."""

    def __init__(self, sample_id: str, protein_name: str):
        """Create an elongation check for the given sample and protein.

        Args:
            sample_id: ID of the sample containing the elongation.
            protein_name: Name of the protein in which elongation was detected.
        """
        message = f'Enlongation detected in sample "{sample_id}" for protein "{protein_name}"'
        super().__init__(message)


class FrameshiftCheck(Check):
    """Quality check raised when a frameshift indel is detected in a CDS."""

    def __init__(self, sample_id: str, protein_name: str, start: int, end: int | None):
        """Create a frameshift check for the given sample, protein, and range.

        Args:
            sample_id: ID of the sample containing the frameshift.
            protein_name: Name of the protein in which the frameshift was detected.
            start: 1-based nucleotide start position of the frameshift region.
            end: 1-based nucleotide end position, or None if the frameshift
                extends to the end of the sequence.
        """
        message = f'Frameshift detected in sample "{sample_id}" for protein "{protein_name}" from position {start} to {end or "end"}'
        super().__init__(message)


class DuplicationCheck(Check):
    """Quality check raised when the same protein appears more than once in a sample."""

    def __init__(self, sample_id: str, protein_name: str):
        """Create a duplication check for the given sample and protein.

        Args:
            sample_id: ID of the sample containing the duplicate protein.
            protein_name: Name of the protein that appears multiple times.
        """
        message = f'Duplicate protein "{protein_name}" detected in sample "{sample_id}"'
        super().__init__(message)


class CDSHasNoValidCodonsCheck(Check):
    """Quality check raised when the CDS has no valid codons."""

    def __init__(self, sample_id: str, sequence_header: str, protein_name: str):
        """Create a check for CDS without valid codons for the given sample and protein.

        Args:
            sample_id: ID of the sample containing the CDS without valid codons.
            sequence_header: Header of the sequence containing the CDS without valid codons.
            protein_name: Name of the protein with the CDS without valid codons.
        """
        message = f'"{sequence_header}" has no valid codon in CDS for protein {protein_name} in sample {sample_id}'
        super().__init__(message)
