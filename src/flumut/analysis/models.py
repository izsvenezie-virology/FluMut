import re
from dataclasses import dataclass, field
from io import TextIOWrapper

from flumut.alignment import NucleotideAlignment
from flumut.alignment.aligner import get_best_alignment, select_candidate_references
from flumut.analysis.parser import parse_header
from flumut.core.models import Alignment
from flumut.flumutdb import Marker, Mutation, Paper, Protein, Reference
from flumut.io.input import read_fasta
from flumut.scan import MarkerScan, PositionScan
from flumut.scan.scanner import scan_markers, scan_positions
from flumut.translation.translator import translate


@dataclass
class ProteinAlignment(Alignment):
    protein: Protein
    reference: Reference

    # From nucleotide alignment properties
    nucleotides: NucleotideAlignment | None = None
    frameshifts: list[tuple[int, int]] = field(default_factory=list)


@dataclass
class Sample:
    id: str
    alignments: list[ProteinAlignment] = field(default_factory=list)

    positions: list[PositionScan] = field(default_factory=list, init=False)
    marker_scans: list[MarkerScan] = field(default_factory=list, init=False)


@dataclass
class Analysis:
    samples: dict[str, Sample] = field(default_factory=dict)
    mutations: set[Mutation] = field(default_factory=set)
    markers: set[Marker] = field(default_factory=set)
    literature: set[Paper] = field(default_factory=set)

    def load_nucleotide_fasta(self, fasta: TextIOWrapper, header_pattern: re.Pattern) -> None:
        for sequence in read_fasta(fasta):
            sample, candidate_hint = parse_header(sequence.name, header_pattern)
            candidates = select_candidate_references(candidate_hint)
            nt_alignment = get_best_alignment(sequence, candidates)

            if sample not in self.samples:
                self.samples[sample] = Sample(sample)

            for protein in nt_alignment.reference.segment.proteins:
                aa_alignment = translate(nt_alignment, protein)
                self.samples[sample].alignments.append(aa_alignment)

    def analyse(self, relaxed: bool = True) -> None:
        self.mutations.clear()
        self.markers.clear()
        self.literature.clear()

        for sample in self.samples.values():
            for alignment in sample.alignments:
                sample.positions += scan_positions(alignment)
            sample.marker_scans = scan_markers(sample.positions, relaxed)

            for position in sample.positions:
                if position.is_detected:
                    self.mutations.add(position.mutation)

            for scan in sample.marker_scans:
                self.markers.add(scan.marker)

        for marker in self.markers:
            for evidence in marker.evidences:
                self.literature.add(evidence.paper)
