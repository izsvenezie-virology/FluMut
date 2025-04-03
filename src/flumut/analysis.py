
from io import TextIOWrapper
from typing import Dict, List, Tuple

from flumut.db_utility.db_data import markers_by_mutations
from flumut.output import write_outputs
from flumut.sequence_utility.aligner import align
from flumut.sequence_utility.fasta_handler import parse_header, read_fasta
from flumut.sequence_utility.models import FastaSequence, Sample
from flumut.sequence_utility.parser import seek_mutations
from flumut.sequence_utility.translator import translate


def analyze(fastas: Tuple[TextIOWrapper], relaxed: bool) -> None:
    sequences = collect_sequences(fastas)
    samples = parse_sequences(sequences)
    seek_markers(samples, relaxed)
    write_outputs(samples)


def collect_sequences(fastas: Tuple[TextIOWrapper]) -> List[FastaSequence]:
    sequences: List[FastaSequence] = []
    for fasta in fastas:
        sequences += read_fasta(fasta)
    return sequences


def parse_sequences(sequences: List[FastaSequence]) -> List[Sample]:
    samples: Dict[str, Sample] = {}

    for sequence in sequences:
        sample, segment = parse_header(sequence.header)
        alignment = align(sequence.sequence, segment)
        sequence.alignment = alignment
        proteins = translate(alignment)
        alignment.proteins = proteins

        if sample not in samples:
            samples[sample] = Sample(sample)
        samples[sample].sequences.append(sequence)

    return list(samples.values())


def seek_markers(samples: List[Sample], relaxed: bool) -> None:
    for sample in samples:
        sample.mutations = seek_mutations(sample)
        sample.markers = markers_by_mutations(sample.mutations, relaxed)
