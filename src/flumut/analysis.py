import logging
from io import TextIOWrapper
from typing import Dict, List, Tuple

import flumutdb

from flumut.db_utility.db_data import markers_by_mutations
from flumut.exceptions import UnmatchingHeaderException
from flumut.output import write_outputs
from flumut.sequence_utility.aligner import align, select_candidate_references
from flumut.sequence_utility.fasta_handler import get_header_pattern, parse_header, read_fasta
from flumut.sequence_utility.models import FastaSequence, Sample
from flumut.sequence_utility.parser import seek_mutations
from flumut.sequence_utility.translator import translate


def analyze(fastas: Tuple[TextIOWrapper], relaxed: bool, allow_unmatching_headers: bool) -> None:
    flumutdb.initialize()

    sequences = collect_sequences(fastas)
    samples = parse_sequences(sequences, allow_unmatching_headers)
    seek_markers(samples, relaxed)
    write_outputs(samples)


def collect_sequences(fastas: Tuple[TextIOWrapper]) -> List[FastaSequence]:
    logging.info(f'Collecting sequences from {len(fastas)} Fasta files')
    sequences: List[FastaSequence] = []
    for fasta in fastas:
        sequences += read_fasta(fasta)
    logging.info(f'Collected {len(sequences)} sequences')
    return sequences


def parse_sequences(sequences: List[FastaSequence], allow_unmatching_headers: bool) -> List[Sample]:
    logging.info('Parsing sequences...')
    samples: Dict[str, Sample] = {}

    for sequence in sequences:
        logging.info(f'Parsing "{sequence.header}" from {sequence.file}.')
        sample, segment = parse_header(sequence.header)

        if sample is None:
            if not allow_unmatching_headers:
                raise UnmatchingHeaderException(sequence.header, get_header_pattern()) from None
            logging.info(f'Cannot extract sample ID from "{sequence.header}". The whole header is used as sample ID')
            sample = sequence.header

        candidates = select_candidate_references(segment)
        alignment = align(sequence.sequence, candidates)
        sequence.alignment = alignment
        logging.info(f'Aligned "{sequence.header}" on "{alignment.reference.name}" with score {alignment.score:.0f}')
        if alignment.score < len(sequence.sequence) * 0.5:
            logging.warning(
                f'Low alignment quality between "{sequence.header}" and reference "{alignment.reference.name}" (score {alignment.score:.0f})'
            )

        proteins = translate(alignment)
        alignment.proteins = proteins
        for protein in proteins:
            for frameshift in protein.frameshifts:
                logging.warning(
                    f'Frameshift in "{sequence.header}" from position {frameshift[0]} to {frameshift[1]} in protein "{protein.name}" (length {len(protein.sequence)})'
                )

        if sample not in samples:
            samples[sample] = Sample(sample)
        samples[sample].sequences.append(sequence)

    return list(samples.values())


def seek_markers(samples: List[Sample], relaxed: bool) -> None:
    logging.info('Seeking markers...')
    for sample in samples:
        sample.mutations = seek_mutations(sample)
        sample.markers = markers_by_mutations(sample.mutations, relaxed)
