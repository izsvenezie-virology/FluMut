import csv
from typing import Dict, List, Tuple
from click.types import File

from DataClass import Mutation, Sample


def matrix_format(mutations: List[Mutation]) -> Tuple[List[str], List[Dict[str, str]]]:
    header = ['Sample']
    samples = {}
    for mutation in mutations:
        if mutation.found and mutation.name not in header:
            header.append(mutation.name)
        for sample in mutation.samples:
            if sample not in samples:
                samples[sample] = {'Sample': sample}
            samples[sample][mutation.name] = mutation.samples[sample]
    return header, samples.values()


def tabular_output(samples: List[Sample]) -> Tuple[List[str], List[Dict[str, str]]]:
    header = ['Sample', 'Marker mutations', 'Found mutations', 'Effect', 'Subtype', 'Papers']
    data = []
    for sample in samples:
        for marker in sample.markers:
            marker['Sample'] = sample.name
            data.append(marker)
    return header, data


def write_csv(output_file: File, header: List[str], data: List[Dict[str, str]]) -> None:
    writer = csv.DictWriter(output_file, header, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(data)    
