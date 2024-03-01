import csv
from typing import Dict, List, Tuple
from click.types import File

from DataClass import Mutation


def matrix_format(mutations: list[Mutation]) -> Tuple[List[str], List[Dict[str, str]]]:
    samples = {}
    header = ['Sample']
    for mutation in mutations:
        if mutation.found and mutation.name not in header:
            header.append(mutation.name)
        for sample in mutation.samples:
            if sample not in samples:
                samples[sample] = {'Sample': sample}
            samples[sample][mutation.name] = mutation.samples[sample]
    return header, samples.values()


def tabular_output(output_file: File, markers_per_sample) -> None:
    lines = []
    lines.append('Sample\tMarker mutations\tFound mutations\tEffect\tPapers\tSubtype')
    for sample in markers_per_sample:
        for marker in markers_per_sample[sample]:
            lst = [sample] + list(marker.values())
            string = '\t'.join(lst)
            lines.append(string)
    out_str = '\n'.join(lines)
    output_file.write(out_str)

def write_csv(output_file: File, header: List[str], data: List[Dict[str, str]]) -> None:
    writer = csv.DictWriter(output_file, header, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(data)    
