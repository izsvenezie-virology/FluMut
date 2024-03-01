from collections import defaultdict
from click.types import File

from DataClass import Mutation


def matrix_output(output_file: File, mutations: list[Mutation]):
    lines = defaultdict(list)
    header = ['Sample']
    for mutation in mutations:
        if not mutation.found:
            continue
        header.append(mutation.name)
        for sample, sequence in mutation.samples.items():
            lines[sample].append(sequence)
    output_file.write('\t'.join(header))
    for sample, sequences in lines.items():
        output_file.write('\t'.join([sample] + sequences)+'\n')

def tabular_output(output_file: File, markers_per_sample) -> None:
    lines = []
    lines.append('Sample\teffect\tpaper\tsubtype\tfound_mutations\tmarker_mutations')
    for sample in markers_per_sample:
        for marker in markers_per_sample[sample]:
            lst = [sample] + list(marker.values())
            string = '\t'.join(lst)
            lines.append(string)
    out_str = '\n'.join(lines)

    output_file.encoding = 'utf-8'
    output_file.write(out_str)
