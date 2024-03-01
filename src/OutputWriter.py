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
