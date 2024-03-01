import csv
from click.types import File

from DataClass import Mutation


def matrix_output(output_file: File, mutations: list[Mutation]):
    samples = {}
    header = ['Sample']
    for mutation in mutations:
        if mutation.found and mutation.name not in header:
            header.append(mutation.name)
        for sample in mutation.samples:
            if sample not in samples:
                samples[sample] = {'Sample': sample}
            samples[sample][mutation.name] = mutation.samples[sample]
    writer = csv.DictWriter(output_file, header, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(samples.values())


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
