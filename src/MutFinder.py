#! /usr/bin/env python

import csv
import re
import subprocess
from io import TextIOWrapper
from pathlib import Path
from sys import stderr
from uuid import uuid4
import click
from collections import defaultdict


@click.command
@click.option('--tmp-path', type=str, default='.mutfinder_tmp')
@click.option('-r', '--name-regex', type=str, default=r'(?P<sample>.+)_(?P<segment>.+)')
@click.option('--markers-file', type=click.File('r'), default='src/data/markers.tsv')
@click.option('--references-fasta', type=click.File('r'), default='src/data/references.fa')
@click.argument('samples-fasta', type=click.File('r'))
def main(tmp_path, name_regex: str, markers_file: click.File, references_fasta: click.File, samples_fasta: click.File):
    tmp_dir = Path(tmp_path)
    tmp_dir.mkdir(exist_ok=True)
    pattern = re.compile(name_regex)

    markers = load_markers(markers_file)
    mutations = load_mutations(markers)  # TODO: extract mutations
    references = load_references(references_fasta)

    muts_per_sample = defaultdict(list)

    for name, seq in read_fasta(samples_fasta):
        sample,  segment = parse_name(name, pattern)

        for ref_protein, ref_sequence in references[segment].items():
            ref_nucl, sample_nucl = align_2_sequences(
                ref_protein, ref_sequence, name, seq, tmp_path)
            ref_nucl, sample_nucl = trim_alignment(ref_nucl, sample_nucl)

            ref_aa, sample_aa = map(translate, [ref_nucl, sample_nucl])
            ref_aa, sample_aa = align_2_sequences(
                ref_protein, ref_aa, name, sample_aa, tmp_path)

            # TODO: checks like len, frameshift, stalk deletion ecc.

            muts_per_sample[sample] += find_mutations(ref_aa, sample_aa, mutations[ref_protein])
    
    markers_per_sample = defaultdict(list)
    for sample in muts_per_sample:
        markers_per_sample[sample] = get_markers_from_mutations(muts_per_sample[sample], markers)
        
    for sample in markers_per_sample:
        for marker in markers_per_sample[sample]:
            lst = [sample] + list(marker.values())
            string = '\t'.join(lst)
            print(string)


    try:
        tmp_dir.rmdir()
    except:
        print(
            f'Tmp folder not removed because is not empty: {tmp_dir.absolute()}', file=stderr)
        pass


def get_markers_from_mutations(muts, markers):
    found_markers = []
    for marker in markers:
        marker['FoundMutations'] = []
        marker_muts = marker['Marker'].split(';')
        for mut in marker_muts:
            if mut in muts:
                marker['FoundMutations'].append(mut)
        if marker['FoundMutations']:
            marker['FoundMutations'] = ';'.join(marker['FoundMutations'])
            found_markers.append(marker)
    return found_markers


def find_mutations(ref_aa, sample_aa, mutations):
    found_mutations = []
    pattern = re.compile(r'^.+:.(?P<pos>\d+)(?P<alt>.)$')
    for mutation in mutations:
        match = pattern.match(mutation)
        pos = adjust_mutation_position(ref_aa, int(match.group('pos')) - 1)
        if sample_aa[pos] == match.group('alt'):
            found_mutations.append(mutation)
    return found_mutations


def adjust_mutation_position(ref_seq, pos):
    adjusted_pos = pos
    while True:
        dash_count = ref_seq.count('-', 0, adjusted_pos)
        adjusted_pos = pos + dash_count
        if ref_seq.count('-', 0, adjusted_pos) == dash_count:
            break
    return pos


def align_2_sequences(ref_name, ref_seq, sample_name, sample_seq, tmp_path):
    tmp_file = Path(tmp_path, f'{str(uuid4())}.fa')
    fasta = f'>{ref_name}\n{ref_seq}\n>{sample_name}\n{sample_seq}'
    try:
        with open(tmp_file, 'w') as fa:
            fa.write(fasta)
        ref_aligned, sample_aligned = mafft_alignment(tmp_file)
    finally:
        tmp_file.unlink()
    return ref_aligned, sample_aligned


def load_markers(makers_file: click.File):
    return list(csv.DictReader(makers_file, delimiter='\t'))


def load_mutations(markers: dict):
    mutations = defaultdict(list)
    for marker in markers:
        muts = marker['Marker'].split(';')
        for mut in muts:
            segment = mut.split(':')[0]
            if mut not in mutations[segment]:
                mutations[segment].append(mut)
    return mutations


def read_fasta(fasta_file: TextIOWrapper):
    '''Create a Fasta reading a file in Fasta format'''
    name = None
    for line in fasta_file:
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(seq)
            name = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def load_references(ref_fasta: click.File):
    ref = defaultdict(dict)
    for name, seq in read_fasta(ref_fasta):
        name_info = name.split('-')
        ref[name_info[0]][name] = seq
    return ref


def parse_name(name, pattern):
    match = pattern.match(name)
    return match.group('sample'), match.group('segment')


def mafft_alignment(fasta_path: str):
    mafft_proc = subprocess.run(
        f'mafft --auto --thread 1 {fasta_path}'.split(), capture_output=True)
    pattern = re.compile(r'>.+?\n(?P<ref>.+)\n>.+?\n(?P<sample>.+)', re.DOTALL)
    stdout = mafft_proc.stdout.decode('utf-8').upper()
    match = pattern.match(stdout)
    return match.group('ref').replace('\n', ''), match.group('sample').replace('\n', '')


def trim_alignment(ref, sample):
    ref_trim = ref.lstrip('-')
    l_offset = len(ref) - len(ref_trim)
    sample_trim = sample[l_offset:]
    return ref_trim, sample_trim


def translate(seq_nucl):
    seq_aa = []
    for i in range(0, len(seq_nucl), 3):
        seq_aa.append(translation_dict.get(seq_nucl[i:i+3], '?'))
    return ''.join(seq_aa)


translation_dict = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


if __name__ == '__main__':
    main()
