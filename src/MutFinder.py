#! /usr/bin/env python

import csv
import re
import sys
from collections import defaultdict
from io import TextIOWrapper
from typing import Dict, Generator, List, Tuple

import click
from Bio.Align import PairwiseAligner
from click import File

PRINT_ALIGNMENT = False


@click.command
@click.option('-n', '--name-regex', type=str, default=r'(?P<sample>.+)_(?P<segment>.+)')
@click.option('-M', '--markers-file', type=File('r'), default='src/data/markers.tsv')
@click.option('-R', '--references-fasta', type=File('r'), default='src/data/references.fa')
@click.option('-A', '--annotation-file', type=File('r'), default='src/data/annotations.tsv')
@click.option('-o', '--output', type=File('w'), default='-')
@click.argument('samples-fasta', type=File('r'))
def main(name_regex: str, markers_file: File, references_fasta: File, annotation_file: File, output: File, samples_fasta: File):
    seq_name_re = re.compile(name_regex)

    markers = load_markers(markers_file)
    mutations = load_mutations(markers)
    references = load_references(references_fasta)
    annotations = load_annotations(annotation_file)

    muts_per_sample = defaultdict(list)

    for name, seq in read_fasta(samples_fasta):
        sample,  segment = parse_name(name, seq_name_re)
        if segment not in references:
            print(f'Segment {segment} is not present in your reference. Sequence id: {name}', file=sys.stderr)
            continue
        ref_nucl, sample_nucl = pairwise_alignment(references[segment], seq)

        for protein, cds in annotations[segment].items():
            ref_coding, sample_coding = get_coding_sequences(
                ref_nucl, sample_nucl, cds)
            ref_aa, _ = translate(ref_coding)
            sample_aa, degenerations = translate(sample_coding)

            muts_per_sample[sample] += find_mutations(
                ref_aa, sample_aa, mutations[protein])

    markers_per_sample = defaultdict(list)
    for sample in muts_per_sample:
        markers_per_sample[sample] = get_markers_from_mutations(
            muts_per_sample[sample], markers)

    for sample in markers_per_sample:
        for marker in markers_per_sample[sample]:
            lst = [sample] + list(marker.values())
            string = '\t'.join(lst)
            print(string)


def load_markers(makers_file: click.File) -> List[Dict[str, str]]:
    return list(csv.DictReader(makers_file, delimiter='\t'))


def load_mutations(markers: dict) -> Dict[str, List[str]]:
    mutations = defaultdict(list)
    for marker in markers:
        muts = marker['Marker'].split(';')
        for mut in muts:
            segment = mut.split(':')[0]
            if mut not in mutations[segment]:
                mutations[segment].append(mut)
    return mutations


def load_references(ref_fasta: File) -> Dict[str, str]:
    refs = defaultdict(str)
    for name, seq in read_fasta(ref_fasta):
        refs[name] = seq
    return refs


def load_annotations(annotation_file: File) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    ann = defaultdict(lambda: defaultdict(list))
    for line in annotation_file:
        if line.startswith('#'):
            continue
        info = line.split('\t')
        ann[info[0]][info[1]].append((int(info[2]), int(info[3])))
    return ann


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
        pos = adjust_position(ref_aa, int(match.group('pos')) - 1)
        if sample_aa[pos] == match.group('alt'):
            found_mutations.append(mutation)
    return found_mutations


def adjust_position(ref_seq, pos):
    adjusted_pos = pos
    while True:
        dash_count = ref_seq.count('-', 0, adjusted_pos + 1)
        adjusted_pos = pos + dash_count
        if ref_seq.count('-', 0, adjusted_pos+1) == dash_count:
            break
    return adjusted_pos


def pairwise_alignment(ref_seq, sample_seq):
    aligner = PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    alignment = aligner.align(ref_seq, sample_seq)[0]
    if PRINT_ALIGNMENT:
        print(alignment, file=sys.stderr)
    return alignment[0], alignment[1]


def read_fasta(fasta_file: TextIOWrapper) -> Generator[str, str, None]:
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


def parse_name(name, pattern):
    match = pattern.match(name)
    return match.group('sample'), match.group('segment')



def translate(seq):
    degenerations = {}
    nucls = list(seq)
    aas = []

    for i in range(0, len(nucls), 3):
        triplet = nucls[i:i+3]
        if '-' in triplet:
            if triplet == ['-']*3:
                aas.append('-')
                continue
            if not aas and triplet[0] == '-':
                aas.append('?')
                continue
            triplet = [n for n in triplet if n != '-']
            while len(triplet) < 3:
                next_nucl = find_next_nucl(nucls, i)
                if not next_nucl:
                    break
                triplet.append(nucls[next_nucl])
                nucls[next_nucl] = '-'
        aa = translation_dict.get(''.join(triplet), '?')
        if aa == '?':
            degenerations[i/3] = ''.join(triplet)
        aas.append(aa)
    return ''.join(aas), degenerations


def find_next_nucl(seq, start):
    for i in range(start + 3, len(seq)):
        if not seq[i] == '-':
            return i
    return None


def get_coding_sequences(ref_seq, sample_seq, cds):
    ref_nucl = ''
    sample_nucl = ''

    for rng in cds:
        start = adjust_position(ref_seq, rng[0] - 1)
        end = adjust_position(ref_seq, rng[1] - 1)
        ref_nucl += ref_seq[start:end + 1]
        sample_nucl += sample_seq[start:end + 1]
    return ref_nucl, sample_nucl


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
