#! /usr/bin/env python

import csv
import itertools
import re
import sys
from collections import defaultdict
from io import TextIOWrapper
from typing import Dict, Generator, List, Tuple

import click
from Bio.Align import PairwiseAligner
from click import File

PRINT_ALIGNMENT = False


@click.command()
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
            print(
                f'Segment {segment} is not present in your reference. Sequence id: {name}', file=sys.stderr)
            continue
        ref_nucl, sample_nucl = pairwise_alignment(references[segment], seq)

        for protein, cds in annotations[segment].items():
            ref_coding, sample_coding = get_coding_sequences(
                ref_nucl, sample_nucl, cds)
            ref_aa: str = ''.join(translate(ref_coding))
            sample_aa = translate(sample_coding)

            muts_per_sample[sample] += find_mutations(
                ref_aa, sample_aa, mutations[protein])

    markers_per_sample = defaultdict(list)
    for sample in muts_per_sample:
        markers_per_sample[sample] = match_markers(
            muts_per_sample[sample], markers)

    lines = []
    lines.append(
        '\t'.join(['Sample'] + list(markers[0].keys()) + ['Mutations present']))
    for sample in markers_per_sample:
        for marker in markers_per_sample[sample]:
            lst = [sample] + list(marker.values())
            string = '\t'.join(lst)
            lines.append(string)
    print('\n'.join(lines), file=output)


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
    return {name: seq for name, seq in read_fasta(ref_fasta)}


def load_annotations(annotation_file: File) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    ann = defaultdict(lambda: defaultdict(list))
    for line in annotation_file:
        if line.startswith('#'):
            continue
        info = line.split('\t')
        ann[info[0]][info[1]].append((int(info[2]), int(info[3])))
    return ann


def match_markers(muts, markers):
    found_markers = []
    for marker in markers:
        found_muts = []
        for mut in marker['Marker'].split(';'):
            if mut in muts:
                found_muts.append(mut)
        if found_muts:
            mrk = marker.copy()
            mrk['FoundMutations'] = ';'.join(found_muts)
            found_markers.append(mrk)
    return found_markers


def find_mutations(ref_aa, sample_aa, mutations):
    found_mutations = []
    pattern = re.compile(r'^.+:.(?P<pos>\d+)(?P<alt>.)$')
    for mutation in mutations:
        match = pattern.match(mutation)
        pos = adjust_position(ref_aa, int(match.group('pos')))
        if match.group('alt') in sample_aa[pos]:
            found_mutations.append(mutation)
    return found_mutations


def pairwise_alignment(ref_seq: str, sample_seq: str) -> Tuple[str, str]:
    '''Align sequence against a reference'''
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
                yield name, ''.join(seq).upper()
            name = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq).upper()


def parse_name(name: str, pattern: re.Pattern) -> Tuple[str, str]:
    '''Get sample and segment information by sequence name'''
    match = pattern.match(name)
    return match.group('sample'), match.group('segment')


def translate(seq: str) -> List[str]:
    '''Translate nucleotidic sequence in AA sequence'''
    nucls = list(seq)
    aas = []
    is_first = True

    for i in range(0, len(nucls), 3):
        codon = get_codon(nucls, i, is_first)
        aa = translate_codon(codon)
        if aa not in ('-', '?'):
            is_first = False
        aas.append(aa)
    return aas


def get_codon(seq: List[str], start: int, is_first: bool) -> [str, str, str]:
    '''Exctract the codon'''
    codon = seq[start:start + 3]
    if codon == ['-', '-', '-']:  # If the codon is a deletion
        return codon
    # If the codon starts from mid codon (to avoid frameshifts in truncated sequences):
    if is_first and codon[0] == '-':
        return codon
    codon = [n for n in codon if n != '-']
    while len(codon) < 3:
        if not (next_nucl := find_next_nucl(seq, start)):
            break
        codon.append(seq[next_nucl])
        seq[next_nucl] = '-'
    return codon


def translate_codon(codon: List[str]) -> str:
    '''Translate a codon into a set of AAs, containing all possible combinations in case of degenerations'''
    if 'N' in codon:
        return '?'
    undegenerated_codon = [degeneration_dict[nucl] for nucl in codon]
    codons = list(itertools.product(*undegenerated_codon))
    aas = [translation_dict.get(''.join(c), '?') for c in codons]
    return ''.join(set(aas))


def find_next_nucl(seq: List[str], start: int):
    '''Returns the position of the next non deleted nucleotide'''
    for i in range(start + 3, len(seq)):
        if not seq[i] == '-':
            return i
    return None


def get_coding_sequences(ref_seq: str, sample_seq: str, cds: List[Tuple[int, int]]):
    '''Cut and assemble the nucleotide sequences based on positions given by the cds'''
    cds.sort(key=lambda x: x[0])
    ref_nucl = ''
    sample_nucl = ''

    for rng in cds:
        start = adjust_position(ref_seq, rng[0])
        end = adjust_position(ref_seq, rng[1]) + 1
        ref_nucl += ref_seq[start:end]
        sample_nucl += sample_seq[start:end]
    return ref_nucl, sample_nucl


def adjust_position(ref_seq: str, pos: int) -> int:
    '''Adjust 1-based position to 0-based, considering reference sequence gaps'''
    pos -= 1    # conversion 1-based to 0-based numeration
    dashes = 0
    adj_pos = pos
    while ref_seq.count('-', 0, adj_pos + 1) != dashes:
        adj_pos = pos + (dashes := ref_seq.count('-', 0, adj_pos + 1))
    return adj_pos


translation_dict = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    '---': '-'
}


degeneration_dict = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'U': ['U'], '-': ['-'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}


if __name__ == '__main__':
    main()
