#! /usr/bin/env python

import csv
import itertools
import re
import sys
from collections import defaultdict
from io import TextIOWrapper
from typing import Dict, Generator, List, Tuple
from importlib.resources import files
import sqlite3
import click
from Bio.Align import PairwiseAligner
from click import File

PRINT_ALIGNMENT = False
SKIP_UNMATCH_NAMES_OPT = '--skip-unmatch-names'
SKIP_UNKNOWN_SEGMENTS_OPT = '--skip-unknown-segments'


class Mutation:
    def __init__(self, name: str, type: str, ref: str, alt: str, pos: int) -> None:
        self.name: str = name
        self.type: str = type
        self.ref: str = ref
        self.alt: str = alt
        self.pos: int = pos


@click.command()
@click.option(SKIP_UNMATCH_NAMES_OPT, is_flag=True, default=False, help='Skips sequences with name that does not match the pattern')
@click.option(SKIP_UNKNOWN_SEGMENTS_OPT, is_flag=True, default=False, help='Skips sequences with name that does not match the pattern')
@click.option('-n', '--name-regex', type=str, default=r'(?P<sample>.+)_(?P<segment>.+)', show_default=True, help='Regular expression to parse sequence name')
@click.option('-D', '--db-file', type=str, default=files('data').joinpath('mutfinderDB.sqlite'), help='Source database')
@click.option('-o', '--output', type=File('w'), default='-', help='The output file [default: stdout]')
@click.argument('samples-fasta', type=File('r'))
def main(name_regex: str, output: File, samples_fasta: File, db_file: str,
         skip_unmatch_names: bool, skip_unknown_segments: bool):
    '''
    Search for markers of interest in the SAMPLES-FASTA file.
    '''

    # Initialization
    pattern = re.compile(name_regex)
    
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    mutations = load_mutations(cur)
    references = load_references(cur)
    annotations = load_annotations(cur)
    conn.close()

    muts_per_sample = defaultdict(list)
    markers_per_sample = defaultdict(list)

    # Per sequence analysis
    for name, seq in read_fasta(samples_fasta):
        sample,  segment = parse_name(name, pattern, skip_unmatch_names)
        if sample is None or segment is None:
            continue
        if segment not in references:
            print(f'Unknown segment {segment} found in {name}', file=sys.stderr)
            if skip_unknown_segments: continue
            sys.exit(f'To force execution use {SKIP_UNKNOWN_SEGMENTS_OPT} option.')

        ref_nucl, sample_nucl = pairwise_alignment(references[segment], seq)

        for protein, cds in annotations[segment].items():
            ref_coding, sample_coding = get_coding_sequences(
                ref_nucl, sample_nucl, cds)
            ref_aa = ''.join(translate(ref_coding))
            sample_aa = translate(sample_coding)

            muts_per_sample[sample] += find_mutations(
                ref_aa, sample_aa, mutations[protein])

    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    for sample in muts_per_sample:
        markers_per_sample[sample] = match_markers([mut.name for mut in muts_per_sample[sample]], cur)
    conn.close()

    lines = []
    lines.append('Sample\teffect\tpaper\tsubtype\tfound_mutations\tmarker_mutations')
    for sample in markers_per_sample:
        for marker in markers_per_sample[sample]:
            lst = [sample] + list(marker.values())
            string = '\t'.join(lst)
            lines.append(string)
    out_str = '\n'.join(lines)
    output.encoding = 'utf-8'
    output.write(out_str)


def load_mutations(cur: sqlite3.Cursor) -> Dict[str, List[str]]:
    mutations = defaultdict(list)
    res = cur.execute("""SELECT reference_name, protein_name, name, type, ref_seq, alt_seq, position
                      FROM mutations_characteristics
                      JOIN mutations ON mutations_characteristics.mutation_name = mutations.name""")
    for mut in res:
        mutations[mut[1]].append(Mutation(*mut[2:]))
    return mutations


def load_references(cur: sqlite3.Cursor) -> Dict[str, str]:
    res = cur.execute("SELECT name, sequence FROM 'references'")
    return {name: sequence for name, sequence in res}


def load_annotations(cur: sqlite3.Cursor) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    res = cur.execute("SELECT reference_name, protein_name, start, end FROM 'annotations'")
    ann = defaultdict(lambda: defaultdict(list))
    for ref, prot, start, end in res:
        ann[ref][prot].append((start, end))
    return ann


def match_markers(muts, cur: sqlite3.Cursor):
    muts_str = ','.join([f"'{mut}'" for mut in muts])
    res = cur.execute(f"""
                        WITH marker_mutation_count AS (SELECT markers_mutations.marker_id, group_concat(markers_mutations.mutation_name) AS marker_total_mutations
                                                    FROM markers_mutations GROUP BY markers_mutations.marker_id),
                            markers_tbl AS (SELECT  markers_mutations.marker_id,
                                                group_concat(mutation_name) AS found_mutations,
                                                marker_mutation_count.marker_total_mutations
                                        FROM markers_mutations
                                        JOIN marker_mutation_count ON marker_mutation_count.marker_id = markers_mutations.marker_id
                                        WHERE mutation_name IN (
                                        {muts_str})
                                        GROUP BY markers_mutations.marker_id)

                        SELECT  group_concat(markers_mutations.mutation_name) AS 'Mutations', 
                                effects.name AS 'Effect', 
                                papers.id AS 'Paper', 
                                markers_effects.subtype AS 'Subtype',
                                markers_tbl.found_mutations AS 'Found mutations',
                                markers_tbl.marker_total_mutations AS 'Marker total mutations'
                        FROM markers_effects
                        JOIN markers_mutations ON markers_mutations.marker_id = markers_effects.marker_id
                        JOIN papers ON papers.id = markers_effects.paper_id
                        JOIN effects ON effects.id = markers_effects.effect_id
                        JOIN markers_tbl ON markers_tbl.marker_id = markers_effects.marker_id
                        WHERE markers_mutations.marker_id IN (
                            SELECT markers_tbl.marker_id 
                            FROM markers_tbl)
                        GROUP BY markers_mutations.marker_id, effects.id, papers.id, markers_effects.subtype
                       """)
    found_markers = []
    for _, effect, paper, subtype, found_mutations, marker_mutations in res:
        found_markers.append({
            'effect': effect,
            'paper' : paper,
            'subtype': subtype,
            'found_mutations': found_mutations,
            'marker_mutations': marker_mutations
        })
    return found_markers


def find_mutations(ref_aa, sample_aa, mutations: List[Mutation]):
    found_mutations = []
    for mutation in mutations:
        pos = adjust_position(ref_aa, mutation.pos)
        if mutation.alt in sample_aa[pos]:
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


def parse_name(name: str, pattern: re.Pattern, force: bool) -> Tuple[str, str]:
    '''Get sample and segment information by sequence name'''
    match = pattern.match(name)
    try:
        sample = match.groupdict().get('sample', match.group(1))
        seg = match.groupdict().get('segment', match.group(2))
    except (IndexError, AttributeError):
        print(f'Failed to parse "{name}" with pattern "{pattern.pattern}".', file=sys.stderr)
        if force: return None, None
        sys.exit(f'To force execution and skip this sequence use {SKIP_UNMATCH_NAMES_OPT} option.')
    else:
        return sample, seg


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


def get_codon(seq: List[str], start: int, is_first: bool) -> List[str]:
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
    return ''.join(sorted(set(aas)))


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
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'U': ['T'], '-': ['-'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}


if __name__ == '__main__':
    main()
