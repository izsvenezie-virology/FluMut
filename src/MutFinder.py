#! /usr/bin/env python

import csv
import itertools
import re
import sys
from collections import defaultdict
from io import TextIOWrapper
from typing import Dict, Generator, List, Optional, Tuple
from importlib.resources import files
import sqlite3
import click
from Bio.Align import PairwiseAligner
from click import File
import OutputFormatter
from DataClass import Mutation

PRINT_ALIGNMENT = False
SKIP_UNMATCH_NAMES_OPT = '--skip-unmatch-names'
SKIP_UNKNOWN_SEGMENTS_OPT = '--skip-unknown-segments'
__version__ = '0.2.3'
__author__ = 'Edoardo Giussani'
__contact__ = 'egiussani@izsvenezie.it'


@click.command()
@click.version_option(__version__, '-v', '--version', message=f'%(prog)s, version %(version)s, by {__author__} ({__contact__})')
@click.option(SKIP_UNMATCH_NAMES_OPT, is_flag=True, default=False, help='Skips sequences with name that does not match the pattern')
@click.option(SKIP_UNKNOWN_SEGMENTS_OPT, is_flag=True, default=False, help='Skips sequences with name that does not match the pattern')
@click.option('-s', '--strict', is_flag=True, help='Reports only markers where all mutations are found in sample')
@click.option('-n', '--name-regex', type=str, default=r'(?P<sample>.+)_(?P<segment>.+)', show_default=True, help='Regular expression to parse sequence name')
@click.option('-D', '--db-file', type=str, default=files('data').joinpath('mutfinderDB.sqlite'), help='Source database')
@click.option('-t', '--tabular-output', type=File('w'), default='-', help='The output file [default: stdout]')
@click.option('-m', '--matrix-output', type=File('w'), default=None, help='Report of sequences found in each mutation')
@click.argument('samples-fasta', type=File('r'))
def main(name_regex: str, tabular_output: File, samples_fasta: File, db_file: str, matrix_output: File,
         strict: bool, skip_unmatch_names: bool, skip_unknown_segments: bool) -> None:
    '''
    Search for markers of interest in the SAMPLES-FASTA file.
    '''

    # Initialization
    pattern = re.compile(name_regex)
    
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    segments = load_segments(cur)
    mutations = load_mutations(cur)
    annotations = load_annotations(cur)
    conn.close()

    muts_per_sample: Dict[str, List] = defaultdict(list)
    markers_per_sample: Dict[str, List] = defaultdict(list)

    # Per sequence analysis
    for name, seq in read_fasta(samples_fasta):
        sample,  segment = parse_name(name, pattern, skip_unmatch_names)
        if sample is None or segment is None:
            continue
        if segment not in segments:
            print(f'Unknown segment {segment} found in {name}', file=sys.stderr)
            if skip_unknown_segments: continue
            sys.exit(f'To force execution use {SKIP_UNKNOWN_SEGMENTS_OPT} option.')

        reference_name, reference_sequence = select_reference(segments[segment], seq)
        ref_nucl, sample_nucl = pairwise_alignment(reference_sequence, seq)

        for protein, cds in annotations[reference_name].items():
            ref_coding, sample_coding = get_coding_sequences(
                ref_nucl, sample_nucl, cds)
            ref_aa = ''.join(translate(ref_coding))
            sample_aa = translate(sample_coding)

            muts_per_sample[sample] += find_mutations(
                ref_aa, sample_aa, sample, mutations[protein])

    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    for sample in muts_per_sample:
        markers_per_sample[sample] = match_markers([mut.name for mut in muts_per_sample[sample]], cur, strict)
    conn.close()

    if matrix_output:
        header, data = OutputFormatter.matrix_format(itertools.chain.from_iterable(mutations.values()))
        OutputFormatter.write_csv(tabular_output, header, data)

    if tabular_output:
        OutputFormatter.tabular_output(tabular_output, markers_per_sample)


def load_mutations(cur: sqlite3.Cursor) -> Dict[str, List[Mutation]]:
    mutations = defaultdict(list)
    res = cur.execute("""SELECT reference_name, protein_name, name, type, ref_seq, alt_seq, position
                      FROM mutations_characteristics
                      JOIN mutations ON mutations_characteristics.mutation_name = mutations.name""")
    for mut in res:
        mutations[mut[1]].append(Mutation(*mut[2:]))
    return mutations

def load_segments(cur: sqlite3.Cursor) -> Dict[str, Dict[str, str]]:
    res = cur.execute("SELECT segment_name, name, sequence FROM 'references'")
    segments = defaultdict(dict)
    for segment, name, sequence in res:
        segments[segment][name] = sequence
    return segments


def load_annotations(cur: sqlite3.Cursor) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    res = cur.execute("SELECT reference_name, protein_name, start, end FROM 'annotations'")
    ann = defaultdict(lambda: defaultdict(list))
    for ref, prot, start, end in res:
        ann[ref][prot].append((start, end))
    return ann


def match_markers(muts: List[Mutation], cur: sqlite3.Cursor, strict: str) -> List[Dict[str, str]]:
    muts_str = ','.join([f"'{mut}'" for mut in muts])
    res = cur.execute(f"""
    WITH markers_tbl AS (SELECT marker_id,
                                group_concat(mutation_name) AS found_mutations,
                                count(mutation_name) AS found_mutations_count
                            FROM markers_mutations
                            WHERE mutation_name IN ({muts_str})
                            GROUP BY markers_mutations.marker_id)

    SELECT  markers_summary.all_mutations AS 'Marker mutations',
            markers_tbl.found_mutations AS 'Found mutations',
            markers_effects.effect_name AS 'Effect', 
            markers_effects.paper_id AS 'Paper', 
            markers_effects.subtype AS 'Subtype'
    FROM markers_effects
    JOIN markers_tbl ON markers_tbl.marker_id = markers_effects.marker_id
    JOIN markers_summary ON markers_summary.marker_id = markers_effects.marker_id
    WHERE markers_effects.marker_id IN (
        SELECT markers_tbl.marker_id 
        FROM markers_tbl) { 'AND markers_summary.all_mutations_count = markers_tbl.found_mutations_count' if strict else '' }
    GROUP BY markers_effects.marker_id, markers_effects.effect_name, markers_effects.paper_id, markers_effects.subtype
    """)
    found_markers = []
    for marker_mutations, found_mutations,  effect, paper, subtype in res:
        found_markers.append({
            'marker_mutations': marker_mutations,
            'found_mutations': found_mutations,
            'effect': effect,
            'paper' : paper,
            'subtype': subtype
        })
    return found_markers


def select_reference(references: Dict[str, str], ref_seq: str) -> Tuple[str, str]:
    if len(references) > 1:
        NotImplementedError('Selection for reference from segments with more than one is not yet implemented')
    (name, sequence), = references.items()
    return name, sequence

def find_mutations(ref_aa: str, sample_aa: List[str], sample_name: str, mutations: List[Mutation]):
    found_mutations = []
    for mutation in mutations:
        pos = adjust_position(ref_aa, mutation.pos)
        mutation.samples[sample_name] = sample_aa[pos]
        if mutation.alt in sample_aa[pos]:
            mutation.found = True
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


def read_fasta(fasta_file: TextIOWrapper) -> Generator[str, None, None]:
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


def find_next_nucl(seq: List[str], start: int) -> Optional[int]:
    '''Returns the position of the next non deleted nucleotide'''
    for i in range(start + 3, len(seq)):
        if not seq[i] == '-':
            return i
    return None


def get_coding_sequences(ref_seq: str, sample_seq: str, cds: List[Tuple[int, int]]) -> Tuple[str, str]:
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
