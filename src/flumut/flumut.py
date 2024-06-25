#! /usr/bin/env python

import re
import sys
import click
import operator
import itertools

from io import TextIOWrapper
from click import File
from typing import Dict, Generator, List, Optional, Tuple

from collections import defaultdict
from importlib_resources import files
from Bio.Align import PairwiseAligner

from flumut.DbReader import close_connection, execute_query, open_connection, to_dict, update_db
from flumut import OutputFormatter
from flumut.DataClass import Mutation, Sample

PRINT_ALIGNMENT = False
SKIP_UNMATCH_NAMES_OPT = '--skip-unmatch-names'
SKIP_UNKNOWN_SEGMENTS_OPT = '--skip-unknown-segments'
DB_FILE = files('flumutdata').joinpath('flumut_db.sqlite')
__version__ = '0.5.3'
__author__ = 'Edoardo Giussani'
__contact__ = 'egiussani@izsvenezie.it'

def update(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    update_db(DB_FILE)
    ctx.exit()

@click.command()
@click.help_option('-h', '--help')
@click.version_option(__version__, '-v', '--version', message=f'%(prog)s, version %(version)s, by {__author__} ({__contact__})')
@click.option('--update', is_flag=True, callback=update, expose_value=False, is_eager=True, help='Update the database to the latest version and exit.')
@click.option(SKIP_UNMATCH_NAMES_OPT, is_flag=True, default=False, help='Skip sequences with name that does not match the regular expression pattern.')
@click.option(SKIP_UNKNOWN_SEGMENTS_OPT, is_flag=True, default=False, help='Skip sequences with segment not present in the database.')
@click.option('-r', '--relaxed', is_flag=True, help='Report markers of which at least one mutation is found.')
@click.option('-n', '--name-regex', type=str, default=r'(?P<sample>.+)_(?P<segment>.+)', show_default=True, help='Set regular expression to parse sequence name.')
@click.option('-D', '--db-file', type=str, default=DB_FILE, help='Set source database.')
@click.option('-m', '--markers-output', type=File('w', 'utf-8'), default=None, help='TSV markers output file.')
@click.option('-M', '--mutations-output', type=File('w', 'utf-8'), default=None, help='TSV mutations output file.')
@click.option('-l', '--literature-output', type=File('w', 'utf-8'), default=None, help='TSV literature output file.')
@click.option('-x', '--excel-output', type=str, default=None, help='Excel complete report.')
@click.argument('fasta-file', type=File('r'))
def main(name_regex: str, fasta_file: File, db_file: str, 
         markers_output: File, mutations_output: File, literature_output: File, excel_output: str,
         relaxed: bool, skip_unmatch_names: bool, skip_unknown_segments: bool) -> None:
    '''
    Find markers of zoonotic interest in H5N1 avian influenza viruses.
    '''

    # Initialization
    samples: Dict[str, Sample] = {}
    pattern = re.compile(name_regex)

    open_connection(db_file)
    segments = load_segments()
    mutations = load_mutations()
    annotations = load_annotations()
    close_connection()

    # Per sequence analysis
    for name, seq in read_fasta(fasta_file):
        sample,  segment = parse_name(name, pattern, skip_unmatch_names)
        if sample is None or segment is None:
            continue
        if segment not in segments:
            print(f'Unknown segment {segment} found in {name}', file=sys.stderr)
            if skip_unknown_segments: continue
            sys.exit(f'To force execution use {SKIP_UNKNOWN_SEGMENTS_OPT} option.')

        if sample not in samples:
            samples[sample] = Sample(sample)

        reference_name, reference_sequence = select_reference(segments[segment], seq)
        ref_align, sample_align = pairwise_alignment(reference_sequence, seq)

        for protein, cds in annotations[reference_name].items():
            ref_cds, sample_cds = get_cds(ref_align, sample_align, cds)
            ref_aa = ''.join(translate(ref_cds))
            sample_aa = translate(sample_cds)

            samples[sample].mutations += find_mutations(
                ref_aa, sample_aa, sample, mutations[protein])

    open_connection(db_file)
    for sample in samples.values():
        sample.markers = match_markers(sample.mutations, relaxed)
    papers = load_papers()
    close_connection()
    found_mutations = list(itertools.chain.from_iterable(mutations.values()))
    found_mutations.sort(key=operator.attrgetter('protein', 'pos', 'alt'))

    # Outputs
    if mutations_output:
        header, data = OutputFormatter.mutations_dict(found_mutations)
        OutputFormatter.write_csv(mutations_output, header, data)

    if markers_output:
        header, data = OutputFormatter.markers_dict(samples.values())
        OutputFormatter.write_csv(markers_output, header, data)
    
    if literature_output:
        header, data = OutputFormatter.papers_dict(papers)
        OutputFormatter.write_csv(literature_output, header, data)

    if excel_output:
        wb = OutputFormatter.get_workbook(excel_output.endswith('.xlsm'))
        header, data = OutputFormatter.mutations_dict(found_mutations)
        wb = OutputFormatter.write_excel_sheet(wb, 'Mutations', header, data)
        header, data = OutputFormatter.markers_dict(samples.values())
        wb = OutputFormatter.write_excel_sheet(wb, 'Markers', header, data)
        header, data = OutputFormatter.papers_dict(papers)
        wb = OutputFormatter.write_excel_sheet(wb, 'Literature', header, data)
        wb = OutputFormatter.save_workbook(wb, excel_output)


def load_mutations() -> Dict[str, List[Mutation]]:
    mutations = defaultdict(list)
    res = execute_query(""" SELECT reference_name, protein_name, name, type, ref_seq, alt_seq, position
                            FROM mutation_mappings
                            JOIN mutations ON mutation_mappings.mutation_name = mutations.name""")
    for mut in res:
        mutations[mut[1]].append(Mutation(*mut[2:]))
    return mutations

def load_segments() -> Dict[str, Dict[str, str]]:
    res = execute_query("SELECT segment_name, name, sequence FROM 'references'")
    segments = defaultdict(dict)
    for segment, name, sequence in res:
        segments[segment][name] = sequence
    return segments


def load_annotations() -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    res = execute_query("SELECT reference_name, protein_name, start, end FROM 'annotations'")
    ann = defaultdict(lambda: defaultdict(list))
    for ref, prot, start, end in res:
        ann[ref][prot].append((start, end))
    return ann

def load_papers() -> List[Dict[str, str]]:
    return execute_query("""SELECT  id AS 'Short name',
                                    title AS 'Title',
                                    authors AS 'Authors',
                                    year AS 'Year',
                                    journal AS 'Journal',
                                    web_address AS 'Link',
                                    doi AS 'DOI'
                            FROM papers""", to_dict).fetchall()


def match_markers(muts: List[Mutation], relaxed: bool) -> List[Dict[str, str]]:
    muts_str = ','.join([f"'{mut.name}'" for mut in muts])
    res = execute_query(f"""
    WITH markers_tbl AS (SELECT marker_id,
                                group_concat(mutation_name) AS found_mutations,
                                count(mutation_name) AS found_mutations_count
                            FROM markers_mutations
                            WHERE mutation_name IN ({muts_str})
                            GROUP BY markers_mutations.marker_id)

    SELECT  markers_summary.all_mutations AS 'Marker',
            markers_tbl.found_mutations AS 'Mutations in your sample',
            markers_effects.effect_name AS 'Effect', 
            group_concat(markers_effects.paper_id, '; ') AS 'Literature', 
            markers_effects.subtype AS 'Subtype'
    FROM markers_effects
    JOIN markers_tbl ON markers_tbl.marker_id = markers_effects.marker_id
    JOIN markers_summary ON markers_summary.marker_id = markers_effects.marker_id
    WHERE markers_effects.marker_id IN (
        SELECT markers_tbl.marker_id 
        FROM markers_tbl) { 'AND markers_summary.all_mutations_count = markers_tbl.found_mutations_count' if not relaxed else '' }
    GROUP BY markers_effects.marker_id, markers_effects.effect_name, markers_effects.subtype
    """, to_dict)
    return res.fetchall()


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
        next_nucl = find_next_nucl(seq, start)
        if not next_nucl:
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


def get_cds(ref_seq: str, sample_seq: str, cds: List[Tuple[int, int]]) -> Tuple[str, str]:
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
        dashes = ref_seq.count('-', 0, adj_pos + 1)
        adj_pos = pos + dashes
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
