import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

from flumut.flumutdb.initializer import initialize
from flumut.flumutdb.models import (
    Annotation,
    DbVersion,
    Effect,
    Evidence,
    Host,
    Mapping,
    Marker,
    MarkerMutation,
    Mutation,
    MutationType,
    Paper,
    Protein,
    Reference,
    Segment,
    Subtype,
    Target,
)

BASE = Path(__file__).parent

# Download from https://github.com/izsvenezie-virology/FluMutDB/releases/latest/download/flumut_db.sqlite
INPUT_DB = BASE / 'data' / 'flumut_db.sqlite'
OUTPUT_DB = BASE / 'flumut_db_v7.sqlite'
NEW_REFERENCE_PATH = BASE / 'data' / 'new_refs.json'

CONNECTION = sqlite3.connect(INPUT_DB)
CURSOR = CONNECTION.cursor()

SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
PROTEIN_NAMES = {
    'NS-1': 'NS1',
    'NS-2': 'NS2',
    'HA1-5': 'HA1',
    'HA2-5': 'HA2',
    'NA-1': 'NA',
}
REFERENCE_NAMES = {
    'HA': 'H5',
    'NA': 'N1',
    'NS': 'NS allele A',
}
REFERENCE_SOURCES = {
    'PB2': 'EPI2414998',
    'PB1': 'EPI2414999',
    'PA': 'EPI2414997',
    'H5': 'EPI2415001',
    'NP': 'EPI2415015',
    'N1': 'EPI2415000',
    'MP': 'EPI2414996',
    'NS allele A': 'EPI2414995',
}
EFFECT_NAMES = {
    'Enhanced binding affinity to mammalian cells': ('Enhanced binding affinity in mammalian cells', None),
    'Reduced inhibition to Zanamiriv': ('Reduced inhibition to Zanamivir', None),
    'Evade human BTN3A3 (inhibitor of avian influenza A viruses replication)': (
        'Evade human BTN3A3',
        'Inhibitor of avian influenza A viruses replication',
    ),
    'Evade human stimulator of interferon genes (STING) (inhibitor of avian influenza A viruses replication)': (
        'Evade human stimulator of interferon genes (STING)',
        'Inhibitor of avian influenza A viruses replication',
    ),
    'Conferred Amantidine resistance': ('Conferred resistance to Amantidine', None),
    'Dual α2–3 and α2–6 binding': ('Dual α2-3 and α2-6 binding', None),
    'Increased pseudovirus binding to α2-6': ('Increased binding to α2-6', None),
    'Loss of binding to α2–3': ('Loss of binding to α2-3', None),
    'No binding to α2–6': ('No binding to α2-6', None),
    'Transmitted via aerosol among ferrets': ('Transmitted via aerosol in ferrets', None),
    'Transmitted via aerosol among guinea pigs': ('Transmitted via aerosol in guinea pigs', None),
    'Normal susceptibility to baloxavir (≤3-fold increase)': ('Normal susceptibility to baloxavir', '≤3-fold increase'),
    'Reduced susceptibility to baloxavir (>3-fold increase) ': ('Reduced susceptibility to baloxavir', '>3-fold increase)'),
}
HOST_NAMES = {
    'host': '',
    'chickens (but not ducks)': 'chickens',
    'ferret': 'ferrets',
}
TARGET_NAMES = {
    'amantadine': 'Amantadine',
    'baloxavir': 'Baloxavir',
    'favipiravir': 'Favipiravir',
    'laninamivir': 'Laninamivir',
    'oseltamivir': 'Oseltamivir',
    'peramivir': 'Peramivir',
    'rimantadine': 'Rimantadine',
    'zanamivir': 'Zanamivir',
}
SUBTYPE_NOTES = {
    '(human isolate)': 'human isolate',
    '(human isolate) backbone with H9N2 HA': 'human isolate backbone with H9N2 HA',
    'with all internal genes from H7N9': 'all internal genes from H7N9',
    '(avian)': '',
}

HOST_RE = re.compile(r'^(.+?)\s+in\s+(.+)$')
TARGET_RE = re.compile(
    r'^(.+?)\s+to\s+(' + '|'.join(map(re.escape, TARGET_NAMES)) + r')$',
    re.IGNORECASE,
)


def load_numbering(path: Path) -> dict[int, dict[str, int]]:
    result = defaultdict(dict)
    with open(path, 'r', encoding='utf-8') as f:
        content = f.readlines()
    header = content[0].strip().split('\t')
    for line in content[1:]:
        parts = list(map(int, line.strip().split('\t')))
        for i, pos in enumerate(parts):
            result[int(parts[0])][header[i]] = int(pos)
    return result


HA1_NUMBERING = load_numbering(BASE / 'data' / 'ha1_numbering.tsv')
HA2_NUMBERING = load_numbering(BASE / 'data' / 'ha2_numbering.tsv')
NA_NUMBERING = load_numbering(BASE / 'data' / 'na_numbering.tsv')


def sanitize_protein_name(name: str) -> str:
    return PROTEIN_NAMES.get(name, name).strip()


def sanitize_reference_name(name: str) -> str:
    return REFERENCE_NAMES.get(name, name).strip()


def sanitize_effect_name(name: str) -> tuple[str, str | None]:
    return EFFECT_NAMES.get(name, (name.strip(), None))


def sanitize_host_name(name: str) -> str:
    return HOST_NAMES.get(name, name).lower().strip()


def sanitize_target_name(name: str) -> str:
    return TARGET_NAMES.get(name.casefold(), '').strip()


def sanitize_subtype_note(note: str) -> str:
    return SUBTYPE_NOTES.get(note, note.strip())


def sanitize_mutation_name(name: str) -> str:
    protein, mutation = name.split(':')
    protein = sanitize_protein_name(protein)
    return f'{protein}:{mutation[1:].strip()}'


def sanitize_proteins(proteins: dict[str, Protein]) -> dict[str, Protein]:
    result = {}
    for name, protein in proteins.items():
        result[name] = protein
        result[sanitize_protein_name(name)] = protein
    return result


def sanitize_references(references: dict[str, Reference]) -> dict[str, Reference]:
    result = {}
    for name, reference in references.items():
        result[name] = reference
        result[sanitize_reference_name(name)] = reference
    return result


def migrate_segments() -> dict[str, Segment]:
    Segment.create_table()
    result = {}
    for i, name in enumerate(SEGMENTS):
        segment = Segment.create(name=name, order=i + 1)
        result[name] = segment
    return result


def migrate_proteins(segments: dict[str, Segment]) -> dict[str, Protein]:
    Protein.create_table()
    result = {}
    indexes = defaultdict(int)

    CURSOR.execute('SELECT name, segment_name FROM proteins')
    for old_name, segment_name in CURSOR.fetchall():
        name = sanitize_protein_name(old_name)
        indexes[segment_name] += 1
        protein = Protein.create(name=name, segment=segments[segment_name], order=indexes[segment_name])
        result[old_name] = protein
    return result


def migrate_references(segments: dict[str, Segment]) -> dict[str, Reference]:
    Reference.create_table()
    result = {}
    indexes = defaultdict(int)

    CURSOR.execute('SELECT name, segment_name, sequence FROM "references"')
    for old_name, segment_name, sequence in CURSOR.fetchall():
        name: str = sanitize_reference_name(old_name)
        if segment_name in ('HA', 'NA'):
            index = int(name[1:])
        else:
            indexes[segment_name] += 1
            index = indexes[segment_name]
        reference = Reference.create(
            name=name, segment=segments[segment_name], sequence=sequence.strip(), order=index, source=REFERENCE_SOURCES.get(name)
        )
        result[old_name] = reference
    return result


def migrate_annotations(proteins: dict[str, Protein], references: dict[str, Reference]) -> None:
    Annotation.create_table()
    indexes = defaultdict(int)

    CURSOR.execute('SELECT start, end, protein_name, reference_name FROM annotations')

    for start, end, protein_name, reference_name in CURSOR.fetchall():
        indexes[(protein_name, reference_name)] += 1
        Annotation.create(
            start=int(start),
            end=int(end),
            protein=proteins[protein_name],
            reference=references[reference_name],
            order=indexes[(protein_name, reference_name)],
        )


def import_new_references(source: Path, segments: dict[str, Segment], proteins: dict[str, Protein]) -> dict[str, Reference]:
    data = json.loads(source.read_text('utf-8'))
    result = {}
    for segment_name, references in data.items():
        for ref_data in references:
            if segment_name in ('HA', 'NA'):
                index = int(ref_data['name'][1:])
            else:
                index = 2  # NS allele B
            reference = Reference.create(
                name=ref_data['name'],
                segment=segments[segment_name],
                sequence=ref_data['sequence'],
                order=index,
                source=ref_data['source'],
            )
            result[reference.name] = reference
            indexes = defaultdict(int)
            for annotation in ref_data['annotations']:
                protein_name = annotation['protein']
                indexes[(protein_name, reference.name)] += 1
                Annotation.create(
                    start=int(annotation['start']),
                    end=int(annotation['end']),
                    protein=proteins[protein_name],
                    reference=reference,
                    order=indexes[(protein_name, reference.name)],
                )
    return result


def migrate_mutations(proteins: dict[str, Protein], references: dict[str, Reference]) -> dict[str, Mutation]:
    Mutation.create_table()
    Mapping.create_table()
    indexes = defaultdict(int)
    result = {}
    existing_names = {}

    CURSOR.execute('SELECT mutation_name, reference_name, position, alt_seq FROM mutation_mappings ORDER BY position')

    for old_name, old_reference_name, position, alt_seq in CURSOR.fetchall():
        mutation_name = sanitize_mutation_name(old_name)
        if mutation_name in existing_names:
            result[old_name] = existing_names[mutation_name]
            continue

        protein = proteins[sanitize_protein_name(old_name.split(':')[0])]
        reference = references[old_reference_name]
        indexes[(protein, reference)] += 1

        if protein.segment.name == 'HA':
            mutation = create_ha_mutations(protein, position, alt_seq, indexes[(protein, reference)], references, existing_names)
        elif protein.segment.name == 'NA':
            mutation = create_na_mutations(mutation_name, protein, position, alt_seq, indexes[(protein, reference)], references)
        elif protein.segment.name == 'NS':
            mutation = create_ns_mutations(mutation_name, protein, position, alt_seq, indexes[(protein, reference)], references)
        else:
            mutation = Mutation.create(
                name=mutation_name, type=MutationType.SNP.value, protein=protein, order=indexes[(protein, reference)]
            )
            Mapping.create(mutation=mutation, reference=reference, position=position, alteration=alt_seq)
        result[old_name] = mutation
        existing_names[mutation.name] = mutation
    return result


def create_ha_mutations(
    protein: Protein, position: int, alt_seq: str, index: int, references: dict[str, Reference], existing_names: dict[str, Mutation]
) -> Mutation:
    mutation = Mutation.create(name='tmp', type=MutationType.SNP.value, protein=protein, order=index)
    if protein.name == 'HA1':
        positions = HA1_NUMBERING[position]
    else:
        positions = HA2_NUMBERING[position]
    for reference_name, pos in positions.items():
        if not pos:
            continue
        reference = references[reference_name]
        mapping = Mapping.create(mutation=mutation, reference=reference, position=pos, alteration=alt_seq)
        if reference_name == 'H3':
            mutation.name = f'{mutation.protein.name}:{mapping.position}{mapping.alteration}'
        elif reference_name == 'H5':
            h5_mapping = mapping
    if mutation.name == 'tmp':
        mutation.name = f'{mutation.protein.name}:{h5_mapping.position}{h5_mapping.alteration}'  # pyright: ignore[reportPossiblyUnboundVariable]
        mutation.notes = 'H5 numbering'
        print(f'Using H5 numbering in mutation {mutation.get_id()}: {mutation.name}')
    if mutation.name in existing_names:
        mutation.delete_instance()
        print(f'Duplicated mutation name {mutation.name}: using mutation {existing_names[mutation.name].get_id()} instead')
        return existing_names[mutation.name]
    mutation.save()
    return mutation


def create_na_mutations(
    mutation_name: str, protein: Protein, position: int, alt_seq: str, index: int, references: dict[str, Reference]
) -> Mutation:
    mutation = Mutation.create(name=mutation_name, type=MutationType.SNP.value, protein=protein, order=index)
    positions = NA_NUMBERING[position]
    for reference_name, pos in positions.items():
        if not pos:
            continue
        reference = references[reference_name]
        Mapping.create(mutation=mutation, reference=reference, position=pos, alteration=alt_seq)
    mutation.save()
    return mutation


def create_ns_mutations(
    mutation_name: str, protein: Protein, position: int, alt_seq: str, index: int, references: dict[str, Reference]
) -> Mutation:
    mutation = Mutation.create(name=mutation_name, type=MutationType.SNP.value, protein=protein, order=index)
    Mapping.create(mutation=mutation, reference=references['NS allele A'], position=position, alteration=alt_seq)
    Mapping.create(mutation=mutation, reference=references['NS allele B'], position=position, alteration=alt_seq)
    return mutation


def migrate_effects() -> dict[str, tuple[Effect, Host | None, Target | None]]:
    Effect.create_table()
    Host.create_table()
    Target.create_table()
    result = {}

    CURSOR.execute('SELECT name FROM effects')

    for (name,) in CURSOR.fetchall():
        sanitized_name, notes = sanitize_effect_name(name)
        effect_name, host_name, target_name = parse_effect(sanitized_name)
        effect, created = Effect.get_or_create(name=effect_name, defaults={'notes': notes})
        if created:
            print(f'Effect "{name}" is merged into effect {effect.get_id()}: {effect.name}')
        host, target = None, None
        if host_name := sanitize_host_name(host_name):
            host, _ = Host.get_or_create(name=host_name)
            print(f'Extracted host "{host_name}" from "{name}"')
        if target_name := sanitize_target_name(target_name):
            target, _ = Target.get_or_create(name=target_name)
            print(f'Extracted target "{target_name}" from "{name}"')
        result[name] = (effect, host, target)
    return result


def parse_effect(name: str) -> tuple[str, str, str]:
    if m := HOST_RE.match(name):
        return m.group(1), m.group(2), ''
    if m := TARGET_RE.match(name):
        return m.group(1), '', m.group(2)
    return name, '', ''


def migrate_subtypes() -> dict[str, tuple[Subtype, str]]:
    Subtype.create_table()
    result = {}

    CURSOR.execute('SELECT DISTINCT subtype FROM markers_effects')
    for (name,) in CURSOR.fetchall():
        new_name, notes = parse_subtype(name)
        subtype, _ = Subtype.get_or_create(name=new_name.strip())
        result[name] = (subtype, sanitize_subtype_note(notes))
    return result


def parse_subtype(name: str) -> tuple[str, str]:
    if name.startswith('A('):
        name = name[2 : name.find(')')]
    if ' ' in name:
        subtype, notes = name.split(' ', 1)
        return subtype, notes
    return name, ''


def migrate_papers() -> dict[str, Paper]:
    Paper.create_table()
    result = {}

    CURSOR.execute('SELECT id, title, authors, year, journal, web_address, doi FROM papers')
    for short_name, title, authors, year, journal, url, doi in CURSOR.fetchall():
        values = {
            'short_name': short_name.strip(),
            'title': title.strip(),
            'authors': authors.strip(),
            'year': int(year),
            'journal': journal.strip() or None,
            'url': url.strip() or None,
        }
        paper, created = Paper.get_or_create(doi=doi, defaults=values)
        if not created:
            print(f'Duplicated doi paper {doi}: using {paper.short_name} (id {paper.get_id()}) instead')
        result[short_name] = paper
    return result


def migrate_markers(mutations: dict[str, Mutation]) -> dict[str, Marker]:
    Marker.create_table()
    MarkerMutation.create_table()
    result = {}
    existing_names = {}

    CURSOR.execute('SELECT id, notes FROM markers')

    for id, notes in CURSOR.fetchall():
        marker = Marker.create(name='tmp', notes=notes or None)
        CURSOR.execute(f'SELECT mutation_name FROM markers_mutations WHERE marker_id == {id}')
        for (mutation_name,) in CURSOR.fetchall():
            MarkerMutation.create(marker=marker, mutation=mutations[mutation_name])
        marker.name = '; '.join(sorted(mutation.name for mutation in marker.mutations))

        if not marker.name:
            print(f'Discarding empty marker: {id} ({notes})')
        if marker.name in existing_names:
            print(f'Found duplicated marker: {id} ({notes}), using marker {marker.get_id()} instead')
            result[id] = existing_names[marker.name]
            marker.delete_instance()
            continue
        result[id] = marker
        existing_names[marker.name] = marker
        marker.save()
    return result


def migrate_evidences(
    markers: dict[str, Marker],
    papers: dict[str, Paper],
    effects: dict[str, tuple[Effect, Host | None, Target | None]],
    subtypes: dict[str, tuple[Subtype, str]],
) -> None:
    Evidence.create_table()
    CURSOR.execute('SELECT marker_id, paper_id, effect_name, subtype FROM markers_effects')
    for marker_id, paper_id, effect_name, subtype_name in CURSOR.fetchall():
        effect, host, target = effects[effect_name]
        subtype, notes = subtypes[subtype_name]
        values = {'marker': markers[marker_id], 'paper': papers[paper_id], 'effect': effect, 'subtype': subtype}
        if host:
            values['host'] = host
        if target:
            values['target'] = target
        if notes:
            values['notes'] = f'Subtype: {notes}'
        Evidence.create(**values)


# Create v7 database
if OUTPUT_DB.exists():
    OUTPUT_DB.unlink()

initialize(str(OUTPUT_DB), read_only=False)

DbVersion.create_table()
DbVersion.create(major=7, minor=0, patch=0)
segments = migrate_segments()
proteins = migrate_proteins(segments)
proteins = sanitize_proteins(proteins)
references = migrate_references(segments)
migrate_annotations(proteins, references)
new_references = import_new_references(NEW_REFERENCE_PATH, segments, proteins)
references.update(new_references)
references = sanitize_references(references)
mutations = migrate_mutations(proteins, references)
effects = migrate_effects()
subtypes = migrate_subtypes()
papers = migrate_papers()
markers = migrate_markers(mutations)
migrate_evidences(markers, papers, effects, subtypes)
