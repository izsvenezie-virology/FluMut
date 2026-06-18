"""
Migrates flumut_db.sql (v6 schema) to a new SQLite file (v7 schema).

Reads the v6 SQL dump, applies all structural and data transformations
from the original sql_updates.py pipeline, and writes a clean SQLite
file using the ORM models defined in src/flumutdb/models.py.

Usage:
    python migrations/update_v7/migrate.py

The result is written next to this script as flumut_db_v7.sqlite.
"""

import json
import re
import sqlite3
from collections import defaultdict
from csv import DictReader
from pathlib import Path

from peewee import SqliteDatabase  # noqa: E402

from flumut.core.globals import DATABASE_PROXY
from flumut.flumutdb.models import (  # noqa: E402
    Annotation,
    DbVersion,
    Effect,
    Evidence,
    Host,
    Mapping,
    Marker,
    Mutation,
    Paper,
    Protein,
    Reference,
    Segment,
    Subtype,
)

# ── constants ──────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent

SEQUENCE_IDS = {
    'PB2': 'EPI2414998',
    'PB1': 'EPI2414999',
    'PA': 'EPI2414997',
    'HA': 'EPI2415001',
    'NP': 'EPI2415015',
    'NA': 'EPI2415000',
    'MP': 'EPI2414996',
    'NS': 'EPI2414995',
}

# Old reference name → new reference name
REF_RENAMES = {'HA': 'H5', 'NA': 'N1', 'NS': 'NS allele A'}

# Old protein name → new protein name
PROTEIN_RENAMES = {'NS-1': 'NS1', 'NS-2': 'NS2'}

# Old mutation name prefix → new prefix (applied in order)
MUT_PREFIX_SUBS = [
    ('HA1-5:', 'HA1:'),
    ('HA2-5:', 'HA2:'),
    ('NA-1:', 'NA:'),
    ('NS-1:', 'NS1:'),
    ('NS-2:', 'NS2:'),
]

# ── helpers ────────────────────────────────────────────────────────────────────


def load_numbering(path: Path, primary_col: str) -> dict[str, dict[str, str]]:
    """Return {primary_value: {other_col: value, ...}} from a TSV file."""
    result: dict[str, dict[str, str]] = {}
    with open(path) as f:
        for row in DictReader(f, delimiter='\t'):
            key = row.pop(primary_col)
            result[key] = dict(row)
    return result


def transform_mutation_name(
    old_name: str,
    ha1_numbering: dict[str, dict[str, str]],
    ha2_numbering: dict[str, dict[str, str]],
) -> tuple[str, int | None]:
    """
    Apply the v7 mutation-name pipeline to a single name.

    Steps (mirror sql_updates.py):
      1. Replace subtype prefixes (HA1-5 → HA1, NA-1 → NA, etc.)
      2. Remove the reference amino-acid that precedes the position
         (PROT:R292K → PROT:292K)
      3. For HA mutations, replace H5-numbered position with H3 position;
         keep H5 numbering when no H3 equivalent exists.
      4. Special-case: HA1:129V → HA1:129V (H5 numbering)

    Returns (new_name, default_position).
    """
    name = old_name

    # Step 1 – prefix substitution
    for old_prefix, new_prefix in MUT_PREFIX_SUBS:
        if name.startswith(old_prefix):
            name = new_prefix + name[len(old_prefix) :]
            break

    # Step 2 – remove ref AA (first character after the colon)
    colon = name.index(':')
    name = name[: colon + 1] + name[colon + 2 :]

    # # Step 3 – H3 renumbering for HA mutations
    # if name.startswith('HA'):
    #     colon = name.index(':')
    #     protein = name[:colon]  # "HA1" or "HA2"
    #     pos_str = name[colon + 1 : -1]  # everything between ':' and last char
    #     alt_aa = name[-1]

    #     numbering = ha1_numbering if protein == 'HA1' else ha2_numbering
    #     h3_pos = numbering.get(pos_str, {}).get('H3', '')

    #     if h3_pos:
    #         name = f'{protein}:{h3_pos}{alt_aa}'
    #     # else: keep H5 name as-is (OR IGNORE behaviour)

    # # Step 4 – special rename for the one H5-only HA position
    # if name == 'HA1:129V':
    #     name = 'HA1:129V (H5 numbering)'

    # Extract default_position
    colon = name.index(':')
    pos_match = re.search(r'\d+', name[colon + 1 :])
    default_position = int(pos_match.group()) if pos_match else None

    return name, default_position


# ── main ───────────────────────────────────────────────────────────────────────


def main() -> None:
    sql_path = BASE / 'data' / 'flumut_db.sql'
    output_path = BASE / 'flumut_db_v7.sqlite'

    print(f'Source : {sql_path}')
    print(f'Output : {output_path}')

    # ── 1. Load v6 dump into an in-memory SQLite ──────────────────────────────
    src = sqlite3.connect(':memory:')
    src.row_factory = sqlite3.Row
    src.executescript(sql_path.read_text(encoding='utf-8'))

    # ── 2. Load auxiliary data ────────────────────────────────────────────────
    with open(BASE / 'data' / 'new_refs.json') as f:
        new_refs: dict[str, list[dict]] = json.load(f)

    # Drop subtypes that are already covered by the renamed originals
    excluded = {'H5', 'N1', 'NS allele A'}
    for refs in new_refs.values():
        refs[:] = [r for r in refs if r['subtype'] not in excluded]

    ha1_numbering = load_numbering(BASE / 'data' / 'ha1_numbering.tsv', 'H5')
    ha2_numbering = load_numbering(BASE / 'data' / 'ha2_numbering.tsv', 'H5')
    na_numbering = load_numbering(BASE / 'data' / 'na_numbering.tsv', 'N1')

    # ── 3. Create new DB ──────────────────────────────────────────────────────
    if output_path.exists():
        output_path.unlink()

    new_db = SqliteDatabase(str(output_path), pragmas={'foreign_keys': 1})
    DATABASE_PROXY.initialize(new_db)

    MarkerMutation = Marker.mutations.get_through_model()
    new_db.create_tables(
        [
            Segment,
            Protein,
            Reference,
            Annotation,
            Mutation,
            Mapping,
            Effect,
            Subtype,
            Host,
            Paper,
            Marker,
            MarkerMutation,
            Evidence,
            DbVersion,
        ]
    )

    with new_db.atomic():
        # ── 4. Segments ───────────────────────────────────────────────────────
        seg_map: dict[str, Segment] = {}
        for row in src.execute('SELECT name FROM segments'):
            seg = Segment.create(name=row['name'])
            seg_map[row['name']] = seg

        # ── 5. Proteins (rename NS-1/NS-2) ────────────────────────────────────
        prot_map: dict[str, Protein] = {}  # old name → Protein
        for row in src.execute('SELECT name, segment_name FROM proteins'):
            old_name = row['name']
            prot = Protein.create(
                name=PROTEIN_RENAMES.get(old_name, old_name),
                segment=seg_map[row['segment_name']],
            )
            prot_map[old_name] = prot

        # ── 6. References (rename + sequence_id → source) + new refs ─────────
        ref_map: dict[str, Reference] = {}  # both old & new names → Reference
        for row in src.execute('SELECT name, segment_name, sequence FROM "references"'):
            old_name = row['name']
            new_name = REF_RENAMES.get(old_name, old_name)
            ref = Reference.create(
                name=new_name,
                segment=seg_map[row['segment_name']],
                sequence=row['sequence'],
                source=SEQUENCE_IDS.get(old_name, ''),
            )
            ref_map[old_name] = ref
            ref_map[new_name] = ref

        for seg_name, refs in new_refs.items():
            for ref_data in refs:
                ref = Reference.create(
                    name=ref_data['subtype'],
                    segment=seg_map[seg_name],
                    sequence=ref_data['sequence'],
                    source=ref_data.get('sequence_id', ''),
                )
                ref_map[ref_data['subtype']] = ref

        # ── 7. Annotations (original + new-ref annotations) ───────────────────
        for row in src.execute('SELECT protein_name, reference_name, start, end FROM annotations'):
            ref_name = REF_RENAMES.get(row['reference_name'], row['reference_name'])
            ref = ref_map.get(ref_name)
            prot = prot_map.get(row['protein_name'])
            if ref and prot:
                Annotation.create(protein=prot, reference=ref, start=row['start'], end=row['end'])

        for seg_name, refs in new_refs.items():
            for ref_data in refs:
                ref = ref_map[ref_data['subtype']]
                for ann in ref_data.get('proteins', []):
                    prot = prot_map.get(ann['name'])
                    if prot:
                        Annotation.create(
                            protein=prot,
                            reference=ref,
                            start=ann['start'],
                            end=ann['end'],
                        )

        # ── 8. Mutations (delete NA-5:D199G, apply name transformations) ──────
        mut_map: dict[str, Mutation] = {}  # old name → Mutation
        for row in src.execute("SELECT name, type, protein_name FROM mutations WHERE name != 'NA-5:D199G'"):
            old_name = row['name']
            new_name, default_pos = transform_mutation_name(old_name, ha1_numbering, ha2_numbering)
            # Two old mutations can collapse to the same new name (e.g. D183G and
            # E183G → HA1:183G) because the ref-AA was stripped.  Mirror the SQL
            # OR REPLACE behaviour by re-using an already-created entry.
            mut, _ = Mutation.get_or_create(
                name=new_name,
                defaults={
                    'type': row['type'],
                    'protein': prot_map[row['protein_name']],
                    'default_position': default_pos,
                },
            )
            mut_map[old_name] = mut

        # ── 9. Mappings ───────────────────────────────────────────────────────

        # 9a. Migrate existing mappings (cascade: ref renames already in ref_map)
        for row in src.execute('SELECT mutation_name, reference_name, position, alt_seq FROM mutation_mappings'):
            mut = mut_map.get(row['mutation_name'])
            if mut is None:
                continue  # deleted mutation (NA-5:D199G)
            ref_name = REF_RENAMES.get(row['reference_name'], row['reference_name'])
            ref = ref_map.get(ref_name)
            if ref is None:
                continue
            Mapping.create(
                mutation=mut,
                reference=ref,
                position=row['position'],
                alteration=row['alt_seq'],
            )

        # 9b. New mappings for new references

        # Fetch old per-segment mappings (using original reference names)
        old_ha_maps = list(src.execute("SELECT mutation_name, position, alt_seq FROM mutation_mappings WHERE reference_name = 'HA'"))
        old_na_maps = list(src.execute("SELECT mutation_name, position, alt_seq FROM mutation_mappings WHERE reference_name = 'NA'"))
        old_ns_maps = list(src.execute("SELECT mutation_name, position, alt_seq FROM mutation_mappings WHERE reference_name = 'NS'"))

        # NS allele B: same positions as old NS allele A
        ns_b_ref = ref_map.get('NS allele B')
        if ns_b_ref:
            for row in old_ns_maps:
                mut = mut_map.get(row['mutation_name'])
                if mut:
                    Mapping.create(
                        mutation=mut,
                        reference=ns_b_ref,
                        position=row['position'],
                        alteration=row['alt_seq'],
                    )

        # NA subtypes: convert from N1 numbering
        for ref_data in new_refs.get('NA', []):
            ref = ref_map[ref_data['subtype']]
            for row in old_na_maps:
                mut = mut_map.get(row['mutation_name'])
                if mut is None:
                    continue
                new_pos = na_numbering.get(str(row['position']), {}).get(ref_data['subtype'], '')
                if new_pos:
                    Mapping.create(
                        mutation=mut,
                        reference=ref,
                        position=int(new_pos),
                        alteration=row['alt_seq'],
                    )

        # HA subtypes: convert from H5 numbering
        for ref_data in new_refs.get('HA', []):
            ref = ref_map[ref_data['subtype']]
            for row in old_ha_maps:
                mut = mut_map.get(row['mutation_name'])
                if mut is None:
                    continue
                # Determine HA1 vs HA2 from the old mutation name prefix
                protein_prefix = row['mutation_name'].split(':')[0][:3]  # "HA1" or "HA2"
                numbering = ha1_numbering if protein_prefix == 'HA1' else ha2_numbering
                new_pos = numbering.get(str(row['position']), {}).get(ref_data['subtype'], '')
                if new_pos:
                    Mapping.create(
                        mutation=mut,
                        reference=ref,
                        position=int(new_pos),
                        alteration=row['alt_seq'],
                    )

        # ── 10. Effects + Hosts ───────────────────────────────────────────────
        _HOST_RE = re.compile(r'^(.+?)\s+in\s+(.+)$')
        eff_map: dict[str, Effect] = {}
        host_map: dict[str, Host] = {}
        host_for_effect: dict[str, Host | None] = {}

        for row in src.execute('SELECT name FROM effects'):
            old_name = row['name']
            m = _HOST_RE.match(old_name)
            if m:
                base, host_name = m.group(1), m.group(2)
                if host_name not in host_map:
                    host_map[host_name] = Host.create(name=host_name)
                host_obj = host_map[host_name]
            else:
                base, host_obj = old_name, None
            eff, _ = Effect.get_or_create(name=base)
            eff_map[old_name] = eff
            host_for_effect[old_name] = host_obj

        # ── 11. Papers ────────────────────────────────────────────────────────
        # Map old paper columns to the new schema (id → short_name,
        # web_address → url) while honouring the new UNIQUE constraint on doi:
        #   * Empty DOIs become NULL: many papers legitimately have no DOI, and
        #     SQLite allows duplicate NULLs but not duplicate empty strings.
        #   * A handful of publications were entered twice under one DOI (e.g.
        #     "Smith, 2006" and "Smith, 2006b"). Merge each pair into a single
        #     paper built from the cleaned-up "...b" record, drop the trailing
        #     "b" from its short_name, and route both old ids to it so every
        #     evidence row is preserved.
        paper_rows = list(src.execute('SELECT id, title, authors, year, journal, web_address, doi FROM papers'))

        # Group the ids of publications that were entered more than once under
        # the same DOI, then pick the "...b" record of each group as the keeper.
        ids_by_doi: dict[str, list[str]] = defaultdict(list)
        for row in paper_rows:
            if row['doi']:
                ids_by_doi[row['doi']].append(row['id'])

        alias_to_keeper: dict[str, str] = {}  # duplicate id → keeper id
        for duplicate_ids in ids_by_doi.values():
            if len(duplicate_ids) < 2:
                continue
            keeper = next(pid for pid in duplicate_ids if pid.endswith('b'))
            for pid in duplicate_ids:
                if pid != keeper:
                    alias_to_keeper[pid] = keeper
        keepers = set(alias_to_keeper.values())

        paper_map: dict[str, Paper] = {}  # old id → Paper (aliases reuse their keeper)
        for row in paper_rows:
            old_id = row['id']
            if old_id in alias_to_keeper:
                continue  # built from its keeper instead; linked up below
            paper_map[old_id] = Paper.create(
                short_name=old_id[:-1] if old_id in keepers else old_id,
                title=row['title'],
                authors=row['authors'],
                year=row['year'],
                journal=row['journal'],
                url=row['web_address'],
                doi=row['doi'] or None,
            )

        # Point every duplicate id at the single paper built from its keeper.
        for alias, keeper in alias_to_keeper.items():
            paper_map[alias] = paper_map[keeper]

        # ── 12. Markers ───────────────────────────────────────────────────────
        marker_map: dict[int, Marker] = {}  # old id → Marker
        for row in src.execute('SELECT id, notes FROM markers'):
            marker_map[row['id']] = Marker.create(notes=row['notes'])

        # ── 13. Marker↔Mutation M2M ───────────────────────────────────────────
        for row in src.execute('SELECT marker_id, mutation_name FROM markers_mutations'):
            mut = mut_map.get(row['mutation_name'])
            marker = marker_map.get(row['marker_id'])
            if mut and marker:
                marker.mutations.add(mut)

        # ── 14. Evidence ──────────────────────────────────────────────────────
        subtype_map: dict[str, Subtype] = {}

        for row in src.execute('SELECT marker_id, paper_id, effect_name, subtype FROM markers_effects'):
            marker = marker_map.get(row['marker_id'])
            paper = paper_map.get(row['paper_id'])
            eff = eff_map.get(row['effect_name'])
            if not (marker and paper and eff):
                continue
            sub_name = row['subtype']
            if sub_name not in subtype_map:
                subtype_map[sub_name], _ = Subtype.get_or_create(name=sub_name)
            Evidence.create(
                marker=marker,
                paper=paper,
                effect=eff,
                subtype=subtype_map[sub_name],
                host=host_for_effect.get(row['effect_name']),
            )

        # ── 15. DbVersion ─────────────────────────────────────────────────────
        DbVersion.create(major=7, minor=0, patch=0)

    src.close()
    print('Done.')


if __name__ == '__main__':
    main()
