#! /usr/bin/env python

import sqlite3
from collections import defaultdict
from itertools import chain

DB_FILE = 'flumut_db.sqlite'

def to_dict(cursor: sqlite3.Cursor, row: sqlite3.Row):
    result = {}
    for idx, col in enumerate(cursor.description):
        result[col[0]] = row[idx]
    return result


def main():
    db = sqlite3.connect(DB_FILE)
    cursor = db.cursor()
    cursor.row_factory = to_dict

    rows = cursor.execute('''
        SELECT  markers_summary.all_mutations AS 'marker',
                markers_effects.effect_name AS 'effect',
                markers_effects.subtype AS 'subtype',
                markers_effects.paper_id AS 'paper',
                papers.doi AS 'doi',
                papers.web_address AS 'url'
        FROM markers_effects
        JOIN markers_summary ON markers_summary.marker_id = markers_effects.marker_id
        JOIN papers ON papers.id = markers_effects.paper_id
        ''').fetchall()
    markers = defaultdict(lambda: defaultdict(list))

    for row in rows:
        link = row['url']
        if row['doi']:
             link = f'https://doi.org/{row["doi"]}'
        markers[(row['marker'], row['effect'])][row['subtype']].append(f'[{row["paper"]}]({link})')

    md_table = """| Marker | Effect | Subtypes | Literature |\n| ------ | ------ | -------- | ---------- |\n"""
    for (marker, effect), subtypes in markers.items():
        literature = {}
        for index, paper in enumerate(set(list(chain(*subtypes.values())))):
            literature[paper] = index + 1
        literature_str = '<br>'.join([f'{id}. {name}' for name, id in literature.items()])

        subtypes_str = []
        for subtype, papers in subtypes.items():
            ids = ', '.join(map(str, sorted([literature[paper] for paper in papers])))
            subtypes_str.append(f'{subtype} <sup>{ids}</sup>')
        md_table += f'| {marker.replace(",", ", ")} | {effect} | {"<br>".join(subtypes_str)} | {literature_str} |\n'

    with open('_includes/markers.md', 'w', encoding='utf-8') as f:
        print(md_table, file=f)

if __name__ == '__main__':
    main()
