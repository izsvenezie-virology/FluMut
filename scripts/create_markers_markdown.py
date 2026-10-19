#!/usr/bin/env python3
"""Render the bundled FluMutDB as the markers table published on the documentation site.

Writes ``docs/_includes/markers.md``, which ``docs/docs/markers.md`` includes and
DataTables turns into a searchable table. The four columns match the widths that
include declares.

Usage::

    python scripts/create_markers_markdown.py
"""

import argparse
from pathlib import Path

from flumut.core.options import DatabaseOptions
from flumut.flumutdb import Marker, open_database

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / 'docs' / '_includes' / 'markers.md'

COLUMNS = ('Marker', 'Effect', 'Subtype', 'Literature')


def paper_link(paper) -> str:
    """A Markdown link to ``paper``, by URL or DOI, falling back to its bare name."""
    target = f'https://doi.org/{paper.doi}' if paper.doi else paper.url or ''
    name = paper.short_name
    return f'[{name}]({target})' if target else name


def marker_sort_key(marker: Marker) -> tuple[int, ...]:
    """Order markers by their first mutation, so the table follows segment order."""
    return min((mutation.sort_key for mutation in marker.mutations), default=(1_000_000,))


def marker_rows(marker: Marker) -> list[tuple[str, str, str, str]]:
    """One row per (effect, subtype) the literature reports for ``marker``."""
    grouped: dict[tuple[str, str], list] = {}
    for evidence in marker.evidences:
        grouped.setdefault((evidence.get_effect_name(), evidence.subtype.name), []).append(evidence.paper)

    rows = []
    for (effect, subtype), papers in grouped.items():
        unique = sorted({paper.short_name: paper for paper in papers}.values(), key=lambda paper: paper.short_name)
        rows.append((marker.name, effect, subtype, '; '.join(paper_link(p) for p in unique)))
    return sorted(rows, key=lambda row: (row[1], row[2]))


def render(markers: list[Marker]) -> str:
    """Render the whole table, markers in segment order."""
    lines = [
        '| ' + ' | '.join(COLUMNS) + ' |',
        '|' + '|'.join(['---'] * len(COLUMNS)) + '|',
    ]
    for marker in sorted(markers, key=marker_sort_key):
        lines += ['| ' + ' | '.join(row) + ' |' for row in marker_rows(marker)]
    return '\n'.join(lines) + '\n'


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--database', default=None, help='Database to read. Defaults to the bundled one.')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT, help='Markdown file to write.')
    args = parser.parse_args()

    open_database(DatabaseOptions(path=args.database, read_only=True))

    markers: list[Marker] = list(Marker.select())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render(markers), encoding='utf-8', newline='\n')


if __name__ == '__main__':
    main()
