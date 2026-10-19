"""End-to-end tests asserting properties that hold for any database version.

Which markers get reported depends entirely on FluMutDB content, so nothing
here asserts a marker name or a row count: a new marker or a new piece of
evidence must not turn these red.
"""

from openpyxl import load_workbook

#: The only sample in ``single_sample.fa``, as parsed by the default name regex.
SAMPLE_ID = 'A/black-headed_gull/Italy/23VIR1099-4/2023_H5N1_2023-02-03'

MARKERS_COLUMNS = [
    'Sample',
    'Marker',
    'Mutations in your sample',
    'Effect',
    'Subtype',
    'Literature',
    'Marker (H5 numbering)',
    'Notes',
]
LITERATURE_COLUMNS = ['Short name', 'Title', 'Authors', 'Year', 'Journal', 'Link', 'DOI']


def test_run_produces_wellformed_outputs(default_run) -> None:
    assert default_run.exception is None
    assert default_run.exit_code == 0

    assert default_run.columns('markers') == MARKERS_COLUMNS
    assert default_run.columns('literature') == LITERATURE_COLUMNS
    assert default_run.columns('mutations')[0] == 'Sample'
    assert {row['Sample'] for row in default_run.rows('mutations')} == {SAMPLE_ID}

    for row in default_run.rows('markers'):
        assert row['Marker'] and row['Effect'] and row['Subtype'], row


def test_run_without_any_output_fails(flumut) -> None:
    assert flumut(outputs=()).exit_code != 0


def test_outputs_are_mutually_consistent(default_run) -> None:
    """Nothing referenced in the markers output is missing from the others."""
    marker_rows = default_run.rows('markers')

    published = {row['Short name'] for row in default_run.rows('literature')}
    cited = {name for row in marker_rows for name in row['Literature'].split('; ') if name}
    assert cited <= published, f'cited but not listed in literature: {sorted(cited - published)}'

    columns = set(default_run.columns('mutations'))
    detected = {name for row in marker_rows for name in row['Mutations in your sample'].split('; ') if name}
    assert detected <= columns, f'reported but not tabulated: {sorted(detected - columns)}'


def test_excel_report_mirrors_the_tsvs(default_run) -> None:
    workbook = load_workbook(default_run.paths['excel'], keep_vba=True)
    assert {'Mutations', 'Markers', 'Literature', 'Checks'} <= set(workbook.sheetnames)

    for sheet_name, output in (('Markers', 'markers'), ('Literature', 'literature')):
        rows = sum(1 for row in workbook[sheet_name].iter_rows(min_row=2) if row[0].value not in (None, ''))
        assert rows == len(default_run.rows(output)), f'{sheet_name} sheet does not match {output}.tsv'


def test_relaxed_reports_a_superset_of_strict_markers(default_run, relaxed_run) -> None:
    """``--relaxed`` loosens the detection threshold, so it can only add markers."""
    strict = {(row['Sample'], row['Marker']) for row in default_run.rows('markers')}
    relaxed = {(row['Sample'], row['Marker']) for row in relaxed_run.rows('markers')}
    assert strict <= relaxed, f'lost under --relaxed: {sorted(strict - relaxed)}'
