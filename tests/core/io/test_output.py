from io import StringIO
from unittest.mock import MagicMock

import pytest

from flumut.core.io.output import get_literature_data, get_markers_data, get_mutations_data, write_tsv

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_analysis(samples: dict | None = None, literature: tuple = (), mutations: tuple = ()) -> MagicMock:
    analysis = MagicMock()
    analysis.samples = samples or {}
    analysis.literature = set(literature)
    analysis.mutations = set(mutations)
    return analysis


def _make_sample(sample_id: str, marker_scans: tuple = (), positions: tuple = ()) -> MagicMock:
    sample = MagicMock()
    sample.id = sample_id
    sample.marker_scans = list(marker_scans)
    sample.positions = list(positions)
    return sample


def _make_paper(short_name: str, **fields) -> MagicMock:
    paper = MagicMock()
    paper.short_name = short_name
    for field, value in fields.items():
        setattr(paper, field, value)
    return paper


def _make_evidence(effect: str, subtype: str, paper: str, host: str | None = None) -> MagicMock:
    evidence = MagicMock()
    evidence.effect.name = effect
    evidence.subtype.name = subtype
    evidence.paper.short_name = paper
    evidence.host = None
    if host is not None:
        evidence.host = MagicMock()
        evidence.host.name = host
    return evidence


def _make_scan(marker_name: str, mutation_names: tuple, evidences: list) -> MagicMock:
    scan = MagicMock()
    scan.marker.name = marker_name
    scan.marker.evidences = evidences
    scan.detected_mutations = [MagicMock(**{'mutation.name': name}) for name in mutation_names]
    return scan


# ---------------------------------------------------------------------------
# write_tsv
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'rows, expected',
    [
        ([{'Name': 'A', 'Value': '1'}], ['Name\tValue', 'A\t1']),
        ([{'X': 'a'}, {'X': 'b'}, {'X': 'c'}], ['X', 'a', 'b', 'c']),
    ],
    ids=['header_from_first_row', 'rows_written_in_order'],
)
def test_write_tsv(rows: list[dict], expected: list[str]) -> None:
    output = StringIO()
    write_tsv(output, rows)  # type: ignore[arg-type]
    assert output.getvalue().splitlines() == expected


# ---------------------------------------------------------------------------
# Extractors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'extractor',
    [get_literature_data, get_markers_data, get_mutations_data],
    ids=['literature', 'markers', 'mutations'],
)
def test_extractors_return_nothing_for_an_empty_analysis(extractor) -> None:
    assert extractor(_make_analysis()) == []


def test_get_literature_data_returns_every_column() -> None:
    paper = _make_paper(
        'Doe2020',
        title='Some Title',
        authors='Doe et al.',
        year=2020,
        journal='Nature',
        url='https://example.com',
        doi='10.1234/example',
    )
    assert get_literature_data(_make_analysis(literature=(paper,))) == [
        {
            'Short name': 'Doe2020',
            'Title': 'Some Title',
            'Authors': 'Doe et al.',
            'Year': 2020,
            'Journal': 'Nature',
            'Link': 'https://example.com',
            'DOI': '10.1234/example',
        }
    ]


def test_get_literature_data_is_sorted_by_short_name() -> None:
    """``analysis.literature`` is a set, so the output must impose its own order to stay reproducible."""
    papers = (_make_paper('Zhang2023'), _make_paper('Doe2020'), _make_paper('Smith2021'))
    rows = get_literature_data(_make_analysis(literature=papers))
    assert [row['Short name'] for row in rows] == ['Doe2020', 'Smith2021', 'Zhang2023']


@pytest.mark.parametrize(
    'evidences, expected_effect, expected_literature',
    [
        ((('Increased replication', 'H5N1', 'Doe2020', None),), 'Increased replication', 'Doe2020'),
        ((('Increased replication', 'H5N1', 'Doe2020', 'Chicken'),), 'Increased replication in Chicken', 'Doe2020'),
        (
            (('Increased replication', 'H5N1', 'Doe2020', None), ('Increased replication', 'H5N1', 'Smith2021', None)),
            'Increased replication',
            'Doe2020; Smith2021',
        ),
    ],
    ids=['single_evidence', 'host_appended_to_effect', 'papers_of_one_effect_joined'],
)
def test_get_markers_data(evidences: tuple, expected_effect: str, expected_literature: str) -> None:
    """One row per sample x marker x (effect, subtype), with that group's papers joined."""
    scan = _make_scan('I97T', ('I97T',), [_make_evidence(*evidence) for evidence in evidences])
    analysis = _make_analysis(samples={'sample1': _make_sample('sample1', marker_scans=(scan,))})

    assert get_markers_data(analysis) == [
        {
            'Sample': 'sample1',
            'Marker': 'I97T',
            'Mutations in your sample': 'I97T',
            'Effect': expected_effect,
            'Subtype': 'H5N1',
            'Literature': expected_literature,
        }
    ]


@pytest.mark.parametrize(
    'scanned, expected_row',
    [
        (True, {'Sample': 'sample1', 'I97T': 'T'}),
        (False, {'Sample': 'sample1'}),
    ],
    ids=['scanned_in_sample', 'not_scanned_in_sample'],
)
def test_get_mutations_data(scanned: bool, expected_row: dict) -> None:
    """A mutation only gets a column for a sample where it was actually scanned."""
    mutation = MagicMock()
    mutation.name = 'I97T'
    positions = (MagicMock(mutation=mutation, ammino_acid='T'),) if scanned else ()
    analysis = _make_analysis(
        samples={'sample1': _make_sample('sample1', positions=positions)},
        mutations=(mutation,),
    )
    assert get_mutations_data(analysis) == [expected_row]


def test_get_mutations_data_writes_one_row_per_sample() -> None:
    samples = {name: _make_sample(name) for name in ('sample1', 'sample2')}
    rows = get_mutations_data(_make_analysis(samples=samples))
    assert [row['Sample'] for row in rows] == ['sample1', 'sample2']
