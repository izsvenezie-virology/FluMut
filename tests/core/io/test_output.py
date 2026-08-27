from io import StringIO
from unittest.mock import MagicMock

from flumut.core.io.output import (
    get_literature_data,
    get_markers_data,
    get_mutations_data,
    write_tsv,
)

# ---------------------------------------------------------------------------
# write_tsv tests
# ---------------------------------------------------------------------------


def test_write_tsv_single_row() -> None:
    output = StringIO()
    write_tsv(output, [{'Name': 'A', 'Value': '1'}])  # type: ignore[arg-type]
    assert output.getvalue() == 'Name\tValue\nA\t1\n'


def test_write_tsv_multiple_rows_written_in_order() -> None:
    output = StringIO()
    write_tsv(output, [{'X': 'a'}, {'X': 'b'}, {'X': 'c'}])  # type: ignore[arg-type]
    lines = output.getvalue().splitlines()
    assert lines == ['X', 'a', 'b', 'c']


# ---------------------------------------------------------------------------
# get_literature_data tests
# ---------------------------------------------------------------------------


def test_get_literature_data_empty_literature_returns_empty_list() -> None:
    analysis = MagicMock()
    analysis.literature = set()
    assert get_literature_data(analysis) == []


def test_get_literature_data_single_paper_all_fields_present() -> None:
    paper = MagicMock()
    paper.short_name = 'Doe2020'
    paper.title = 'Some Title'
    paper.authors = 'Doe et al.'
    paper.year = 2020
    paper.journal = 'Nature'
    paper.url = 'https://example.com'
    paper.doi = '10.1234/example'
    analysis = MagicMock()
    analysis.literature = {paper}
    result = get_literature_data(analysis)
    assert len(result) == 1
    assert result[0] == {
        'Short name': 'Doe2020',
        'Title': 'Some Title',
        'Authors': 'Doe et al.',
        'Year': 2020,
        'Journal': 'Nature',
        'Link': 'https://example.com',
        'DOI': '10.1234/example',
    }


# ---------------------------------------------------------------------------
# get_markers_data tests
# ---------------------------------------------------------------------------


def _make_evidence(effect: str, subtype: str, paper: str, host: str | None = None) -> MagicMock:
    ev = MagicMock()
    ev.effect.name = effect
    ev.subtype.name = subtype
    ev.paper.short_name = paper
    if host:
        ev.host = MagicMock()
        ev.host.name = host
    else:
        ev.host = None
    return ev


def _make_scan(marker_name: str, mutation_names: list[str], evidences: list) -> MagicMock:
    scan = MagicMock()
    scan.marker.name = marker_name
    scan.marker.evidences = evidences
    detected = []
    for n in mutation_names:
        dm = MagicMock()
        dm.mutation.name = n
        detected.append(dm)
    scan.detected_mutations = detected
    return scan


def test_get_markers_data_no_samples_returns_empty_list() -> None:
    analysis = MagicMock()
    analysis.samples = {}
    assert get_markers_data(analysis) == []


def test_get_markers_data_basic_row() -> None:
    evidence = _make_evidence('Increased replication', 'H5N1', 'Doe2020')
    scan = _make_scan(marker_name='I97T', mutation_names=['I97T'], evidences=[evidence])
    sample = MagicMock()
    sample.id = 'sample1'
    sample.marker_scans = [scan]
    analysis = MagicMock()
    analysis.samples = {'sample1': sample}
    result = get_markers_data(analysis)
    assert len(result) == 1
    assert result[0]['Sample'] == 'sample1'
    assert result[0]['Marker'] == 'I97T'
    assert result[0]['Effect'] == 'Increased replication'
    assert result[0]['Subtype'] == 'H5N1'
    assert result[0]['Literature'] == 'Doe2020'
    assert result[0]['Mutations in your sample'] == 'I97T'


def test_get_markers_data_evidence_with_host_appends_host_to_effect() -> None:
    evidence = _make_evidence('Increased replication', 'H5N1', 'Doe2020', host='Chicken')
    scan = _make_scan(marker_name='I97T', mutation_names=['I97T'], evidences=[evidence])
    sample = MagicMock()
    sample.id = 'sample1'
    sample.marker_scans = [scan]
    analysis = MagicMock()
    analysis.samples = {'sample1': sample}
    result = get_markers_data(analysis)
    assert result[0]['Effect'] == 'Increased replication in Chicken'


def test_get_markers_data_multiple_papers_same_effect_joined_with_semicolon() -> None:
    ev1 = _make_evidence('Increased replication', 'H5N1', 'Doe2020')
    ev2 = _make_evidence('Increased replication', 'H5N1', 'Smith2021')
    scan = _make_scan(marker_name='I97T', mutation_names=['I97T'], evidences=[ev1, ev2])
    sample = MagicMock()
    sample.id = 'sample1'
    sample.marker_scans = [scan]
    analysis = MagicMock()
    analysis.samples = {'sample1': sample}
    result = get_markers_data(analysis)
    assert len(result) == 1
    assert result[0]['Literature'] == 'Doe2020; Smith2021'


# ---------------------------------------------------------------------------
# get_mutations_data tests
# ---------------------------------------------------------------------------


def test_get_mutations_data_no_samples_returns_empty_list() -> None:
    analysis = MagicMock()
    analysis.mutations = set()
    analysis.samples = {}
    assert get_mutations_data(analysis) == []


def test_get_mutations_data_detected_mutation_stored_in_row() -> None:
    mutation = MagicMock()
    mutation.name = 'I97T'
    mutation.default_position = 97
    pos = MagicMock()
    pos.mutation = mutation
    pos.ammino_acid = 'T'
    sample = MagicMock()
    sample.id = 'sample1'
    sample.positions = [pos]
    analysis = MagicMock()
    analysis.mutations = {mutation}
    analysis.samples = {'sample1': sample}
    result = get_mutations_data(analysis)
    assert len(result) == 1
    assert result[0]['Sample'] == 'sample1'
    assert result[0]['I97T'] == 'T'


def test_get_mutations_data_undetected_mutation_absent_from_row() -> None:
    mutation = MagicMock()
    mutation.name = 'D701N'
    mutation.default_position = 701
    sample = MagicMock()
    sample.id = 'sample1'
    sample.positions = []  # no positions detected
    analysis = MagicMock()
    analysis.mutations = {mutation}
    analysis.samples = {'sample1': sample}
    result = get_mutations_data(analysis)
    assert 'D701N' not in result[0]


def test_get_mutations_data_each_sample_gets_a_row() -> None:
    analysis = MagicMock()
    analysis.mutations = set()
    s1, s2 = MagicMock(positions=[]), MagicMock(positions=[])
    s1.id, s2.id = 'sample1', 'sample2'
    analysis.samples = {'sample1': s1, 'sample2': s2}
    result = get_mutations_data(analysis)
    assert len(result) == 2
    assert {r['Sample'] for r in result} == {'sample1', 'sample2'}
