"""Tests for the option groups shared by every frontend.

Both the CLI and the GUI build a :class:`FluMutOptions` and hand it to the
workflow, so these cover the contract between them: which group owns which
option, what a missing option falls back to, and what happens to an option no
group declares.
"""

from dataclasses import fields

import pytest

from flumut.cli import cli
from flumut.core.options import (
    DEFAULT_NAME_REGEX,
    AnalysisOptions,
    DatabaseOptions,
    FluMutOptions,
    InputOptions,
    OutputOptions,
)

OPTION_GROUPS = (InputOptions, OutputOptions, AnalysisOptions, DatabaseOptions)

#: One value per option a frontend can set, flat, the way Click hands them over.
FLAT_VALUES = {
    'fasta_files': ('first.fa', 'second.fa'),
    'name_regex': r'(?P<sample>.+)-(?P<segment>.+)',
    'relaxed': True,
    'markers_output': 'markers.tsv',
    'mutations_output': 'mutations.tsv',
    'literature_output': 'literature.tsv',
    'excel_output': 'report.xlsm',
    'path': 'custom.sqlite',
}


# ---------------------------------------------------------------------------
# from_flat
# ---------------------------------------------------------------------------


def test_from_flat_routes_every_value_to_its_group() -> None:
    """A flat mapping is split across the groups that declare its keys."""
    options = FluMutOptions.from_flat(FLAT_VALUES)

    assert options.input == InputOptions(FLAT_VALUES['fasta_files'], FLAT_VALUES['name_regex'])
    assert options.output == OutputOptions('markers.tsv', 'mutations.tsv', 'literature.tsv', 'report.xlsm')
    assert options.analysis == AnalysisOptions(relaxed=True)
    assert options.database == DatabaseOptions(path='custom.sqlite')


def test_from_flat_falls_back_to_the_defaults() -> None:
    """A frontend provides only the options it exposes; the rest keep their default."""
    options = FluMutOptions.from_flat({'fasta_files': ('only.fa',)})

    assert options.input.name_regex == DEFAULT_NAME_REGEX
    assert options.analysis.relaxed is False
    assert options.database == DatabaseOptions(path=None, read_only=True)
    assert options.output.any_requested() is False


def test_from_flat_rejects_an_option_no_group_declares() -> None:
    """An option named after nothing would be dropped in silence, so it raises instead."""
    with pytest.raises(TypeError, match='relaxd'):
        FluMutOptions.from_flat({'relaxd': True})


# ---------------------------------------------------------------------------
# OutputOptions
# ---------------------------------------------------------------------------


def test_no_output_is_requested_by_default() -> None:
    assert OutputOptions().any_requested() is False


@pytest.mark.parametrize('output', [field.name for field in fields(OutputOptions)])
def test_any_single_output_counts_as_requested(output: str) -> None:
    """Every field is an output, so setting any one of them must be enough."""
    assert OutputOptions(**{output: 'somewhere.out'}).any_requested() is True


# ---------------------------------------------------------------------------
# CLI contract
# ---------------------------------------------------------------------------


def test_every_cli_parameter_belongs_to_an_option_group() -> None:
    """Adding a CLI option means adding the matching field, or nothing carries its value.

    ``all_versions`` is the one exception: the command takes it as a named
    argument and answers it before the workflow starts, so it never reaches a group.
    """
    exposed = {param.name for param in cli.params if param.expose_value} - {'all_versions'}
    declared = {field.name for group in OPTION_GROUPS for field in fields(group)}

    assert exposed <= declared, f'no option group declares: {sorted(exposed - declared)}'
