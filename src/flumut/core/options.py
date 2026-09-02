from collections.abc import Mapping
from dataclasses import dataclass, fields
from io import TextIOWrapper
from typing import Any, TypeVar

from typing_extensions import Self

# Default options
DEFAULT_NAME_REGEX = r'(?P<sample>.+)_(?P<segment>.+)'
DEFAULT_RELAXED = False

OptionGroupT = TypeVar('OptionGroupT', bound='OptionGroup')


class OptionGroup:
    """Base of the option groups, giving each one a constructor from a flat mapping."""

    @classmethod
    def from_flat(cls, values: Mapping[str, Any]) -> Self:
        """Build the group from ``values``, taking only the keys it declares.

        Keys the group does not declare are left to the other groups, and keys
        absent from ``values`` keep their default, so a frontend has to provide
        only the options it actually exposes.
        """
        declared = {field.name for field in fields(cls)}  # type: ignore[arg-type]
        return cls(**{name: value for name, value in values.items() if name in declared})


@dataclass(frozen=True)
class InputOptions(OptionGroup):
    """Where the sequences come from and how their names are read."""

    fasta_files: tuple[TextIOWrapper, ...] = ()
    name_regex: str = DEFAULT_NAME_REGEX


@dataclass(frozen=True)
class AnalysisOptions(OptionGroup):
    """How markers are called."""

    relaxed: bool = DEFAULT_RELAXED


@dataclass(frozen=True)
class OutputOptions(OptionGroup):
    """The reports to write, each an open file handle or None when not requested."""

    markers_output: TextIOWrapper | None = None
    mutations_output: TextIOWrapper | None = None
    literature_output: TextIOWrapper | None = None
    excel_output: TextIOWrapper | None = None

    def any_requested(self) -> bool:
        """Whether at least one report was asked for."""
        return any(getattr(self, field.name) for field in fields(self))


@dataclass(frozen=True)
class DatabaseOptions(OptionGroup):
    """Which FluMutDB to open, and how. ``path`` of None means the bundled database."""

    path: str | None = None
    read_only: bool = True


#: Attribute name of each group on :class:`FluMutOptions`, mapped to its type.
_GROUPS: dict[str, type[OptionGroup]] = {
    'input': InputOptions,
    'output': OutputOptions,
    'analysis': AnalysisOptions,
    'database': DatabaseOptions,
}


@dataclass(frozen=True)
class FluMutOptions:
    """Everything one run needs, grouped by concern."""

    input: InputOptions = InputOptions()
    output: OutputOptions = OutputOptions()
    analysis: AnalysisOptions = AnalysisOptions()
    database: DatabaseOptions = DatabaseOptions()

    @classmethod
    def from_flat(cls, values: Mapping[str, Any]) -> 'FluMutOptions':
        """Build every group from one flat mapping, such as Click's parsed parameters.

        Args:
            values: Option values keyed by field name. Missing keys keep their default.

        Raises:
            TypeError: If ``values`` holds a key no group declares, which would
                otherwise be dropped without a trace.
        """
        declared = {field.name for group in _GROUPS.values() for field in fields(group)}  # type: ignore[arg-type]
        unknown = set(values) - declared
        if unknown:
            raise TypeError(f'No option group declares: {", ".join(sorted(unknown))}')
        return cls(**{name: group.from_flat(values) for name, group in _GROUPS.items()})  # type: ignore
