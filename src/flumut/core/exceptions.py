class FluMutError(Exception):
    """Exception base class for FluMut errors."""


class MissingFastaFilesError(FluMutError):
    """Raised when no FASTA files are provided to the workflow."""


class NoOutputRequestedError(FluMutError):
    """Raised when the workflow is asked to run without a single output file."""


class UnwritableOutputError(FluMutError):
    """Raised when an output file cannot be written where it was requested."""

    def __init__(self, path: str, reason: str) -> None:
        super().__init__(f'Cannot write output file {path}: {reason}')


class InvalidNameRegexError(FluMutError):
    """Raised when the regular expression used to parse sequence names does not compile."""

    def __init__(self, name_regex: str, reason: object) -> None:
        super().__init__(f'Invalid name regex "{name_regex}": {reason}')


class UpdateCheckError(FluMutError):
    """Raised when the latest release cannot be retrieved or understood."""
