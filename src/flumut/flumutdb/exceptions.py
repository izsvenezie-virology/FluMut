class FluMutDBError(Exception):
    pass


class DatabaseNotFoundError(FluMutDBError, FileNotFoundError):
    """Raised when the database file does not exist."""


class InvalidDatabaseError(FluMutDBError):
    """Raised when the database file cannot be read as a FluMutDB database."""

    def __init__(self, path: str, reason: object) -> None:
        self.message = f'Cannot read database {path}: {reason}. Is it a FluMutDB file?'
        super().__init__(self.message)


class MissingVersionError(FluMutDBError):
    """Raised when the database version is missing."""

    def __init__(self, *args: object) -> None:
        self.message = 'Database version is missing. Is the database correctly formatted?'
        super().__init__(self.message, *args)


class IncompatibleVersionError(FluMutDBError):
    """Raised when the database version is incompatible with the application."""

    def __init__(self, version, required_version, *args: object) -> None:
        self.message = f'Incompatible database version: {version}. Expected major version: {required_version}.'
        super().__init__(self.message, *args)
