class MissingVersionError(Exception):
    """Raised when the database version is missing."""

    def __init__(self, *args: object) -> None:
        self.message = (
            "Database version is missing. Is the database correctly formatted?"
        )
        super().__init__(self.message, *args)


class IncompatibleVersionError(Exception):
    """Raised when the database version is incompatible with the application."""

    def __init__(self, version, required_version, *args: object) -> None:
        self.message = (
            f"Incompatible database version: {version}. "
            f"Expected major version: {required_version}."
        )
        super().__init__(self.message, *args)
