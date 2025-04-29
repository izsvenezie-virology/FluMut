class UnmatchingHeaderException(Exception):
    def __init__(self, header, regex) -> None:
        self.header = header
        self.regex = regex
        self.message = f'Unable to parse "{header}" with regular expression "{regex}".'
        super().__init__(self.message)


class PermissionDeniedException(Exception):
    def __init__(self, file_name) -> None:
        self.file_name = file_name
        self.message = f'Permission denied while trying to write "{file_name}".'
        super().__init__(self.message)
