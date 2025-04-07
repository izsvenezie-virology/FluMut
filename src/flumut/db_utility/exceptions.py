class DBVersionError(Exception):
    def __init__(self, actual_version: str, expected_version: int, is_higher: bool) -> None:
        self.actual_version = actual_version
        self.expected_version = expected_version
        if is_higher:
            self.message = f'FluMutDB {self.actual_version} is too recent.\n'
            self.message += f'This version of FluMut works only with FluMutDB v.{expected_version}.x.\n'
            self.message += f'Please, search for FluMut updates.'
        else:
            self.message = f'FluMutDB {self.actual_version} is too old.\n'
            self.message += f'This version of FluMut works only with FluMutDB v.{expected_version}.x.\n'
            self.message += f'Please, update FluMutDB with `flumut --update`.'

        super().__init__(self.message)


class DBFileError(FileNotFoundError):
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.message = f'Unable to open FluMutDB file: {file_path}.'
        super().__init__(self.message)
