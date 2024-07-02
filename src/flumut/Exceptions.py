import sys

class UnmatchNameException(Exception):
    def __init__(self, name, regex) -> None:
        self.name = name
        self.regex = regex
        self.message = f'Unable to parse "{name}" with regular expression "{regex}".'
        super().__init__(self.message)


class UnknownSegmentException(Exception):
    def __init__(self, name, regex, segment) -> None:
        self.name = name
        self.regex = regex
        self.message = f'Unrecognized segment "{segment}", found in "{name}" parsed with "{regex}".'
        super().__init__(self.message)


class UnknownNucleotideException(Exception):
    def __init__(self, codon) -> None:
        self.codon = codon
        self.message = f'Unexpected nucleotide in codon "{codon}".'
        super().__init__(self.message)

class MalformedFastaException(Exception):
    def __init__(self) -> None:
        self.message = 'Provided FASTA file does not start whit ">".'
        super().__init__(self.message)

class PermissionDeniedException(Exception):
    def __init__(self, file_name) -> None:
        self.file_name = file_name
        self.message = f'Permission denied while trying to write "{file_name}".'
        super().__init__(self.message)
