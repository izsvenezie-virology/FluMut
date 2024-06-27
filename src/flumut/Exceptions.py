
class UnmatchNameException(Exception):
    def __init__(self, name, regex) -> None:
        self.name = name
        self.regex = regex
        self.message = f'Unable to parse "{name}" with regular expression "{regex}".'
        super().__init__(self.message)


class UnknownSegmentExeption(Exception):
    def __init__(self, name, regex, segment) -> None:
        self.name = name
        self.regex = regex
        self.message = f'Unrecognized segment "{segment}", found in "{name}" parsed with "{regex}".'
        super().__init__(self.message)

class UnknownNucleotideExeption(Exception):
    def __init__(self, codon) -> None:
        self.codon = codon
        self.message = f'Unexpected nucleotide in codon "{codon}".'
        super().__init__(self.message)
