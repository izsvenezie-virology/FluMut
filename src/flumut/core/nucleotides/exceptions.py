class CDSHasNoValidCodonsError(Exception):
    pass


class UnknownNucleotideError(Exception):
    def __init__(self, nucleotides: str, header: str | None = None) -> None:
        self.nucleotides = nucleotides
        message = f'Unknown nucleotide(s) "{nucleotides}"' + f' in sequence {header}.' if header else '.'
        super().__init__(message)
