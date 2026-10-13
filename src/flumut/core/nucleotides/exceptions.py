from flumut.core.exceptions import FluMutError


class CDSHasNoValidCodonsError(FluMutError):
    pass


class UnknownNucleotideError(FluMutError):
    def __init__(self, nucleotides: str, header: str | None = None) -> None:
        self.nucleotides = nucleotides
        message = f'Unknown nucleotide(s) "{nucleotides}"' + (f' in sequence {header}.' if header else '.')
        super().__init__(message)


class NoReferenceFoundError(FluMutError):
    pass
