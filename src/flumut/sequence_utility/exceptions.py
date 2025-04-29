class UnknownNucleotideException(Exception):
    def __init__(self, codon) -> None:
        self.codon = codon
        self.message = f'Unexpected nucleotide in codon "{codon}".'
        super().__init__(self.message)


class MalformedFastaException(Exception):
    def __init__(self, fasta) -> None:
        self.fasta = fasta
        self.message = f'{fasta} does not start whit ">".'
        super().__init__(self.message)
