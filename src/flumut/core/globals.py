# Sequence symbols
from importlib.resources import files

from peewee import DatabaseProxy

import flumut

WILDCARD = 'N'  # Wildcard character for not known nucleotides
GAP_SYMBOL = '-'  # Symbol for gaps in the alignment
STOP_CODON_SYMBOL = '*'  # Symbol for stop codons in the alignment
UNKNOWN_AA_SYMBOL = '?'  # Symbol for unknown amino acids in the alignment

# Alignment scores
MISMATCH_SCORE = -1  # Score for a mismatch between the reference and the sample sequence
GAP_OPEN_SCORE = -5  # Score for opening a gap in the alignment
GAP_EXTEND_SCORE = -1  # Score for extending a gap in the alignment
GAP_END_SCORE = -0.1  # Score for gaps at the end of the alignment, must be low to allow for partial sequences

# Paths
DB_FILE: str = str(files(flumut) / 'data' / 'flumut_db.sqlite')
EXCEL_TEMPLATE: str = str(files(flumut) / 'data' / 'flumut_output.xlsm')

# DB major version required
DB_MAJOR_VERSION = 7

# DB proxy
DATABASE_PROXY = DatabaseProxy()
