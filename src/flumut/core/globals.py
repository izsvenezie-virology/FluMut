WILDCARD = 'N'  # Wildcard character for not known nucleotides
GAP_SYMBOL = '-'  # Symbol for gaps in the alignment
MISMATCH_SCORE = -1  # Score for a mismatch between the reference and the sample sequence
GAP_OPEN_SCORE = -5  # Score for opening a gap in the alignment
GAP_EXTEND_SCORE = -1  # Score for extending a gap in the alignment
GAP_END_SCORE = -0.1  # Score for gaps at the end of the alignment, must be low to allow for partial sequences
