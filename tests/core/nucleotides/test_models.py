import pytest

from flumut.core.globals import GAP_SYMBOL
from flumut.core.nucleotides.models import Alignment

G = GAP_SYMBOL


@pytest.mark.parametrize(
    'reference, expected',
    [
        ([], []),
        (['A'], [1]),
        ([G], [None]),
        (['A', 'C', 'T'], [1, 2, 3]),
        ([G, 'A', 'C'], [None, 1, 2]),
        (['A', 'C', G], [1, 2, None]),
        (['A', G, 'C'], [1, None, 2]),
        (['A', G, G, 'C'], [1, None, None, 2]),
        ([G, G, G], [None, None, None]),
    ],
    ids=[
        'empty',
        'single_residue',
        'single_gap',
        'no_gaps',
        'gap_at_start',
        'gap_at_end',
        'gap_in_middle',
        'consecutive_gaps',
        'all_gaps',
    ],
)
def test_get_positions(reference: list[str], expected: list) -> None:
    """Reference columns are numbered from 1; gaps map to None."""
    assert Alignment(reference=reference).get_positions() == expected
