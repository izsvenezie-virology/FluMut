from flumutdb.models import MutationType

from flumut.scan.models import PositionScan
from flumut.translation.models import ProteinAlignment


def scan(alignment: ProteinAlignment) -> list[PositionScan]:
    positions = alignment.get_positions()
    scans = []

    for mutation in alignment.protein.mutations:
        for mapping in mutation.mappings:
            if not mapping.reference == alignment.reference:
                continue
            index = positions.index(mapping.position)
            aa = alignment.aligned_query[index]
            scans.append(PositionScan(mapping, aa))
    return scans


def is_positive(scan: PositionScan) -> bool:
    match scan.mapping.mutation.type:
        case MutationType.SNP.value:
            return str(scan.mapping.alteration) in scan.ammino_acid
        case _:
            raise NotImplementedError(f'Mutation type {scan.mapping.mutation.type} not supported')
