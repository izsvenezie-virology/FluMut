import logging
import re

from flumut.exceptions import UnmatchingHeaderException


def parse_header(header: str, pattern: re.Pattern, allow_unmatching_headers: bool = False) -> tuple[str, str | None]:
    """
    Parse the header with the header pattern.
    It searches for two groups in order to retrieve sample and segment from the header.
    The header pattern can be set with set_header_pattern function.

    :param `str` header: The FASTA header.
    :return `Optional[str]` sample: The sample name. Returns `None` if  not found.
    :return `Optional[str]`segment: The segment name. Returns `None` if not found.
    """
    logging.debug(f'Parsing "{header}" with "{pattern.pattern}"')
    sample, segment = None, None
    match = pattern.match(header)

    if not match:
        if not allow_unmatching_headers:
            raise UnmatchingHeaderException(header, pattern.pattern) from None
        logging.info(f'Cannot extract sample ID from "{header}". The whole header is used as sample ID.')
        return header, None

    try:
        sample = match.groupdict().get('sample', match.group(1))
        segment = match.groupdict().get('segment', match.group(2))
    except IndexError:
        pass

    logging.debug(f'Sample:  "{sample}" - Segment: "{segment}"')
    return sample or header, segment
