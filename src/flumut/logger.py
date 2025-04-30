import logging
import sys


def initialize_logging(logging_level: str) -> None:
    '''
    Initialize the logger.

    :param `str` logging_level: The verbosity level of the logger.
    '''
    logging.basicConfig(
        level=_logging_levels[logging_level],
        format='[%(levelname)s] %(message)s',
        stream=sys.stderr,
        force=True
    )


_logging_levels = {
    'dbg': logging.DEBUG,
    'inf': logging.INFO,
    'wrn': logging.WARN,
    'err': logging.ERROR
}
