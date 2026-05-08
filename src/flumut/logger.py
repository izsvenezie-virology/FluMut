import logging
import sys

LEVELS = {
    'dbg': logging.DEBUG,
    'inf': logging.INFO,
    'wrn': logging.WARNING,
    'err': logging.ERROR,
}


def initialize_logging(logging_level: str) -> None:
    logging.basicConfig(
        level=LEVELS[logging_level],
        format='[%(levelname)s] %(message)s',
        stream=sys.stderr,
        force=True,
    )
