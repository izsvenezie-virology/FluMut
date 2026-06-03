import logging
import sys

LOGGER = logging.getLogger(__name__)

LEVELS = {
    'dbg': logging.DEBUG,
    'inf': logging.INFO,
    'wrn': logging.WARNING,
    'err': logging.ERROR,
}


def initialize_logging() -> None:
    logging.basicConfig(
        level=LEVELS['wrn'],
        format='[%(levelname)s] %(message)s',
        stream=sys.stderr,
        force=True,
    )


def set_level(level: str) -> None:
    LOGGER.setLevel(LEVELS[level])
