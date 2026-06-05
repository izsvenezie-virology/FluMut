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
    """Configure root logging with WARNING level output to stderr."""
    logging.basicConfig(
        level=LEVELS['wrn'],
        format='[%(levelname)s] %(message)s',
        stream=sys.stderr,
        force=True,
    )


def set_level(level: str) -> None:
    """Set the log level for the flumut logger.

    Args:
        level: Short level key; one of ``'dbg'``, ``'inf'``, ``'wrn'``, ``'err'``.
    """
    LOGGER.setLevel(LEVELS[level])
