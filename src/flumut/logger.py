import logging
import sys
from typing import Final

TRACE: Final[int] = logging.DEBUG + 5
logging.addLevelName(TRACE, 'TRACE')

_logging_levels = {
    'dbg': logging.DEBUG,
    'trc': TRACE,
    'inf': logging.INFO,
    'wrn': logging.WARNING,
    'err': logging.ERROR,
}


def initialize_logging(logging_level: str) -> None:
    logging.basicConfig(
        level=_logging_levels[logging_level],
        format='[%(levelname)s] %(message)s',
        stream=sys.stderr,
        force=True,
    )


class Logger(logging.Logger):
    def trace(self, msg: object, *args: object) -> None:
        if self.isEnabledFor(TRACE):
            self._log(TRACE, msg, args)


logging.setLoggerClass(Logger)
