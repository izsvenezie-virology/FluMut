"""Check whether a newer FluMut release is published on GitHub.

Releases are tagged ``v.MAJOR.MINOR.PATCH``, plus a suffix for a pre-release.
The check is a convenience, never a requirement: :func:`check_for_update`
reports any failure as "no update known", so it cannot cost more than a log
message.
"""

import json
from dataclasses import dataclass
from urllib.request import Request, urlopen

from flumut import __version__
from flumut.core.logger import LOGGER

LATEST_RELEASE_API = 'https://api.github.com/repos/izsvenezie-virology/FluMut/releases/latest'
LATEST_RELEASE_PAGE = 'https://github.com/izsvenezie-virology/FluMut/releases/latest'
DEFAULT_TIMEOUT = 5.0


class UpdateCheckError(RuntimeError):
    """Raised when the latest release cannot be retrieved or understood."""


@dataclass(frozen=True)
class Release:
    """A published FluMut release: its version, and the page presenting it."""

    version: str
    url: str


def check_for_update(current_version: str = __version__, timeout: float = DEFAULT_TIMEOUT) -> Release | None:
    """Return the latest release if it is newer than ``current_version``, else None."""
    try:
        release = fetch_latest_release(timeout)
    except UpdateCheckError as e:
        LOGGER.debug(f'Update check failed: {e}')
        return None
    return release if is_newer(release.version, current_version) else None


def fetch_latest_release(timeout: float = DEFAULT_TIMEOUT) -> Release:
    """Ask GitHub for the latest published release, newer or not.

    Raises:
        UpdateCheckError: If GitHub cannot be reached, or does not answer with a tagged release.
    """
    request = Request(LATEST_RELEASE_API)
    try:
        with urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode('utf-8'))
        return Release(version=payload['tag_name'], url=payload.get('html_url', LATEST_RELEASE_PAGE))
    except OSError as e:
        raise UpdateCheckError(f'Cannot reach {LATEST_RELEASE_API}: {e}') from e
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        raise UpdateCheckError(f'Unreadable answer from {LATEST_RELEASE_API}: {e}') from e
    except KeyError as e:
        raise UpdateCheckError(f'No tag_name in answer from {LATEST_RELEASE_API}: {e}') from e


def is_newer(candidate: str, current: str) -> bool:
    """Whether ``candidate`` is a later version than ``current``, False if either cannot be read.

    A suffix marks a pre-release, ranking below the version it leads to: 1.0.0rc1 before 1.0.0.
    """
    candidate_keys = list(map(int, candidate.replace('v.', '').split('.')))
    current_keys = list(map(int, current.replace('v.', '').split('.')))
    return candidate_keys > current_keys
