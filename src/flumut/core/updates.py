"""Check whether a newer FluMut release is published on GitHub.

Releases are tagged ``v.MAJOR.MINOR.PATCH``; there are no pre-release tags.
The check is a convenience, never a requirement: :func:`check_for_update`
reports any failure as "no update known", so it cannot cost more than a log
message.
"""

import json
import re
from dataclasses import dataclass
from urllib.request import Request, urlopen

from flumut import __version__
from flumut.core.logger import LOGGER

LATEST_RELEASE_API = 'https://api.github.com/repos/izsvenezie-virology/FluMut/releases/latest'
LATEST_RELEASE_PAGE = 'https://github.com/izsvenezie-virology/FluMut/releases/latest'
DEFAULT_TIMEOUT = 5.0

VERSION_PATTERN = re.compile(r'^v?\.?(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+)$')


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


def parse_version(version: str) -> tuple[int, int, int] | None:
    """Return the release numbers of ``version`` as a sortable key, or None if it cannot be read."""
    match = VERSION_PATTERN.match((version or '').strip())
    if match is None:
        return None
    return int(match['major']), int(match['minor']), int(match['patch'])


def is_newer(candidate: str, current: str) -> bool:
    """Whether ``candidate`` is a later version than ``current``, False if either cannot be read."""
    candidate_key = parse_version(candidate)
    current_key = parse_version(current)
    if candidate_key is None or current_key is None:
        LOGGER.debug(f'Cannot compare version "{candidate}" with "{current}"')
        return False
    return candidate_key > current_key
