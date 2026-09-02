"""Tests for the GitHub release check.

The check must be silent about its own failures: whatever GitHub answers, the
caller only ever learns that an update exists or that it does not, so these
cover both the comparison of versions and every way the network call can go
wrong.
"""

import json
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from flumut import __version__
from flumut.cli import cli
from flumut.core.updates import (
    LATEST_RELEASE_PAGE,
    Release,
    check_for_update,
    fetch_latest_release,
    is_newer,
)


def _mock_urlopen(payload: object = None, body: bytes | None = None, error: Exception | None = None) -> MagicMock:
    """Patch ``urlopen`` so it answers with ``payload``/``body``, or raises ``error``."""
    mock = MagicMock()
    if error is not None:
        mock.side_effect = error
    else:
        response = MagicMock()
        response.read.return_value = json.dumps(payload).encode('utf-8') if body is None else body
        mock.return_value.__enter__.return_value = response
    return mock


# ---------------------------------------------------------------------------
# is_newer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ('candidate', 'current'),
    [
        ('1.0.1', '1.0.0'),
        ('1.1.0', '1.0.9'),
        ('2.0.0', '1.9.9'),
        ('1.0.0', '0.6.5'),
        ('1.0.10', '1.0.9'),
    ],
    ids=['patch', 'minor', 'major', 'across-zero', 'two-digits'],
)
def test_is_newer(candidate: str, current: str) -> None:
    """A later version is recognised whichever component grew."""
    assert is_newer(candidate, current) is True
    assert is_newer(current, candidate) is False


def test_is_newer_is_false_for_the_same_version() -> None:
    """A tag prefix does not make a version newer than itself."""
    assert is_newer('v.1.0.0', 'v.1.0.0') is False


# ---------------------------------------------------------------------------
# fetch_latest_release
# ---------------------------------------------------------------------------


def test_fetch_latest_release_reads_tag_and_page() -> None:
    """The release is built from the tag and the page GitHub reports."""
    payload = {'tag_name': 'v.1.2.0', 'html_url': 'https://github.com/izsvenezie-virology/FluMut/releases/tag/v.1.2.0'}
    with patch('flumut.core.updates.urlopen', _mock_urlopen(payload)):
        release = fetch_latest_release()

    assert release == Release(version='v.1.2.0', url=payload['html_url'])


def test_fetch_latest_release_falls_back_to_the_releases_page() -> None:
    """A release without its own page still points somewhere useful."""
    with patch('flumut.core.updates.urlopen', _mock_urlopen({'tag_name': 'v.1.2.0'})):
        assert fetch_latest_release().url == LATEST_RELEASE_PAGE


def test_fetch_latest_release_identifies_itself_and_gives_up_in_time() -> None:
    """The request names FluMut to GitHub, and never waits without a limit."""
    urlopen = _mock_urlopen({'tag_name': '1.0.0'})
    with patch('flumut.core.updates.urlopen', urlopen):
        fetch_latest_release(timeout=1.5)

    request, kwargs = urlopen.call_args
    assert kwargs['timeout'] == 1.5
    assert request[0].get_header('User-agent').startswith('FluMut/')


# ---------------------------------------------------------------------------
# check_for_update
# ---------------------------------------------------------------------------


def test_check_for_update_reports_a_newer_release() -> None:
    """A later release is handed back to the caller."""
    with patch('flumut.core.updates.urlopen', _mock_urlopen({'tag_name': 'v.2.0.0'})):
        release = check_for_update('1.0.0')

    assert release is not None
    assert release.version == 'v.2.0.0'


@pytest.mark.parametrize('tag', ['v.1.0.0', 'v.0.6.5'], ids=['same', 'older'])
def test_check_for_update_is_silent_when_up_to_date(tag: str) -> None:
    """Nothing is reported unless the published release is actually newer."""
    with patch('flumut.core.updates.urlopen', _mock_urlopen({'tag_name': tag})):
        assert check_for_update('1.0.0') is None


@pytest.mark.parametrize(
    'kwargs',
    [{'error': OSError('unreachable')}, {'body': b'not json'}, {'payload': {'message': 'rate limited'}}],
    ids=['network', 'body', 'no-tag'],
)
def test_check_for_update_never_raises(kwargs: dict) -> None:
    """A failed check costs a log line, never the run that asked for it."""
    with patch('flumut.core.updates.urlopen', _mock_urlopen(**kwargs)):
        assert check_for_update('1.0.0') is None


def test_check_for_update_defaults_to_the_running_version() -> None:
    """Called with no argument, the check compares against the installed FluMut."""
    with patch('flumut.core.updates.urlopen', _mock_urlopen({'tag_name': __version__})):
        assert check_for_update() is None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ('tag', 'expected'),
    [('9.9.9', 'is available'), ('0.0.1', 'is up to date')],
    ids=['newer', 'up-to-date'],
)
def test_check_update_flag_reports_the_outcome(tag: str, expected: str) -> None:
    """``--check-update`` says what it found and exits without running an analysis."""
    with patch('flumut.core.updates.urlopen', _mock_urlopen({'tag_name': tag})):
        result = CliRunner().invoke(cli, ['--check-update'])

    assert result.exit_code == 0
    assert expected in result.output


def test_check_update_flag_reports_a_failed_check() -> None:
    """A check the user asked for fails loudly instead of claiming to be up to date."""
    with patch('flumut.core.updates.urlopen', _mock_urlopen(error=OSError('unreachable'))):
        result = CliRunner().invoke(cli, ['--check-update'])

    assert result.exit_code == 1
    assert 'up to date' not in result.output
