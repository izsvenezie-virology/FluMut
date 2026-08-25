"""Bibliographic helpers for `Paper` records: Crossref lookup and short name generation."""

import html
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from string import ascii_lowercase

from flumut.flumutdb.models import Paper

CROSSREF_API = 'https://api.crossref.org/works/'
USER_AGENT = 'FluMut-db-editor (https://github.com/izsvenezie-virology/FluMut)'
REQUEST_TIMEOUT = 15

# Crossref serves dates under several keys, from the most to the least specific.
DATE_KEYS = ('issued', 'published', 'published-print', 'published-online', 'created')

DOI_PREFIX = re.compile(r'^(?:https?://)?(?:dx\.)?(?:doi\.org/|doi:)\s*', re.IGNORECASE)
MARKUP = re.compile(r'<[^>]+>')
WHITESPACE = re.compile(r'\s+')


class CrossrefError(Exception):
    """Raised when the metadata of a DOI cannot be retrieved."""


@dataclass
class PaperMetadata:
    """Crossref metadata mapped onto the fields of `Paper`."""

    doi: str
    title: str
    authors: str
    year: int | None
    journal: str
    url: str


def normalize_doi(raw: str) -> str:
    """Strip the `doi:` and `doi.org/` prefixes users paste along with a DOI."""
    return DOI_PREFIX.sub('', (raw or '').strip()).strip()


def fetch_metadata(doi: str) -> PaperMetadata:
    """Retrieve the metadata of `doi` from Crossref.

    Raises `CrossrefError` when the DOI is unknown or Crossref cannot be reached.
    """
    doi = normalize_doi(doi)
    if not doi:
        raise CrossrefError('No DOI to look up.')

    request = urllib.request.Request(CROSSREF_API + urllib.parse.quote(doi), headers={'User-Agent': USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise CrossrefError(f'DOI "{doi}" is not registered on Crossref.') from error
        raise CrossrefError(f'Crossref answered with HTTP {error.code} for DOI "{doi}".') from error
    except urllib.error.URLError as error:
        raise CrossrefError(f'Could not reach Crossref: {error.reason}') from error
    except (TimeoutError, json.JSONDecodeError) as error:
        raise CrossrefError(f'Crossref did not answer properly: {error}') from error

    message = payload.get('message') or {}
    return PaperMetadata(
        doi=doi,
        title=_clean(_first(message.get('title'))),
        authors=_format_authors(message.get('author') or []),
        year=_extract_year(message),
        journal=_clean(_first(message.get('container-title'))),
        url=message.get('URL') or f'https://doi.org/{doi}',
    )


def build_short_name(authors: str, year: int | None) -> str:
    """Build the citation label used across the database, e.g. `Chutinimitkul S. et al., 2010`.

    `authors` follows the stored format: `Family, Given` entries joined by `;`.
    Returns an empty string when there is not enough data to build a label.
    """
    entries = [entry.strip() for entry in (authors or '').split(';') if entry.strip()]
    if not entries or not year:
        return ''

    family, _, given = entries[0].partition(',')
    family = family.strip()
    given = given.strip()
    if not family:
        return ''

    # Consortia and institutions are stored as a single name without a given part.
    label = f'{family} {given[0]}.' if given else family
    if len(entries) > 1:
        label += ' et al.'
    return f'{label}, {year}'


def disambiguate_short_name(short_name: str, instance: Paper | None = None) -> str:
    """Append the `b`, `c`, ... year suffix the database uses when a short name is already taken.

    `instance` is the paper being edited, so that it does not collide with itself.
    """
    if not short_name:
        return short_name
    for suffix in ('', *ascii_lowercase[1:]):
        candidate = short_name + suffix
        taken = Paper.get_or_none(Paper.short_name == candidate)
        if taken is None or taken == instance:
            return candidate
    return short_name


def _first(values: list | None) -> str:
    return values[0] if values else ''


def _clean(text: str) -> str:
    """Drop the JATS/HTML markup and line breaks Crossref keeps inside titles."""
    return WHITESPACE.sub(' ', html.unescape(MARKUP.sub('', text or ''))).strip()


def _format_authors(entries: list[dict]) -> str:
    names = []
    for entry in entries:
        family = _clean(entry.get('family', ''))
        given = _clean(entry.get('given', ''))
        if family:
            names.append(f'{family}, {given}' if given else family)
        elif name := _clean(entry.get('name', '')):
            names.append(name)
    return ';'.join(names)


def _extract_year(message: dict) -> int | None:
    for key in DATE_KEYS:
        parts = (message.get(key) or {}).get('date-parts') or []
        if parts and parts[0] and parts[0][0]:
            return int(parts[0][0])
    return None
