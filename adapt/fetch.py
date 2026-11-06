# Author: Michel Sanner, PhD. Copyright (c) 2026. Portions adapted from ADCP_docking (BSD 3-Clause); see THIRD_PARTY_NOTICES.md.
"""Fetching an entry from the RCSB by its PDB id, and caching it.

Ported from ADCP_docking/adcp_dock/fetch.py by Aadhil Haq (BSD 3-Clause).
The PDB-format file is fetched; an entry the RCSB distributes only as mmCIF
(too large for the PDB format, or with chain ids of several characters) is
fetched as mmCIF, which ADFR 0.0.23 reads.

Downloads are cached under <outdir>/_rcsb/, so trying the same entry again
does not touch the network.
"""
import gzip
import os
import re
import urllib.error
import urllib.request

CODE = re.compile(r'^[0-9][A-Za-z0-9]{3}$')
URL = 'https://files.rcsb.org/download/%s.%s'
TIMEOUT = 30


class FetchError(Exception):
    pass


def looks_like_code(text):
    return bool(CODE.match((text or '').strip()))


def fetch(code, cachedir, timeout=TIMEOUT):
    """(path, url, was_cached) for entry `code`: its PDB-format file, else its
    mmCIF file."""
    code = (code or '').strip().upper()
    if not looks_like_code(code):
        raise FetchError('%r is not a PDB id. They are four characters, the '
                         'first a digit, such as 6GGD, 1A3R or 4A1X.' % code)
    for ext in ('pdb', 'cif'):
        out = os.path.join(cachedir, '%s.%s' % (code, ext))
        if os.path.isfile(out) and os.path.getsize(out) > 200:
            return out, URL % (code, ext), True
    os.makedirs(cachedir, exist_ok=True)
    for ext in ('pdb', 'cif'):
        url = URL % (code, ext)
        data = _download(code, url, timeout)
        if data is not None:
            break
    else:
        raise FetchError('%s is not at the RCSB (no PDB-format or mmCIF file): '
                         'check the id.' % code)
    if data[:2] == b'\x1f\x8b':
        data = gzip.decompress(data)
    if len(data) < 200:
        raise FetchError('the RCSB returned an empty file for %s' % code)
    out = os.path.join(cachedir, '%s.%s' % (code, ext))
    tmp = out + '.part'
    with open(tmp, 'wb') as fh:
        fh.write(data)
    os.replace(tmp, out)
    return out, url, False


def _download(code, url, timeout):
    """The bytes at `url`, or None if the RCSB has no such file (404)."""
    from . import APP, __version__
    req = urllib.request.Request(url, headers={'User-Agent': '%s/%s' % (APP, __version__)})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise FetchError('could not fetch %s from the RCSB: HTTP %s' % (code, exc.code))
    except Exception as exc:                                  # noqa: BLE001
        raise FetchError('could not fetch %s from the RCSB (%s). If this machine '
                         'has no route to the internet, download the file '
                         'elsewhere and open it with Open PDB.' % (code, exc))
