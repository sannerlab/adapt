"""Fetching an entry from the RCSB by its PDB id, and caching it.

Ported from ADCP_docking/adcp_dock/fetch.py by Aadhil Haq (BSD 3-Clause).
One difference: only the PDB format is fetched. ADapt, like AGFRgui's PDB
panel, reads PDB files (runAGFR.profilePDBfile takes .pdb / .ent), so an entry
the RCSB distributes only as mmCIF (the largest structures) is reported rather
than fetched in a format ADapt cannot open.

Downloads are cached under <outdir>/_rcsb/, so trying the same entry again
does not touch the network.
"""
import gzip
import os
import re
import urllib.error
import urllib.request

CODE = re.compile(r'^[0-9][A-Za-z0-9]{3}$')
URL = 'https://files.rcsb.org/download/%s.pdb'
TIMEOUT = 30


class FetchError(Exception):
    pass


def looks_like_code(text):
    return bool(CODE.match((text or '').strip()))


def fetch(code, cachedir, timeout=TIMEOUT):
    """(path, url, was_cached) for the PDB-format file of entry `code`."""
    code = (code or '').strip().upper()
    if not looks_like_code(code):
        raise FetchError('%r is not a PDB id. They are four characters, the '
                         'first a digit, such as 6GGD, 1A3R or 4A1X.' % code)
    url = URL % code
    out = os.path.join(cachedir, '%s.pdb' % code)
    if os.path.isfile(out) and os.path.getsize(out) > 200:
        return out, url, True
    os.makedirs(cachedir, exist_ok=True)
    from . import APP, __version__
    req = urllib.request.Request(url, headers={'User-Agent': '%s/%s' % (APP, __version__)})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise FetchError(
                '%s has no PDB-format file at the RCSB. Either the id does not '
                'exist, or the entry is too large for the PDB format and is '
                'only distributed as mmCIF, which ADapt does not read (it reads '
                'PDB files, as AGFRgui\'s PDB panel does).' % code)
        raise FetchError('could not fetch %s from the RCSB: HTTP %s' % (code, exc.code))
    except Exception as exc:                                  # noqa: BLE001
        raise FetchError('could not fetch %s from the RCSB (%s). If this machine '
                         'has no route to the internet, download the file '
                         'elsewhere and open it with Open PDB.' % (code, exc))
    if data[:2] == b'\x1f\x8b':
        data = gzip.decompress(data)
    if len(data) < 200:
        raise FetchError('the RCSB returned an empty file for %s' % code)
    tmp = out + '.part'
    with open(tmp, 'wb') as fh:
        fh.write(data)
    os.replace(tmp, out)
    return out, url, False
