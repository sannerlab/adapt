# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Fetch from RCSB. The network tests skip (with a message) when offline."""
import json
import os

from common import check, outdir
from adapt import fetch


def _online():
    try:
        fetch.fetch('1CRN', outdir('rcsb_probe'))
        return True
    except fetch.FetchError as e:
        if 'could not fetch' in str(e):
            print('      skip: no route to the RCSB (%s)' % str(e)[:80])
            return False
        raise


def test_ids_are_checked_before_any_download():
    for bad in ('', 'ABCD', '1AB', '1ABCD', '1A-B', '../x'):
        try:
            fetch.fetch(bad, outdir('rcsb_bad'))
            raise AssertionError('%r accepted' % bad)
        except fetch.FetchError as e:
            check('is not a PDB id' in str(e), str(e))
    check(os.listdir(outdir('rcsb_bad')) == [])


def test_fetch_and_cache():
    if not _online():
        return
    cache = outdir('rcsb_cache')
    path, url, cached = fetch.fetch('1crn', cache)
    check(url == 'https://files.rcsb.org/download/1CRN.pdb' and not cached, (url, cached))
    check(open(path).read().startswith('HEADER'), 'a PDB file')
    path2, _, cached2 = fetch.fetch('1CRN', cache)
    check(path2 == path and cached2, 'second fetch served from the cache')


def test_mmcif_only_entry_is_explained():
    if not _online():
        return
    try:
        fetch.fetch('4V6X', outdir('rcsb_big'))           # a ribosome: mmCIF only
        raise AssertionError('4V6X has no PDB-format file')
    except fetch.FetchError as e:
        check('only distributed as mmCIF' in str(e), str(e))


def test_fetch_through_the_server():
    if not _online():
        return
    import threading
    import urllib.request
    from adapt.server import build
    httpd = build('127.0.0.1', 0, outdir('server_fetch'))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = 'http://127.0.0.1:%d' % httpd.server_address[1]
    try:
        req = urllib.request.Request(base + '/api/fetch?code=1crn&t=' + httpd.app.token,
                                     data=b'{}', method='POST')
        d = json.loads(urllib.request.urlopen(req).read().decode())
        check('error' not in d, d.get('error'))
        check(d['pdbid'] == '1crn' and d['tabs'], (d['pdbid'], d['tabs']))
        rec = json.load(open(os.path.join(d['folder'], 'session.json')))
        src = rec['input']['source']
        check(src['rcsb'].endswith('/1CRN.pdb') and rec['input']['name'] == '1crn.pdb', src)
        check(os.path.isfile(os.path.join(httpd.app.outdir, '_rcsb', '1CRN.pdb')))
    finally:
        httpd.shutdown()
