"""The HTTP layer, end to end on a real socket."""
import json
import os
import threading
import urllib.error
import urllib.request

from common import DATA, check, oracle_receptor, outdir, read


def _server():
    from adapt.server import build
    httpd = build('127.0.0.1', 0, outdir('server'))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd, 'http://127.0.0.1:%d' % httpd.server_address[1], httpd.app.token


def _get(url, raw=False):
    with urllib.request.urlopen(url) as r:
        body = r.read()
        return (r.status, body, dict(r.headers)) if raw else json.loads(body.decode())


def _post(url, obj=None, data=None):
    body = data if data is not None else json.dumps(obj or {}).encode()
    req = urllib.request.Request(url, data=body, method='POST')
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read().decode())


def _code(url, data=None):
    try:
        req = urllib.request.Request(url, data=data, method='POST' if data is not None else 'GET')
        with urllib.request.urlopen(req) as r:
            return r.status
    except urllib.error.HTTPError as e:
        return e.code


def test_token_is_required():
    httpd, base, tok = _server()
    try:
        check(_code(base + '/') == 403)
        check(_code(base + '/api/status') == 403)
        check(_code(base + '/api/upload?name=x.pdb', b'x') == 403)
        check(_code(base + '/?t=wrong') == 403)
        st, body, hdr = _get(base + '/?t=' + tok, raw=True)
        check(st == 200 and b'3Dmol-min.js' in body and b'app.js' in body)
        check(_code(base + '/static/app.js') == 200)
        check(_code(base + '/static/../server.py') == 404)
    finally:
        httpd.shutdown()


def test_upload_edit_save_download():
    httpd, base, tok = _server()
    t = '?t=' + tok
    try:
        st = _get(base + '/api/status' + t)
        check(st['missing'] == [], st['missing'])
        d = _post(base + '/api/upload' + t + '&name=1a3r.pdb',
                  data=read(os.path.join(DATA, '1a3r.pdb')))
        check('error' not in d, d.get('error'))
        sid = d['session']
        check(d['tabs'] == ['ASU and Biomol 1'] and d['current'] == d['tabs'][0])
        check(d['buttons'] == {'receptor': True, 'ligand': False, 'other': True}, d['buttons'])

        # the viewer's atom indices are the server's: one PDB atom line per atom
        s = _get(base + '/api/structure' + t + '&s=' + sid)
        n = sum(1 for l in s['pdb'].splitlines() if l.startswith(('ATOM', 'HETATM')))
        check(n == s['n'] == len(s['atoms']), (n, s['n'], len(s['atoms'])))

        # a residue made ligand and back: the ligand button follows
        def first(nodes, kind):
            for x in nodes:
                if x['kind'] == kind:
                    return x
                y = first(x.get('children', []), kind)
                if y:
                    return y
        res = first(d['tree']['roots'], 'residue')
        v = _post(base + '/api/click' + t, {'s': sid, 'node': res['id'], 'col': 'lig'})
        check(v['states'][str(res['id'])][:2] == ['unchecked', 'checked'], v['states'][str(res['id'])])
        check(v['buttons']['ligand'] is True)
        f = _get(base + '/api/focus' + t + '&s=%s&node=%d' % (sid, res['id']))
        check(len(f['indices']) > 0)
        v = _post(base + '/api/click' + t, {'s': sid, 'node': res['id'], 'col': 'rec'})
        check(v['buttons']['ligand'] is False)

        bad = _code(base + '/api/click' + t, json.dumps({'s': sid, 'node': 99999999, 'col': 'rec'}).encode())
        check(bad == 400, bad)

        sv = _post(base + '/api/save' + t, {'s': sid, 'kind': 'receptor'})
        check('error' not in sv, sv.get('error'))
        st, body, hdr = _get(base + '/api/download' + t + '&s=%s&save=%d' % (sid, sv['save']), raw=True)
        check('1a3r_rec.pdbqt' in hdr.get('Content-Disposition', ''), hdr)
        ref = oracle_receptor('1a3r.pdb', outdir('oracle_server'))
        check(body == read(ref), 'downloaded receptor differs from agfr --toPdbqt')

        no = _code(base + '/api/save' + t, json.dumps({'s': sid, 'kind': 'ligand'}).encode())
        check(no == 400, 'saving a ligand with no ligand atoms is refused (%s)' % no)
    finally:
        httpd.shutdown()
