# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The HTTP layer, end to end on a real socket."""
import json
import os
import re
import threading
import urllib.error
import urllib.request

from common import DATA, check, oracle_receptor, outdir, read, without_biomolecule


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
        names = [f['name'] for f in sv['files']]
        check(names == ['1a3r_rec.pdbqt', '1a3r_rec_H.pdb'], names)
        st, body, hdr = _get(base + '/api/download' + t + '&s=%s&save=%d&file=1a3r_rec.pdbqt'
                             % (sid, sv['save']), raw=True)
        check('1a3r_rec.pdbqt' in hdr.get('Content-Disposition', ''), hdr)
        ref = oracle_receptor('1a3r.pdb', outdir('oracle_server'))
        check(without_biomolecule(body) == read(ref), 'downloaded receptor differs from agfr --toPdbqt')
        check(b'REMARK 850 AGFR biomolecule ASU and Biomol 1: chains H L P' in body, 'biomolecule recorded')

        no = _code(base + '/api/save' + t, json.dumps({'s': sid, 'kind': 'ligand'}).encode())
        check(no == 400, 'saving a ligand with no ligand atoms is refused (%s)' % no)
    finally:
        httpd.shutdown()


def test_docs_are_served():
    """The documentation is public (no token); only docs/ pages, style and images."""
    httpd, base, tok = _server()
    try:
        st, body, hdr = _get(base + '/docs/', raw=True)
        check(st == 200 and b'<h1>ADapt</h1>' in body and 'text/html' in hdr['Content-Type'])
        for path in ('/docs/case-biomolecule.html', '/docs/style.css', '/docs/img/overview.png'):
            check(_code(base + path) == 200, path)
        for bad in ('/docs/../README.md', '/docs/%2e%2e/README.md', '/docs/build_docs.py',
                    '/docs/src/index.html', '/docs/tools/cdp.py', '/docs/img/../index.html'):
            check(_code(base + bad) == 404, bad)
        st, body, hdr = _get(base + '/?t=' + tok, raw=True)
        check(b'href="docs/index.html"' in body, 'the page links the docs')
    finally:
        httpd.shutdown()


def test_open_a_session_from_elsewhere():
    """A session made through /api/upload (by another tool) can be opened:
    /api/session returns what /api/upload returned, plus its saves;
    /api/sessions lists it; &file= downloads one file of a save."""
    httpd, base, tok = _server()
    t = '?t=' + tok
    try:
        up = _post(base + '/api/upload' + t + '&name=6ggd.pdb', data=read(os.path.join(DATA, '6ggd.pdb')))
        sid = up['session']
        sv = _post(base + '/api/save' + t, {'s': sid, 'kind': 'receptor'})
        d = _get(base + '/api/session' + t + '&s=' + sid)
        check(d['session'] == sid and d['tabs'] == up['tabs'] and d['tree'] == up['tree'], 'same session')
        check(d['saves'] == [{'save': sv['save'], 'kind': 'receptor', 'tab': up['current'],
                              'files': sv['files']}], d['saves'])
        lst = _get(base + '/api/sessions' + t)['sessions']
        check([x['session'] for x in lst] == [sid] and lst[0]['pdbid'] == '6ggd', lst)
        check(_code(base + '/api/session' + t + '&s=nosuchsession') == 400)
        for f in sv['files']:
            st, body, hdr = _get(base + '/api/download' + t + '&s=%s&save=%d&file=%s'
                                 % (sid, sv['save'], f['name']), raw=True)
            check(f['name'] in hdr['Content-Disposition'] and body, f['name'])
        st, body, hdr = _get(base + '/api/download' + t + '&s=%s&save=%d' % (sid, sv['save']), raw=True)
        check(hdr['Content-Type'] == 'application/zip', 'no file: all of them, zipped')
        check(_code(base + '/api/download' + t + '&s=%s&save=%d&file=../x' % (sid, sv['save'])) == 400)
    finally:
        httpd.shutdown()


# ADCP_docking's proxy (adcp_dock/server.py, _rebase_adapt) rewrites these to
# serve ADapt under /adapt/: with relative paths there must be nothing to rewrite
_ABS_ATTR = re.compile(rb"""((?:href|src|action)=["'])/(?!/)""")
_ABS_JS = re.compile(rb"""(['"`])/(api|static|docs)/""")
_ABS_CSS = re.compile(rb"""url\((["']?)/(?!/)""")


def test_page_uses_relative_paths():
    httpd, base, tok = _server()
    try:
        st, html, h = _get(base + '/?t=' + tok, raw=True)
        check(not _ABS_ATTR.findall(html) and not _ABS_JS.findall(html), _ABS_ATTR.findall(html))
        for f in ('app.js', 'app.css'):
            st, body, h = _get(base + '/static/' + f, raw=True)
            check(not _ABS_JS.findall(body) and not _ABS_CSS.findall(body), (f, _ABS_JS.findall(body)))
    finally:
        httpd.shutdown()


def _prefix_proxy(target):
    """A reverse proxy serving `target` under /adapt/, WITHOUT rewriting bodies
    (it strips the prefix from the request path, as ADCP_docking's does)."""
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

    class P(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _go(self, method):
            if not (self.path.startswith('/adapt/') or self.path.startswith('/adapt?')):
                self.send_response(404); self.send_header('Content-Length', '0'); self.end_headers()
                return
            rest = self.path[len('/adapt'):]
            if rest.startswith('?'):
                rest = '/' + rest
            n = int(self.headers.get('Content-Length') or 0)
            data = self.rfile.read(n) if method == 'POST' else None
            hd = {k: self.headers[k] for k in ('Content-Type',) if self.headers.get(k)}
            req = urllib.request.Request(target + rest, data=data, method=method, headers=hd)
            try:
                r = urllib.request.urlopen(req)
                code, body, rh = r.status, r.read(), r.headers
            except urllib.error.HTTPError as e:
                code, body, rh = e.code, e.read(), e.headers
            self.send_response(code)
            for k in ('Content-Type', 'Content-Disposition'):
                if rh.get(k):
                    self.send_header(k, rh.get(k))
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            self._go('GET')

        def do_POST(self):
            self._go('POST')

    p = ThreadingHTTPServer(('127.0.0.1', 0), P)
    threading.Thread(target=p.serve_forever, daemon=True).start()
    return p, 'http://127.0.0.1:%d' % p.server_address[1]


def test_served_under_a_prefix_by_a_plain_proxy():
    """Behind a proxy that maps /adapt/... to ADapt and rewrites nothing, the
    page's relative paths reach ADapt: every asset the page names loads, and
    upload / save / download work through /adapt/api/..."""
    from urllib.parse import urljoin
    httpd, base, tok = _server()
    proxy, pbase = _prefix_proxy(base)
    try:
        page_url = pbase + '/adapt/?t=' + tok
        st, html, h = _get(page_url, raw=True)
        check(st == 200 and b'static/app.js' in html, 'page through the proxy')
        for ref in re.findall(rb'(?:href|src)="([^"#]+)"', html):
            url = urljoin(page_url, ref.decode())
            check(url.startswith(pbase + '/adapt/'), url)
            check(_code(url) == 200, url)
        api = urljoin(page_url, 'api/')                 # what app.js resolves to
        up = _post(api + 'upload?t=%s&name=1a3r.pdb' % tok, data=read(os.path.join(DATA, '1a3r.pdb')))
        sv = _post(api + 'save?t=' + tok, {'s': up['session'], 'kind': 'receptor'})
        st, body, hdr = _get(api + 'download?t=%s&s=%s&save=%d&file=1a3r_rec.pdbqt'
                             % (tok, up['session'], sv['save']), raw=True)
        check(body.startswith(b'REMARK 850 AGFR version'), body[:40])
    finally:
        proxy.shutdown()
        httpd.shutdown()


def test_base_path_option():
    """--base-path /adapt: prefixed requests are served, /adapt redirects to
    /adapt/ (relative paths need the slash), and / still works."""
    from adapt.server import build, clean_base_path
    check(clean_base_path('adapt/') == '/adapt' and clean_base_path('') == '')
    httpd = build('127.0.0.1', 0, outdir('server_base'), base_path='/adapt')
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base, tok = 'http://127.0.0.1:%d' % httpd.server_address[1], httpd.app.token
    try:
        check(_code(base + '/adapt/?t=' + tok) == 200)
        check(_code(base + '/adapt/static/app.js') == 200)
        check(_get(base + '/adapt/api/status?t=' + tok)['missing'] == [])
        check(_code(base + '/?t=' + tok) == 200, 'no prefix still works')
        import http.client
        c = http.client.HTTPConnection('127.0.0.1', httpd.server_address[1])
        c.request('GET', '/adapt?t=' + tok)
        r = c.getresponse()
        check(r.status == 301 and r.getheader('Location') == '/adapt/?t=' + tok, (r.status, r.getheader('Location')))
    finally:
        httpd.shutdown()
