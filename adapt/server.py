"""ADapt's web front end: a local page that does what AGFRgui's PDB panel does.

The server pattern (standard library only, bound to localhost, a random token
in the URL, JSON under /api/*, bundled static files) follows
ADCP_docking/adcp_dock/server.py by Aadhil Haq (BSD 3-Clause). A second
person on the same machine cannot drive someone else's session by guessing
the port, and nothing is exposed beyond localhost unless --host says so;
remote use goes through an SSH tunnel.
"""
import json
import os
import secrets
import threading
import traceback
import urllib.parse
import zipfile
import io
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from . import APP, APP_LONG, __version__, env, fetch, page
from .flags import Refused

MAX_BODY = 1 << 20
MAX_UPLOAD = 256 << 20      # a large assembly PDB can be tens of MB

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
# the documentation (docs/ of the repository): public pages, served without the token
DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'docs')
DOC_TYPES = {'.html': 'text/html; charset=utf-8', '.css': 'text/css', '.png': 'image/png'}


class App(object):
    def __init__(self, outdir):
        self.outdir = os.path.abspath(os.path.expanduser(outdir))
        os.makedirs(self.outdir, exist_ok=True)
        self.token = secrets.token_urlsafe(18)
        self.sessions = {}
        self._lock = threading.Lock()
        self._prov = None

    def provenance(self):
        if self._prov is None:
            self._prov = env.provenance()
        return self._prov

    def status(self):
        return {'version': __version__, 'missing': env.check(),
                'outdir': self.outdir, 'environment': self.provenance()}

    def upload(self, name, data, source=None):
        from .session import Session
        if not data:
            raise Refused('the file is empty')
        sid = secrets.token_hex(8)
        s = Session(sid, data, name, self.outdir, self.provenance(), source)
        with self._lock:
            self.sessions[sid] = s
        d = s.info()
        d.update(s.view(full=True))
        return d

    def fetch(self, code):
        """Download an entry from the RCSB (cached) and open it like an upload."""
        import datetime
        try:
            path, url, cached = fetch.fetch(code, os.path.join(self.outdir, '_rcsb'))
        except fetch.FetchError as e:
            raise Refused(str(e))
        with open(path, 'rb') as fh:
            data = fh.read()
        # lower case, as the RCSB names its files: PDBProfiler takes the PDB
        # id from the file name and looks the entry's ligands up by it
        name = os.path.basename(path).lower()
        d = self.upload(name, data, source={
            'rcsb': url, 'cache': path, 'from_cache': cached,
            'at': datetime.datetime.now().isoformat(timespec='seconds')})
        d['fetched'] = {'url': url, 'from_cache': cached}
        return d

    def session(self, sid):
        s = self.sessions.get(sid or '')
        if s is None:
            raise Refused('no such session; load the PDB file again')
        return s


class Handler(BaseHTTPRequestHandler):
    server_version = '%s/%s' % (APP, __version__)

    def log_message(self, *a):
        pass

    # ---- plumbing --------------------------------------------------------
    def _send(self, obj, code=200, ctype='application/json', headers=()):
        body = (obj if isinstance(obj, bytes)
                else json.dumps(obj).encode() if ctype == 'application/json'
                else obj.encode())
        self.send_response(code)
        self.send_header('Content-Type', ctype + ('; charset=utf-8' if not isinstance(obj, bytes) else ''))
        self.send_header('Content-Length', str(len(body)))
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Cache-Control', 'no-store')
        for k, v in headers:
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _auth(self, q):
        return secrets.compare_digest(
            q.get('t', [''])[0] or self.headers.get('X-Token', ''),
            self.server.app.token)

    def _raw(self, limit):
        n = int(self.headers.get('Content-Length') or 0)
        if n <= 0:
            return b''
        if n > limit:
            raise Refused('the request is larger than %d MB' % (limit >> 20))
        return self.rfile.read(n)

    def _json(self):
        raw = self._raw(MAX_BODY)
        try:
            return json.loads(raw.decode()) if raw else {}
        except ValueError:
            raise Refused('the request body is not JSON')

    def _run(self, fn):
        try:
            return self._send(fn())
        except Refused as e:
            return self._send({'error': str(e)}, 400)
        except ValueError as e:
            return self._send({'error': str(e)}, 400)
        except Exception as e:
            return self._send({'error': '%s: %s' % (type(e).__name__, e),
                               'trace': traceback.format_exc()}, 500)

    # ---- GET ---------------------------------------------------------------
    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        app = self.server.app
        if u.path.startswith('/static/'):
            return self._static(u.path)
        if u.path == '/docs' or u.path.startswith('/docs/'):
            return self._docs(u.path)
        if not self._auth(q):
            if u.path == '/':
                return self._send('open the URL this tool printed; it carries the token',
                                  403, 'text/plain')
            return self._send({'error': 'bad or missing token'}, 403)
        g = lambda k, d=None: q.get(k, [d])[0]
        if u.path == '/':
            return self._send(page.render(), 200, 'text/html')
        if u.path == '/api/status':
            return self._run(app.status)
        if u.path == '/api/structure':
            return self._run(lambda: app.session(g('s')).structure(g('tab')))
        if u.path == '/api/focus':
            return self._run(lambda: app.session(g('s')).focus(g('node')))
        if u.path == '/api/mutateopts':
            return self._run(lambda: app.session(g('s')).mutate_options(g('node')))
        if u.path == '/api/rotamer':
            return self._run(lambda: app.session(g('s')).rotamer(g('node'), g('choice', '-1')))
        if u.path == '/api/download':
            return self._download(app, g('s'), g('save'))
        return self._send({'error': 'no such endpoint'}, 404)

    def _static(self, path):
        name = os.path.basename(path)
        full = os.path.join(STATIC, name)
        if os.path.dirname(os.path.abspath(full)) != os.path.abspath(STATIC) \
                or not os.path.isfile(full):
            return self._send({'error': 'no such file'}, 404)
        ctype = ('application/javascript' if name.endswith('.js')
                 else 'text/css' if name.endswith('.css') else 'text/plain')
        with open(full, 'rb') as fh:
            body = fh.read()
        self.send_response(200)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'max-age=3600')
        self.end_headers()
        self.wfile.write(body)

    def _docs(self, path):
        """docs/*.html, docs/style.css and docs/img/*.png; nothing else under docs/."""
        rel = path[len('/docs'):].lstrip('/') or 'index.html'
        parts = rel.split('/')
        ok = (len(parts) == 1 and parts[0].endswith(('.html', '.css'))) or \
             (len(parts) == 2 and parts[0] == 'img' and parts[1].endswith('.png'))
        full = os.path.abspath(os.path.join(DOCS, *parts))
        if not ok or not full.startswith(os.path.abspath(DOCS) + os.sep) \
                or '..' in parts or not os.path.isfile(full):
            return self._send({'error': 'no such page'}, 404)
        with open(full, 'rb') as fh:
            body = fh.read()
        self.send_response(200)
        self.send_header('Content-Type', DOC_TYPES[os.path.splitext(full)[1]])
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'max-age=600')
        self.end_headers()
        self.wfile.write(body)

    def _download(self, app, sid, n):
        try:
            s = app.session(sid)
            entry, paths = s.save_paths(n)
        except (Refused, ValueError) as e:
            return self._send({'error': str(e)}, 400)
        if len(paths) == 1:
            with open(paths[0], 'rb') as fh:
                body = fh.read()
            name = os.path.basename(paths[0])
            ctype = 'chemical/x-pdbqt'
        else:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as z:
                for p in paths:
                    z.write(p, os.path.basename(p))
            body = buf.getvalue()
            name = '%s_%s.zip' % (s.profile._pdbid, entry['kind'][:3])
            ctype = 'application/zip'
        return self._send(body, 200, ctype, headers=[
            ('Content-Disposition', 'attachment; filename="%s"' % name)])

    # ---- POST --------------------------------------------------------------
    def do_POST(self):
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        if not self._auth(q):
            return self._send({'error': 'bad or missing token'}, 403)
        app = self.server.app
        if u.path == '/api/upload':
            return self._run(lambda: app.upload(q.get('name', [''])[0],
                                                self._raw(MAX_UPLOAD)))
        if u.path == '/api/fetch':
            return self._run(lambda: app.fetch(q.get('code', [''])[0]))
        try:
            body = self._json()
        except Refused as e:
            return self._send({'error': str(e)}, 400)
        s = lambda: app.session(body.get('s'))
        routes = {
            '/api/tab': lambda: _tab(s(), body.get('tab')),
            '/api/click': lambda: s().click(body.get('node'), body.get('col')),
            '/api/altloc': lambda: s().altloc(body.get('node'), body.get('index')),
            '/api/charge': lambda: s().charge(body.get('node'), body.get('value')),
            '/api/sidechain': lambda: s().sidechain(body.get('node'), body.get('act')),
            '/api/mutate': lambda: s().mutate(body.get('node'), body.get('choice', '')),
            '/api/select': lambda: s().select(body.get('what')),
            '/api/save': lambda: _save(s(), body.get('kind')),
        }
        fn = routes.get(u.path)
        if fn is None:
            return self._send({'error': 'no such endpoint'}, 404)
        return self._run(fn)


def _tab(s, name):
    s.switch(name)
    return s.view(full=True)


def _save(s, kind):
    d = s.save(kind)
    d.update(s.view())      # its log: the save's lines, already in s.log
    return d


def build(host='127.0.0.1', port=8766, outdir='~/adapt_runs'):
    httpd = ThreadingHTTPServer((host, port), Handler)
    httpd.daemon_threads = True
    httpd.app = App(outdir)
    return httpd


def serve(host='127.0.0.1', port=8766, outdir='~/adapt_runs'):
    httpd = build(host, port, outdir)
    bound = httpd.server_address[1]
    url = 'http://%s:%d/?t=%s' % (host, bound, httpd.app.token)
    print('%s %s (%s) is at\n\n    %s\n' % (APP, __version__, APP_LONG, url), flush=True)
    print('outputs are kept in %s' % httpd.app.outdir, flush=True)
    miss = env.check()
    if miss:
        print('\n%s cannot run yet. NOT FOUND: %s' % (APP, '; '.join(miss)), flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print('\n%s stopped' % APP)
