"""A minimal Chrome DevTools Protocol client, standard library only.

Used by docs/make_screenshots.py to drive a headless Chrome against a running
ADapt server. There is no websocket package in the environments this runs in,
so the websocket handshake and framing (RFC 6455, client side) are done here.
"""
import base64
import json
import os
import socket
import struct
import subprocess
import time
import urllib.request


class WS(object):
    def __init__(self, url):
        assert url.startswith('ws://')
        hostport, path = url[5:].split('/', 1)
        host, port = hostport.split(':')
        self.s = socket.create_connection((host, int(port)))
        key = base64.b64encode(os.urandom(16)).decode()
        self.s.sendall(('GET /%s HTTP/1.1\r\nHost: %s\r\nUpgrade: websocket\r\n'
                        'Connection: Upgrade\r\nSec-WebSocket-Key: %s\r\n'
                        'Sec-WebSocket-Version: 13\r\n\r\n' % (path, hostport, key)).encode())
        buf = b''
        while b'\r\n\r\n' not in buf:
            buf += self.s.recv(4096)
        self.rest = buf.split(b'\r\n\r\n', 1)[1]
        self.id = 0
        self.events = []

    def _exact(self, n):
        while len(self.rest) < n:
            chunk = self.s.recv(1 << 20)
            if not chunk:
                raise EOFError('chrome closed the connection')
            self.rest += chunk
        out, self.rest = self.rest[:n], self.rest[n:]
        return out

    def recv(self):
        data = b''
        while True:
            h = self._exact(2)
            fin, ln = h[0] & 0x80, h[1] & 0x7f
            if ln == 126:
                ln = struct.unpack('>H', self._exact(2))[0]
            elif ln == 127:
                ln = struct.unpack('>Q', self._exact(8))[0]
            data += self._exact(ln)
            if fin:
                return data.decode()

    def send(self, obj):
        payload = json.dumps(obj).encode()
        mask = os.urandom(4)
        n = len(payload)
        hdr = bytes([0x81])
        if n < 126:
            hdr += bytes([0x80 | n])
        elif n < 65536:
            hdr += bytes([0x80 | 126]) + struct.pack('>H', n)
        else:
            hdr += bytes([0x80 | 127]) + struct.pack('>Q', n)
        self.s.sendall(hdr + mask + bytes(b ^ mask[i % 4] for i, b in enumerate(payload)))

    def call(self, method, **params):
        self.id += 1
        me = self.id
        self.send({'id': me, 'method': method, 'params': params})
        while True:
            m = json.loads(self.recv())
            if m.get('id') == me:
                if 'error' in m:
                    raise RuntimeError('%s: %s' % (method, m['error']))
                return m.get('result', {})
            self.events.append(m)


def start_chrome(profile, port=9345, width=1440, height=960):
    """Headless Chrome with software WebGL (the 3D view needs it)."""
    p = subprocess.Popen(
        ['google-chrome', '--headless=new', '--remote-debugging-port=%d' % port,
         '--user-data-dir=%s' % profile, '--window-size=%d,%d' % (width, height),
         '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--hide-scrollbars',
         '--no-first-run', '--no-default-browser-check', 'about:blank'],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(100):
        try:
            tabs = json.load(urllib.request.urlopen('http://127.0.0.1:%d/json' % port))
            page = [t for t in tabs if t['type'] == 'page'][0]
            return p, WS(page['webSocketDebuggerUrl'])
        except Exception:
            time.sleep(0.2)
    p.terminate()
    raise RuntimeError('chrome did not start')


def js(ws, expr):
    r = ws.call('Runtime.evaluate', expression=expr, awaitPromise=True, returnByValue=True)
    if 'exceptionDetails' in r:
        raise RuntimeError('JS error in %r: %s' % (expr[:120], r['exceptionDetails']))
    return r['result'].get('value')


def wait(ws, expr, timeout=120, what=None):
    t = time.time()
    while time.time() - t < timeout:
        if js(ws, expr):
            return True
        time.sleep(0.25)
    raise RuntimeError('timeout waiting for %s' % (what or expr))


def screenshot(ws, path, clip=None):
    kw = {'format': 'png'}
    if clip:
        kw['clip'] = dict(clip, scale=1)
    r = ws.call('Page.captureScreenshot', **kw)
    with open(path, 'wb') as fh:
        fh.write(base64.b64decode(r['data']))


def js_errors(ws):
    return [e for e in ws.events if e.get('method') == 'Runtime.exceptionThrown' or
            (e.get('method') == 'Runtime.consoleAPICalled' and e['params']['type'] == 'error')]
