# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Desktop use: open the page in the browser, reuse a server that is already
running, and the Linux menu entry (a .desktop file) that starts ADapt with one
click.

A server started on this machine records itself in ~/.adapt/server.json (mode
600: the record holds the page's token; $ADAPT_STATE_DIR moves it). `adapt
serve --open` reads that record first: when the server it names still answers,
a new browser tab is opened on it and no second server is started.
"""
import json
import os
import sys
import urllib.parse
import urllib.request
import webbrowser

LAUNCHER = 'adapt.desktop'


def state_dir():
    return os.path.expanduser(os.environ.get('ADAPT_STATE_DIR', '~/.adapt'))


def state_file():
    return os.path.join(state_dir(), 'server.json')


def record(url, pid, outdir, version):
    """Write the record of the server now running (atomically, mode 600)."""
    d = state_dir()
    os.makedirs(d, mode=0o700, exist_ok=True)
    tmp = '%s.%d' % (state_file(), pid)
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, 'w') as f:
        json.dump({'url': url, 'pid': pid, 'outdir': outdir, 'version': version}, f)
    os.replace(tmp, state_file())


def forget(pid):
    """Remove the record if it is this server's (another may have replaced it)."""
    try:
        with open(state_file()) as f:
            rec = json.load(f)
    except (OSError, ValueError):
        return
    if rec.get('pid') == pid:
        try:
            os.remove(state_file())
        except OSError:
            pass


def _no_proxy_get(url, timeout):
    # the server is on this machine: never through an http_proxy
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(url, timeout=timeout) as r:
        return json.loads(r.read().decode())


def running(outdir=None):
    """The page URL of the recorded server if it still answers (and keeps its
    outputs in `outdir`, when given); else None."""
    try:
        with open(state_file()) as f:
            rec = json.load(f)
    except (OSError, ValueError):
        return None
    url = rec.get('url') or ''
    if outdir and os.path.abspath(os.path.expanduser(outdir)) != rec.get('outdir'):
        return None
    u = urllib.parse.urlparse(url)
    if not u.scheme or not u.query:
        return None
    probe = urllib.parse.urlunparse(u._replace(path=u.path.rstrip('/') + '/api/ping'))
    try:
        d = _no_proxy_get(probe, 3)
    except Exception:
        return None
    return url if d.get('ok') else None


def can_open_browser():
    """A graphical session to open a browser in. Without one (an SSH login, a
    cluster node) webbrowser could start a text browser in the terminal."""
    if sys.platform in ('darwin', 'win32'):
        return True
    return bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY'))


def open_browser(url):
    """Open `url` in a new tab of the default browser; False when there is no
    graphical session or no browser."""
    if not can_open_browser():
        return False
    try:
        return webbrowser.open(url, new=2)
    except Exception:
        return False


# ---------------------------------------------------------------- launcher
def _q(arg):
    """One argument of a .desktop Exec line, quoted as the spec asks."""
    if arg and not any(c in arg for c in ' \t\n"\'\\><~|&;$*?#()`'):
        return arg
    return '"%s"' % arg.replace('\\', '\\\\').replace('"', '\\"').replace('`', '\\`').replace('$', '\\$')


def launcher_text(adapt_cmd, env_root, outdir, log, icon):
    """The .desktop file: starts the server (or reuses the running one) and
    opens the page. The environment and output folder found now are written
    in, since a desktop session does not read the shell's start-up files."""
    exe = ['env', 'ADCPSUITE=' + env_root, adapt_cmd, 'serve', '--open',
           '--outdir', outdir, '--log', log]
    return '\n'.join([
        '[Desktop Entry]',
        'Type=Application',
        'Version=1.0',
        'Name=ADapt',
        'GenericName=AutoDock Analysis and Preparation Tool',
        'Comment=Prepare PDB entries for docking with AutoDock, in a web browser',
        'Exec=' + ' '.join(_q(a) for a in exe),
        'Icon=' + icon,
        'Terminal=false',
        'StartupNotify=false',
        'Categories=Science;Chemistry;Biology;',
        'Keywords=docking;AutoDock;PDB;receptor;PDBQT;',
        ''])


def applications_dir():
    base = os.environ.get('XDG_DATA_HOME') or os.path.join(os.path.expanduser('~'), '.local', 'share')
    return os.path.join(base, 'applications')


def desktop_dir():
    """$XDG_DESKTOP_DIR, else the desktop folder named in ~/.config/user-dirs.dirs
    (a localized one, e.g. ~/Bureau), else ~/Desktop."""
    home = os.path.expanduser('~')
    d = os.environ.get('XDG_DESKTOP_DIR')
    if not d:
        conf = os.path.join(os.environ.get('XDG_CONFIG_HOME') or os.path.join(home, '.config'),
                            'user-dirs.dirs')
        try:
            with open(conf) as f:
                for line in f:
                    if line.startswith('XDG_DESKTOP_DIR='):
                        d = line.split('=', 1)[1].strip().strip('"').replace('$HOME', home)
        except OSError:
            pass
    return d or os.path.join(home, 'Desktop')


def install_launcher(text, on_desktop=False):
    """Write the menu entry (and a copy on ~/Desktop); the paths written."""
    paths = [os.path.join(applications_dir(), LAUNCHER)]
    if on_desktop:
        paths.append(os.path.join(desktop_dir(), LAUNCHER))
    for p in paths:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, 'w') as f:
            f.write(text)
        os.chmod(p, 0o755)
    return paths


def remove_launcher():
    gone = []
    for p in (os.path.join(applications_dir(), LAUNCHER), os.path.join(desktop_dir(), LAUNCHER)):
        if os.path.exists(p):
            os.remove(p)
            gone.append(p)
    return gone
