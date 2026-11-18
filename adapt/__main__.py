# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Command line: python -m adapt serve | check | examples | test (bin/adapt wraps it)."""
import argparse
import os
import sys

DEFAULT_OUTDIR = os.environ.get('ADAPT_OUTDIR', '~/adapt_runs')


def main(argv=None):
    from . import APP, APP_LONG, COPYRIGHT, __version__
    p = argparse.ArgumentParser(
        prog='adapt',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='%s %s: %s\n\n'
                    'Prepares receptor and ligand PDBQT files from a PDB entry in a web\n'
                    'browser, the way AGFRgui\'s PDB panel does, with ADFR\'s own code.'
                    % (APP, __version__, APP_LONG),
        epilog='examples:\n'
               '  %(prog)s serve                   start %(app)s on port 8766 and print its URL\n'
               '  %(prog)s serve --port 9000       another fixed port (the browser remembers\n'
               '                                folders and theme per port)\n'
               '  %(prog)s serve --open            open the page in the browser (or a new tab\n'
               '                                on the server already running); Quit in the\n'
               '                                page, or 10 idle minutes, stops it\n'
               '  %(prog)s launcher --desktop      add ADapt to the applications menu and the\n'
               '                                desktop (Linux): one click starts it\n'
               '  %(prog)s ligands --update        rebuild the ligand table from the RCSB\n'
               '  %(prog)s check                   "ready", or what %(app)s is missing\n'
               '  %(prog)s examples                where ADFR keeps its example PDB files\n'
               '  %(prog)s test -k parity          run part of the test suite\n\n'
               '%(app)s runs ADFR from the ADCP distribution (the adcpsuite env); the\n'
               'wrapper bin/adapt sets that env up. README.md has a first session.\n\n'
               '%(copyright)s'
               % {'prog': 'adapt', 'app': APP, 'copyright': COPYRIGHT})
    p.add_argument('--version', action='version', version='%s %s\n%s' % (APP, __version__, COPYRIGHT))
    sub = p.add_subparsers(dest='cmd')
    s = sub.add_parser('serve', help='start %s in a browser page (prints its URL)' % APP)
    s.add_argument('--host', default='127.0.0.1',
                   help='address to bind (default 127.0.0.1: this machine only)')
    s.add_argument('--port', type=int, default=8766, help='0 picks a free port')
    s.add_argument('--base-path', default=os.environ.get('ADAPT_BASE_PATH', ''),
                   help='serve under this path prefix as well, e.g. /adapt when a '
                        'proxy forwards /adapt/... unchanged (the page itself uses '
                        'relative paths, so a proxy that strips the prefix needs '
                        'nothing); default none, or $ADAPT_BASE_PATH')
    s.add_argument('--outdir', default=DEFAULT_OUTDIR,
                   help='where every session keeps its input, outputs and '
                        'session.json (default %(default)s); keep it on NFS, '
                        'never /tmp')
    s.add_argument('--open', action='store_true',
                   help='open the page in the default browser; when a server '
                        'started on this machine for the same output folder still '
                        'runs, open a new tab on it instead of starting another; '
                        'a busy port is replaced by a free one; the page gets a '
                        'Quit button')
    s.add_argument('--idle-exit', type=int, default=None, metavar='MIN',
                   help='stop after MIN minutes without a request (an open page '
                        'checks in every minute); default 10 with --open, else never')
    s.add_argument('--log', metavar='FILE',
                   help='append the server\'s messages to FILE instead of the terminal '
                        '(for a start from the desktop, which has no terminal)')
    la = sub.add_parser('launcher', help='add %s to the applications menu (Linux)' % APP)
    la.add_argument('--desktop', action='store_true', help='also put the icon on ~/Desktop')
    la.add_argument('--remove', action='store_true', help='remove the menu entry and the icon')
    la.add_argument('--outdir', default=DEFAULT_OUTDIR,
                    help='output folder of the sessions it starts (default %(default)s)')
    lg = sub.add_parser('ligands', help="the ligand table: which components each PDB entry holds")
    lg.add_argument('--update', action='store_true',
                    help='rebuild it from the RCSB PDB (a few minutes), for entries released '
                         'since it was built; kept in ~/.adapt/ligands/')
    sub.add_parser('check', help='report what %s needs that is missing' % APP)
    sub.add_parser('examples', help="list the example PDB files that ship with ADFR")
    t = sub.add_parser('test', help='run the %s test suite (tests/run_tests.py)' % APP)
    t.add_argument('-k', default='', help='only tests whose name contains this')
    a = p.parse_args(argv)
    if a.cmd == 'test':
        import subprocess
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        runner = os.path.join(here, 'tests', 'run_tests.py')
        return subprocess.call([sys.executable, runner] + (['-k', a.k] if a.k else []))
    from . import env
    if a.cmd == 'examples':
        # ADFR's own test inputs, in whichever ADFR installation is in use
        import ADFR
        d = os.path.join(os.path.dirname(os.path.abspath(ADFR.__file__)), 'tests', 'Data')
        print(d)
        for f in sorted(os.listdir(d)):
            if f.endswith('.pdb'):
                print('  ' + f)
        return 0
    if a.cmd == 'launcher':
        return launcher(a)
    if a.cmd == 'ligands':
        from . import ligtable
        if a.update:
            t = ligtable.build(progress=lambda m: print('  ' + m, flush=True))
            print('rebuilt: %s' % t.path)
        t = ligtable.current()
        d = t.describe()
        print('ligand table built %s: %d PDB entries, %d components (%s)'
              % (d['date'], d['entries'], d['codes'], 'shipped with ADapt' if d['shipped'] else t.path))
        return 0
    if a.cmd == 'check':
        miss = env.check()
        print('%s is ready' % APP if not miss else '%s cannot run: missing %s' % (APP, '; '.join(miss)))
        return 1 if miss else 0
    if a.cmd != 'serve':
        p.print_help()
        return 2
    if os.path.realpath(os.path.expanduser(a.outdir)).startswith('/tmp'):
        p.error('--outdir must not be under /tmp; it is wiped')
    if a.log:   # before the imports below, so their warnings go there too
        os.makedirs(os.path.dirname(os.path.abspath(os.path.expanduser(a.log))), exist_ok=True)
        out = open(os.path.expanduser(a.log), 'a', buffering=1)
        sys.stdout = sys.stderr = out
        import datetime
        print('---- %s' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), flush=True)
    from .server import serve
    from .server import clean_base_path
    try:
        base = clean_base_path(a.base_path)
    except ValueError as e:
        p.error(str(e))
    idle = a.idle_exit if a.idle_exit is not None else (10 if a.open else 0)
    serve(a.host, a.port, a.outdir, base, open_browser=a.open, idle_exit=idle)
    return 0


def launcher(a):
    """Write (or remove) the Linux menu entry that runs `adapt serve --open`
    with this environment and output folder."""
    from . import APP, desktop
    if a.remove:
        gone = desktop.remove_launcher()
        print('\n'.join('removed ' + g for g in gone) or 'no %s launcher was installed' % APP)
        return 0
    if not sys.platform.startswith('linux'):
        print('the launcher is for Linux desktops; elsewhere run: bin/adapt serve --open')
        return 2
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env_root = os.environ.get('ADS_ROOT') or sys.prefix
    outdir = os.path.abspath(os.path.expanduser(a.outdir))
    if outdir.startswith('/tmp'):
        print('--outdir must not be under /tmp; it is wiped')
        return 2
    text = desktop.launcher_text(os.path.join(here, 'bin', 'adapt'), env_root, outdir,
                                 os.path.join(desktop.state_dir(), 'server.log'),
                                 os.path.join(here, 'adapt', 'static', 'adapt-icon.png'))
    for p in desktop.install_launcher(text, a.desktop):
        print('wrote ' + p)
    if a.desktop:
        print('On GNOME, right-click the desktop icon once and choose "Allow Launching".')
    print('%s now starts from the applications menu: it opens in the browser, and Quit '
          'in the page stops it.' % APP)
    return 0


if __name__ == '__main__':
    sys.exit(main())
