"""Command line: python -m agfr_web serve | check | examples | test (bin/agfr-web wraps it)."""
import argparse
import os
import sys

DEFAULT_OUTDIR = os.environ.get('AGFR_WEB_OUTDIR', '~/agfr_web_runs')


def main(argv=None):
    from . import APP, APP_LONG, __version__
    p = argparse.ArgumentParser(
        prog='agfr-web',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='%s %s: %s\n\n'
                    'Prepares receptor and ligand PDBQT files from a PDB entry in a web\n'
                    'browser, the way AGFRgui\'s PDB panel does, with ADFR\'s own code.'
                    % (APP, __version__, APP_LONG),
        epilog='examples:\n'
               '  %(prog)s serve                   start %(app)s on port 8766 and print its URL\n'
               '  %(prog)s serve --port 9000       another fixed port (the browser remembers\n'
               '                                   folders and theme per port)\n'
               '  %(prog)s check                   "ready", or what %(app)s is missing\n'
               '  %(prog)s examples                where ADFR keeps its example PDB files\n'
               '  %(prog)s test -k parity          run part of the test suite\n\n'
               '%(app)s runs ADFR from the ADCP distribution (the adcpsuite env); the\n'
               'wrapper bin/agfr-web sets that env up. README.md has a first session.'
               % {'prog': 'agfr-web', 'app': APP})
    p.add_argument('--version', action='version', version='%s %s' % (APP, __version__))
    sub = p.add_subparsers(dest='cmd')
    s = sub.add_parser('serve', help='start %s in a browser page (prints its URL)' % APP)
    s.add_argument('--host', default='127.0.0.1',
                   help='address to bind (default 127.0.0.1: this machine only)')
    s.add_argument('--port', type=int, default=8766, help='0 picks a free port')
    s.add_argument('--outdir', default=DEFAULT_OUTDIR,
                   help='where every session keeps its input, outputs and '
                        'session.json (default %(default)s); keep it on NFS, '
                        'never /tmp')
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
    if a.cmd == 'check':
        miss = env.check()
        print('%s is ready' % APP if not miss else '%s cannot run: missing %s' % (APP, '; '.join(miss)))
        return 1 if miss else 0
    if a.cmd != 'serve':
        p.print_help()
        return 2
    if os.path.realpath(os.path.expanduser(a.outdir)).startswith('/tmp'):
        p.error('--outdir must not be under /tmp; it is wiped')
    from .server import serve
    serve(a.host, a.port, a.outdir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
