"""python -m agfr_web serve [--host H] [--port P] [--outdir DIR] | check | examples | test [-k NAME]"""
import argparse
import os
import sys

DEFAULT_OUTDIR = os.environ.get('AGFR_WEB_OUTDIR', '~/agfr_web_runs')


def main(argv=None):
    p = argparse.ArgumentParser(prog='agfr-web', description=__doc__)
    sub = p.add_subparsers(dest='cmd')
    s = sub.add_parser('serve', help='start the page (prints its URL)')
    s.add_argument('--host', default='127.0.0.1',
                   help='address to bind (default 127.0.0.1: this machine only)')
    s.add_argument('--port', type=int, default=8766, help='0 picks a free port')
    s.add_argument('--outdir', default=DEFAULT_OUTDIR,
                   help='where every session keeps its input, outputs and '
                        'session.json (default %(default)s); keep it on NFS, '
                        'never /tmp')
    sub.add_parser('check', help='report what is missing')
    sub.add_parser('examples', help="list the example PDB files that ship with ADFR")
    t = sub.add_parser('test', help='run the test suite (tests/run_tests.py)')
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
        print('ready' if not miss else 'missing: ' + '; '.join(miss))
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
