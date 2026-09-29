"""python -m agfr_web serve [--host H] [--port P] [--outdir DIR]"""
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
    a = p.parse_args(argv)
    from . import env
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
