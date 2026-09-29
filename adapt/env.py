# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The ADFR installation this tool drives, and what it records about it.

The computation is ADFR's own (the AGFRgui code paths), run inside the ADCP
distribution's Python (the `adcpsuite` env, Python 3.7). This module must be
imported before any ADFR module is used, because ADFR finds `reduce` through
$ADS_ROOT (ADFR/utils/MakeGrids.py:findBinary), which the agfr/agfrgui launcher
scripts export and a bare `python` does not.
"""
import glob
import hashlib
import os
import subprocess
import sys

# the launcher scripts compute ADS_ROOT as <bin>/.. of the interpreter
ADS_ROOT = os.path.realpath(os.path.join(os.path.dirname(sys.executable), '..'))
os.environ.setdefault('ADS_ROOT', ADS_ROOT)
os.environ.setdefault('PYTHONNOUSERSITE', 'disable')

# the ADFR source files whose md5 identifies the computation that ran
_TRACKED = ('utils/runAGFR.py', 'recFromPDB.py', 'utils/ADprep.py',
            'GUI/pdbAnalyser.py', '__init__.py')


def _md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def check():
    """Names of what is missing for the PDB panel to work; empty when ready."""
    missing = []
    try:
        import ADFR  # noqa: F401
    except ImportError:
        missing.append('ADFR (run with the adcpsuite python)')
    reduce_bin = os.path.join(os.environ['ADS_ROOT'], 'bin', 'reduce')
    if not os.access(reduce_bin, os.X_OK):
        missing.append('reduce (%s)' % reduce_bin)
    if not os.path.isfile(reduce_bin + '_wwPDB_het_dict.txt'):
        missing.append('reduce het dictionary (%s_wwPDB_het_dict.txt)' % reduce_bin)
    return missing


def provenance():
    """What produced a file: interpreter, ADFR package version, and the md5 of
    the ADFR sources that do the work, so any output can be tied to its code."""
    import ADFR
    pkg = os.path.dirname(os.path.abspath(ADFR.__file__))
    dist = sorted(glob.glob(os.path.join(os.path.dirname(pkg), 'adfr-*.dist-info')))
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        commit = subprocess.run(['git', '-C', here, 'rev-parse', 'HEAD'],
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                universal_newlines=True).stdout.strip() or None
        dirty = subprocess.run(['git', '-C', here, 'status', '--porcelain'],
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                               universal_newlines=True).stdout.strip() != ''
    except OSError:
        commit, dirty = None, None
    return {
        'python': sys.version.split()[0],
        'executable': sys.executable,
        'ADS_ROOT': os.environ.get('ADS_ROOT'),
        'adfr_package': pkg,
        'adfr_dist': os.path.basename(dist[-1])[:-len('.dist-info')] if dist else None,
        'adfr_md5': {f: _md5(os.path.join(pkg, f)) for f in _TRACKED
                     if os.path.isfile(os.path.join(pkg, f))},
        'adapt_commit': commit,
        'adapt_dirty': dirty,
    }
