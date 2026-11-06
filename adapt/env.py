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

# the ADFR and MolKit2 source files whose md5 identifies the computation that ran
_TRACKED = ('utils/runAGFR.py', 'recFromPDB.py', 'utils/ADprep.py',
            'GUI/pdbAnalyser.py', '__init__.py')
_TRACKED_MOLKIT2 = ('openBabelInterface.py', 'AARotamer.py', 'molecule.py')


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
    else:
        from ADFR import recFromPDB
        if not (hasattr(recFromPDB.PDBProfiler, 'chooseAltloc')
                and hasattr(recFromPDB.PDBProfiler, 'completeSidechain')
                and hasattr(recFromPDB, 'ccdIonCharge')):
            missing.append('ADFR 0.0.23 or later (this one, %s, predates the fixes '
                           'ADapt relies on)' % _dist(os.path.dirname(os.path.abspath(
                               ADFR.__file__)), 'adfr'))
    reduce_bin = os.path.join(os.environ['ADS_ROOT'], 'bin', 'reduce')
    if not os.access(reduce_bin, os.X_OK):
        missing.append('reduce (%s)' % reduce_bin)
    if not os.path.isfile(reduce_bin + '_wwPDB_het_dict.txt'):
        missing.append('reduce het dictionary (%s_wwPDB_het_dict.txt)' % reduce_bin)
    return missing


def _git(path, *args):
    """Output of `git -C path args`, or None (not a checkout, no git)."""
    try:
        out = subprocess.run(['git', '--no-optional-locks', '-C', path] + list(args),
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             universal_newlines=True)
    except OSError:
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def _dist(pkg, name):
    """The installed distribution of package folder `pkg` ('adfr-0.0.22'), or,
    for a source checkout used in development (ADAPT_DEV_PATH), its git
    revision ('adfr-v0.0.22-3-g1a2b3c4-dirty (source)')."""
    dist = sorted(glob.glob(os.path.join(os.path.dirname(pkg), name + '-*.dist-info')))
    if dist:
        return os.path.basename(dist[-1])[:-len('.dist-info')]
    rev = _git(os.path.realpath(pkg), 'describe', '--tags', '--always', '--dirty')
    return '%s-%s (source)' % (name, rev) if rev else None


def _sums(pkg, files):
    return {f: _md5(os.path.join(pkg, f)) for f in files
            if os.path.isfile(os.path.join(pkg, f))}


def _snapshot_commit():
    """The commit `git archive` stamped into adapt/_build.py, or None in a git
    checkout (where the placeholders are not expanded)."""
    from . import _build
    c = getattr(_build, 'COMMIT', '')
    return None if not c or c.startswith('$Format') else c


def provenance():
    """What produced a file: interpreter, ADFR and MolKit2 versions, and the md5
    of their sources that do the work, so any output can be tied to its code."""
    import ADFR
    import MolKit2
    pkg = os.path.dirname(os.path.abspath(ADFR.__file__))
    mk2 = os.path.dirname(os.path.abspath(MolKit2.__file__))
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # git only speaks for ADapt when ADapt's folder is the top of a checkout:
    # a snapshot unpacked inside some other repository must not report that
    # repository's commit
    top = _git(here, 'rev-parse', '--show-toplevel')
    checkout = top is not None and os.path.realpath(top) == os.path.realpath(here)
    status = _git(here, 'status', '--porcelain') if checkout else None
    commit = _git(here, 'rev-parse', 'HEAD') if checkout else None
    snapshot = None if commit else _snapshot_commit()
    return {
        'python': sys.version.split()[0],
        'executable': sys.executable,
        'ADS_ROOT': os.environ.get('ADS_ROOT'),
        'adfr_package': pkg,
        'adfr_dist': _dist(pkg, 'adfr'),
        'adfr_md5': _sums(pkg, _TRACKED),
        'molkit2_package': mk2,
        'molkit2_dist': _dist(mk2, 'molkit2'),
        'molkit2_md5': _sums(mk2, _TRACKED_MOLKIT2),
        'adapt_commit': commit or snapshot,
        # a snapshot (git archive) cannot say whether it was edited after unpacking
        'adapt_source': 'git' if commit else ('snapshot' if snapshot else None),
        'adapt_dirty': None if status is None else status != '',
    }
