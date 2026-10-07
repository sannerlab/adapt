# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Shared helpers: the ADFR test PDBs, a per-run output folder, the headless
ADFR oracle, and tree-walking helpers that click the way a user would."""
import contextlib
import datetime
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from adapt import env  # noqa: E402  (sets ADS_ROOT)
from adapt.session import LOCK, Session  # noqa: E402
from adapt.tree import CHECKED, LIG, PARTIAL, REC, UNCHECKED  # noqa: E402,F401

import ADFR  # noqa: E402

DATA = os.path.join(os.path.dirname(ADFR.__file__), 'tests', 'Data')
PDBS = sorted(f for f in os.listdir(DATA) if f.endswith('.pdb'))

# outputs of one test run, inside the repo (git-ignored), never /tmp
OUT = os.path.join(HERE, '_out', datetime.datetime.now().strftime('%Y%m%d-%H%M%S'))

_n = [0]


def outdir(name):
    d = os.path.join(OUT, name)
    os.makedirs(d, exist_ok=True)
    return d


def session(pdb, name=None):
    _n[0] += 1
    path = os.path.join(DATA, pdb)
    with open(path, 'rb') as fh:
        data = fh.read()
    return Session('t%015d' % _n[0], data, pdb, outdir(name or 'sessions'),
                   env.provenance())


@contextlib.contextmanager
def quiet():
    buf = io.StringIO()
    with LOCK, contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        yield buf


def oracle_receptor(pdb, folder, **extra):
    """ADFR's own headless receptor path (what ADFR/tests/testAGFR.py runs):
    runAGFR()(receptorFile=..., toPdbqt=True, ...) -> <pdbid>_rec.pdbqt."""
    from ADFR.utils.runAGFR import runAGFR
    from ADFR.tests.defaultkw import kw_base
    kw = kw_base.copy()
    kw.update({'receptorFile': os.path.join(DATA, pdb), 'toPdbqt': True})
    kw.update(extra)
    from adapt import gaps
    cwd = os.getcwd()
    os.chdir(folder)
    gaps.uninstall()        # ADFR exactly as installed, including its _gaps
    try:
        with quiet():
            runAGFR()(**kw)
    finally:
        gaps.install()
        os.chdir(cwd)
    return os.path.join(folder, os.path.splitext(pdb)[0] + '_rec.pdbqt')


def read(path):
    with open(path, 'rb') as fh:
        return fh.read()


def walk(node):
    yield node
    for c in node.children:
        for x in walk(c):
            yield x


def all_nodes(tree):
    for r in tree.roots:
        for n in walk(r):
            yield n


def find(tree, label=None, kind=None, chain=None, key=None, parent=None):
    out = []
    for r in tree.roots:
        if chain is not None and r.label != chain:
            continue
        for n in walk(r):
            if label is not None and n.label.split(' ')[0] != label and n.label != label:
                continue
            if kind is not None and n.kind != kind:
                continue
            if key is not None and n.key != key:
                continue
            if parent is not None and (n.parent is None or n.parent.label != parent):
                continue
            out.append(n)
    return out


def set_state(s, node, col, want):
    """Click a node's checkbox until it shows `want` (partial -> checked)."""
    for _ in range(3):
        if node.check_state(col) == want:
            return
        s.click(node.id, 'rec' if col == REC else 'lig')
    raise AssertionError('%s never reached %s' % (node.label, want))


def check(cond, msg=''):
    if not cond:
        raise AssertionError(msg)


BIOMOL = b'REMARK 850 AGFR biomolecule '


def without_biomolecule(data):
    """A PDBQT file as ADFR writes it: ADapt's bytes minus the one line it adds
    (REMARK 850 AGFR biomolecule ...), which must be there exactly once, right
    after the "source" line. ADFR's own files pass through unchanged."""
    lines = data.split(b'\n')
    at = [i for i, l in enumerate(lines) if l.startswith(BIOMOL)]
    if not at:
        return data
    check(len(at) == 1, 'one biomolecule line, found %d' % len(at))
    check(lines[at[0] - 1].startswith((b'REMARK 850 AGFR source', b'REMARK 850 source')),
          'biomolecule line follows source')
    return b'\n'.join(lines[:at[0]] + lines[at[0] + 1:])
