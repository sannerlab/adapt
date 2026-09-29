# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Missing residues: ADapt's fixed PDBProfiler._gaps, checked against the PDB
files' own REMARK 465 records (the residues the depositors list as missing)."""
import collections
import os

import numpy

from common import DATA, PDBS, check, find, quiet, session
from adapt import gaps


def remark465(path):
    """{chain: [(resname, number)]} from REMARK 465 (model number optional)."""
    out = collections.defaultdict(list)
    started = False
    for l in open(path):
        if not l.startswith('REMARK 465'):
            continue
        if 'RES C SSSEQI' in l:
            started = True
            continue
        f = l[10:].split()
        if started and len(f) in (3, 4):
            num = int(''.join(ch for ch in f[-1] if ch.isdigit() or ch == '-'))
            out[f[-2]].append((f[-3], num))
    return out


def _profile(pdb):
    from ADFR.utils.runAGFR import runAGFR
    from adapt.session import gui_options
    with quiet():
        a = runAGFR()
        a.proccessPreLoadCmdLineOptions(gui_options())
        a.loadReceptor(os.path.join(DATA, pdb), processReceptor=False)
    return a.PDBprofile


def _count(v):
    mrn, gl, mrc = v
    return (len(mrn[1]) if mrn else 0) + sum(len(g) - 2 for g in gl) + (len(mrc[1]) if mrc else 0)


def _placed(v, num):
    """Is residue number `num` inside one of the reported missing stretches?"""
    mrn, gl, mrc = v
    if mrn and num < mrn[0][1]:
        return True
    if mrc and num > mrc[0][1]:
        return True
    return any(g[0][1] < num < g[-1][1] for g in gl)


def test_adfr_original_reports_no_gap():
    """The defect: ADFR's own _gaps (adfr-0.0.22) reports nothing, even for
    files with REMARK 465 records."""
    gaps.uninstall()
    try:
        for pdb in ('4a1x.pdb', '6ggd.pdb', '1h6w.pdb'):
            check(remark465(os.path.join(DATA, pdb)), pdb + ' lists missing residues')
            check(_profile(pdb)._gapsPerChain == {}, pdb)
    finally:
        gaps.install()


def test_gaps_match_remark465():
    """Every chain of the ten ADFR test PDBs: the same number of missing
    residues as REMARK 465, each at the right place; the one exception is the
    stretch ADFR skips on purpose (1jff A 35-60, whose ends are 2.68 A apart)."""
    for pdb in PDBS:
        p = _profile(pdb)
        want = remark465(os.path.join(DATA, pdb))
        for chid in set(want) | set(p._gapsPerChain):
            v = p._gapsPerChain.get(chid, ['', [], ''])
            listed = want.get(chid, [])
            if (pdb, chid) == ('1jff.pdb', 'A'):
                skipped = [n for _, n in listed if not _placed(v, n)]
                check(skipped == list(range(35, 61)), skipped)
                ag = p._receptor._ag
                c = ag.select('chid A resnum 34 name C').getCoords()[0]
                n = ag.select('chid A resnum 61 name N').getCoords()[0]
                check(numpy.linalg.norm(c - n) < 3.5, 'ADFR: no physical gap below 3.5 A')
                listed = [(r, k) for r, k in listed if k not in skipped]
            check(_count(v) == len(listed), '%s %s: %d reported, REMARK 465 lists %d'
                  % (pdb, chid, _count(v), len(listed)))
            for resname, num in listed:
                check(_placed(v, num), '%s %s %s%d not inside a reported stretch'
                      % (pdb, chid, resname, num))


def test_modified_residues_are_not_missing():
    """FIX 4: MLZ (5fq2 B505) and ALY (6bhg B14) are HETATM residues in the
    chain, and NH2 (1a3r P171) a cap: present, so not reported missing."""
    for pdb, chid, flank in (('5fq2.pdb', 'B', 504), ('6bhg.pdb', 'B', 13)):
        v = _profile(pdb)._gapsPerChain.get(chid, ['', [], ''])
        check(not [g for g in v[1] if g[0][1] == flank], '%s %s: no gap after %d' % (pdb, chid, flank))
    check('P' not in _profile('1a3r.pdb')._gapsPerChain, '1a3r P: nothing missing')


def test_his_tag_terminus():
    """FIX 5: 5hht B lacks its last His (669) of HHHHHH; it is a missing C-terminus."""
    v = _profile('5hht.pdb')._gapsPerChain['B']
    check(v[2] and v[2][0][1] == 668 and v[2][1] == 'H', v[2])


def test_tree_and_3d_markers_4a1x():
    s = session('4a1x.pdb', 'gaps_4a1x')
    s.switch('ASU')
    labels = dict((g.parent.label, [c.label for c in g.children])
                  for g in find(s.tree(), kind='gaps'))
    check(labels == {'A': ['gap: SER-3--GKGG--PRO2', 'C term: GSHHHHHH'],
                     'B': ['gap: SER-3--GK--GLY0', 'C term: SHHHHHH'],
                     'C': ['N term: GEL', 'C term: D'],
                     'D': ['N term: GELGR', 'C term: D']}, labels)
    geo = s.structure()['gaps']
    check(len(geo) == 8, len(geo))
    for g in geo:
        p, q = numpy.array(g['line'][0]), numpy.array(g['line'][1])
        if g['kind'] == 'gap':
            check(numpy.allclose(g['sphere'], (p + q) / 2), 'sphere at the middle')
            check(numpy.linalg.norm(p - q) > 3.5, 'a physical gap')
        else:
            check(abs(numpy.linalg.norm(p - q) - 3.0) < 1e-6, 'terminus sphere 3 A out')
    s.switch('Biomol 1')                       # only chains A and C
    check(sorted(g.parent.label for g in find(s.tree(), kind='gaps')) == ['A', 'C'])
    check(len(s.structure()['gaps']) == 4)
