# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Missing residues as ADapt shows them: tree nodes and 3D markers.

ADFR finds them (PDBProfiler._gapsPerChain); its own tests check them against
the files' REMARK 465 records (ADFR/tests/testGaps.py)."""
import numpy

from common import check, find, session


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
