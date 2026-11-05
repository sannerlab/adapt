# Author: Michel Sanner, PhD. Copyright (c) 2026. Derived from ADFR / AGFRgui code; see THIRD_PARTY_NOTICES.md.
"""Missing residues in the 3D view.

ADFR finds the missing residues of each chain (PDBProfiler._gapsPerChain:
missing N-terminus, gaps, missing C-terminus) by aligning its SEQRES sequence
with the residues present; its tests (ADFR/tests/testGaps.py) check them
against the files' REMARK 465 records. geometry() makes what the view draws
for them.
"""
import numpy as np


# ---- what the 3D view draws for them -------------------------------------
BOND = 1.9     # A: neighbours of a terminal atom, found by distance (no bond table needed)


def _short(seq):
    """gridGUI's label rule: more than 10 residues -> first 3 ... last 3."""
    return seq if len(seq) <= 10 else seq[:3] + '...' + seq[-3:]


def _res(ag, chid, num, icode):
    sel = 'segment A and chid %s resnum `%d`' % (chid if chid != ' ' else '_', num)
    if icode:
        sel += ' icode %s' % icode
    return ag.select(sel)


def geometry(gaps_per_chain, mol):
    """Port of GridGUI.addGeomsForGaps (gridGUI.py): for each gap a dotted line
    from the C of the residue before it to the N of the residue after it, with
    a sphere at its middle labelled with the missing residues; for a missing
    terminus a sphere 3 A beyond the terminal N (or C) atom, opposite its bonds,
    joined to it by a dotted line. Returns [{'line': [p, q], 'sphere': c,
    'label': str, 'kind': 'gap'|'N'|'C'}], coordinates as lists of floats."""
    import numpy as np
    ag = mol._ag
    out = []
    for chid, chgaps in gaps_per_chain.items():
        for gap in chgaps[1]:
            r1 = _res(ag, chid, *gap[0][1:])
            r2 = _res(ag, chid, *gap[-1][1:])
            if not r1 or not r2:
                continue
            c, n = r1.select('name C'), r2.select('name N')
            if c is None or n is None:
                continue
            p, q = c.getCoords()[0], n.getCoords()[0]
            out.append({'kind': 'gap', 'line': [p.tolist(), q.tolist()],
                        'sphere': ((p + q) / 2.).tolist(),
                        'label': _short(''.join(gap[1:-1]))})
        for term, atname in ((chgaps[0], 'N'), (chgaps[-1], 'C')):
            if not len(term):
                continue
            res = _res(ag, chid, *term[0][1:])
            if not res:
                continue
            atoms = res.select('name %s' % atname) or res.select('name CA')
            if atoms is None:
                continue
            at = atoms.getCoords()[0]
            near = ag.select('segment A and chid %s and within %.2f of pt' % (
                chid if chid != ' ' else '_', BOND), pt=at)
            vecs = [x - at for x in (near.getCoords() if near is not None else [])
                    if 0.1 < np.linalg.norm(x - at)]
            if not vecs:
                continue
            v = np.add.reduce(vecs)
            if np.linalg.norm(v) < 1e-6:
                continue
            sp = at - 3 * v / np.linalg.norm(v)
            out.append({'kind': atname, 'line': [sp.tolist(), at.tolist()],
                        'sphere': sp.tolist(), 'label': _short(term[1])})
    return out
