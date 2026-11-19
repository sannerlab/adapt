# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The receptor and ligand PDBQT files ADapt writes can be read by AutoDock Vina.

Vina 1.2 (parse_pdbqt_atom_string) reads an atom's charge from columns 69-76, skipping
leading blanks but refusing a trailing blank, and its AD type from columns 78-79. ADFR up
to 0.0.23 wrote the segment name in 67-70 and the charge as %5.2f in 71-75, so Vina 1.2.7
refused every file ('Charge "A-0.13 " is not valid'); ADFR 0.0.24 writes AutoDock's columns."""
import os
import re

from common import CHECKED, LIG, check, session, set_state

# the atom types Vina 1.2.7 accepts (src/lib/atom_constants.h: atom_kind_data, non_ad_metal_names)
VINA_TYPES = set('C A N O P S H F I NA OA SA HD Mg Mn Zn Ca Fe Cl Br Si At G0 G1 G2 G3 '
                 'CG0 CG1 CG2 CG3 W Cu Na K Hg Co U Cd Ni'.split())


def _save(s, kind):
    d = s.save(kind)
    folder = os.path.join(s.folder, 'save_%03d_%s' % (d['save'], kind[:3]))
    return [os.path.join(folder, f['name']) for f in d['files'] if f['name'].endswith('.pdbqt')]


def _vina_refusals(path):
    with open(path) as fh:
        atoms = [l.rstrip('\n') for l in fh if l.startswith(('ATOM', 'HETATM'))]
    check(atoms, 'no atom lines in %s' % path)
    bad = []
    for l in atoms:
        charge = l[68:76].lstrip()
        if len(l) < 78:
            bad.append(('short', l))
        elif l[66:70] != '    ':
            bad.append(('columns 67-70 not blank', l))
        elif not re.match(r'^-?\d+\.\d{3}$', charge):
            bad.append(('charge %r' % charge, l))
        elif l[77:79].strip() not in VINA_TYPES:
            bad.append(('type %r' % l[77:79], l))
    return atoms, bad


def test_receptor_and_ligand_readable_by_vina():
    """6bhg: the default receptor (segment A in ADFR's atoms) and chain B as ligand."""
    s = session('6bhg.pdb', 'vina_6bhg')
    t = s.tree()
    root = [r for r in t.roots if r.label == 'B'][0]
    targets = []
    for n in root.children:
        if n.kind == 'category' and n.label != 'water':
            targets.append(n)
        if n.kind == 'group' and n.label == 'altloc':
            targets.extend(c for c in n.children if c.label != 'water')
    for n in targets:
        set_state(s, n, LIG, CHECKED)
    for path in _save(s, 'ligand') + _save(s, 'receptor'):
        atoms, bad = _vina_refusals(path)
        check(not bad, '%s: %d of %d atom lines, first: %s'
              % (os.path.basename(path), len(bad), len(atoms), bad[:1]))
