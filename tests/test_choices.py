# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Choices that reach the files: the rotamer chosen is the one built, an
unsupported ion can be replaced (and the replacement undone), and residues
left out for an atom type AD4.1 lacks are in the save record."""
import json
import os

import numpy

from common import CHECKED, REC, check, find, outdir, quiet, read, session
from adapt import env
from adapt.session import Session

HERE = os.path.dirname(os.path.abspath(__file__))
_n = [0]


def _session_of(name, label):
    _n[0] += 1
    with open(os.path.join(HERE, 'data', name), 'rb') as fh:
        data = fh.read()
    with quiet():
        return Session('h%015d' % _n[0], data, name, outdir(label), env.provenance())


def _saved(s, kind):
    d = s.save(kind)
    folder = os.path.join(s.folder, 'save_%03d_%s' % (d['save'], kind[:3]))
    text = read(os.path.join(folder, d['files'][0]['name'])).decode()
    rec = json.load(open(os.path.join(s.folder, 'session.json')))
    return text, rec['saves'][-1], rec


def _groups(node):
    out = []
    while node.parent is not None:
        node = node.parent
        out.append(node.label)
    return out


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, (p2 - p1) / numpy.linalg.norm(p2 - p1), p3 - p2
    v, w = b0 - numpy.dot(b0, b1) * b1, b2 - numpy.dot(b2, b1) * b1
    return numpy.degrees(numpy.arctan2(numpy.dot(numpy.cross(b1, v), w), numpy.dot(v, w)))


def test_chosen_rotamer_is_built():
    """6GGD A:LYS:120: rotamer 41 chosen in the params column (browsing 40
    then 41: one event). It was a preview: the best-scoring rotamer was built."""
    from MolKit2.rotamerLib import RotamerLib
    s = session('6ggd.pdb', 'choice_rotamer')
    s.switch('ASU')
    node = find(s.tree(), key='A:LYS:120:_')[0]
    s.rotamer(node.id, 39)
    s.rotamer(node.id, 40)
    text, entry, rec = _saved(s, 'receptor')
    check('REMARK 850 AGFR rotamer A:LYS:120:_ LYS 41 of 81' in text)
    check(entry['rotamers'] == {'A:LYS:120:_': 41}, entry['rotamers'])
    check([e['rotamer'] for e in rec['events'] if e['what'] == 'rotamer'] == [41], rec['events'])
    xyz = dict((l[12:16].strip(), numpy.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
               for l in text.splitlines() if l.startswith('ATOM') and l[17:26] == 'LYS A 120')
    chis = [_dihedral(*[xyz[a] for a in d[0]]) for d in RotamerLib().getAngleDef('LYS')]
    want = RotamerLib().getAngles('LYS')[40]
    check(all(abs((x - y + 180) % 360 - 180) < 1.0 for x, y in zip(chis, want)), (chis, want))
    s.rotamer(node.id, -1)                    # back to auto
    check(s.profile._rotamerChoice == {}, s.profile._rotamerChoice)


def test_copper_replaced_and_restored():
    """1PLC's copper, under 'unknown atom types', replaced by zinc from its
    menu: it moves to 'manual charges', checked, and the receptor has a zinc
    at +2.00. Undone: the copper is out again."""
    s = _session_of('1plc.pdb.gz', 'choice_cu')
    cu = find(s.tree(), key='A:CU:100:_')[0]
    check(cu.menu == 'substitute' and 'unknown atom types' in _groups(cu), (cu.menu, _groups(cu)))
    s.substitute(cu.id, 'ZN')
    zn = find(s.tree(), key='A:CU:100:_')[0]
    check(zn.parent.label == 'manual charges' and zn.check_state(REC) == CHECKED,
          (zn.parent.label, zn.check_state(REC)))
    text, entry, rec = _saved(s, 'receptor')
    zn_lines = [l for l in text.splitlines() if l.startswith(('ATOM', 'HETATM')) and l[17:20] == ' ZN']
    check(len(zn_lines) == 1 and float(zn_lines[0][70:76]) == 2.0, zn_lines)
    check('REMARK 850 AGFR substituted A:CU:100:_ CU -> ZN charge 2.00' in text)
    check(entry['substitutions'] == {'A:CU:100:_': ['ZN', 2.0]}, entry['substitutions'])
    s.substitute(zn.id, None)
    back = find(s.tree(), key='A:CU:100:_')[0]
    check('unknown atom types' in _groups(back), _groups(back))
    text, entry, rec = _saved(s, 'receptor')
    check(not [l for l in text.splitlines() if l[17:20] in (' ZN', ' CU')], 'no copper, no zinc')


def test_left_out_in_the_save_record():
    """A residue with an atom AD4.1 cannot type (1DDY's B12, cobalt) put in the
    receptor is dropped when the file is written; the save says so."""
    s = _session_of('1ddy.pdb.gz', 'choice_b12')
    s.switch('ASU')
    s.profile._currentBM._ag.select('chid A resname B12').setFlags('receptor', True)
    text, entry, rec = _saved(s, 'receptor')
    check(entry['left_out'] == {'A:B12:101:_': ['CO']}, entry['left_out'])
    check('REMARK 850 AGFR left out A:B12:101:_ (atoms of type CO)' in text)
