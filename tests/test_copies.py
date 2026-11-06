# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Symmetry copies follow their original chain (adapt/copies.py): receptor
checkboxes, alternate locations and ion charges reach every copy, a ligand
is taken from the original only, atoms on a symmetry axis are written once,
and a chain first made by a later REMARK 350 block is in the tree.

Copies used to keep the defaults: with a residue unchecked, its copies were
still written (1HHP A:ASP25)."""
import os

from common import (CHECKED, DATA, LIG, REC, UNCHECKED, check, find, oracle_receptor, outdir,
                    quiet, read, session, set_state)
from adapt import env
from adapt.session import Session

HERE = os.path.dirname(os.path.abspath(__file__))
_n = [0]


def _session_of(path, name, label):
    _n[0] += 1
    with open(path, 'rb') as fh:
        data = fh.read()
    with quiet():
        return Session('c%015d' % _n[0], data, name, outdir(label), env.provenance())


def _with_assemblies(name, blocks):
    """6ggd.pdb whose REMARK 350 has one assembly made of `blocks`:
    [(chains, [x shift of each copy])]."""
    lines = [l for l in open(os.path.join(DATA, '6ggd.pdb')) if not l.startswith('REMARK 350')]
    at = next(i for i, l in enumerate(lines) if l.startswith('SEQRES'))
    new = ['REMARK 350 BIOMOLECULE: 1\n']
    for chains, shifts in blocks:
        new.append('REMARK 350 APPLY THE FOLLOWING TO CHAINS: %s\n' % ', '.join(chains))
        for k, dx in enumerate(shifts):
            for r in range(3):
                row = ''.join('%10.6f' % (1.0 if c == r else 0.0) for c in range(3))
                new.append('REMARK 350   BIOMT%d %3d%s%15.5f\n' % (r + 1, k + 1, row, dx if r == 0 else 0.0))
    lines[at:at] = new
    path = os.path.join(outdir('copies_input'), name)
    with open(path, 'w') as fh:
        fh.writelines(lines)
    return path


def _two_copies(label):
    """A session of 6GGD whose Biomol 1 is chain A twice (in place, +80 A)."""
    s = _session_of(_with_assemblies(label + '.pdb', [(['A'], [0.0, 80.0])]), '6ggd.pdb', label)
    s.switch('Biomol 1')
    return s


def _atoms(s, kind):
    d = s.save(kind)
    folder = os.path.join(s.folder, 'save_%03d_%s' % (d['save'], kind[:3]))
    return [l for f in d['files'] if f['name'].endswith('.pdbqt')
            for l in read(os.path.join(folder, f['name'])).decode().splitlines()
            if l.startswith(('ATOM', 'HETATM'))], d


def test_unchecked_residue_left_out_of_every_copy():
    """1H6W Biomol 1 (chains A, B; three copies) with A:LYS:... unchecked:
    the receptor equals agfr's with --protToRemove for that residue, which
    removes it from every copy."""
    s = session('1h6w.pdb', 'copies_uncheck')
    s.switch('Biomol 1')
    res = [n for n in find(s.tree(), kind='residue', chain='A', parent='protein')][10]
    set_state(s, res, REC, UNCHECKED)
    chid, resname, resnum, icode = res.key.split(':')
    d = s.save('receptor')
    web = os.path.join(s.folder, 'save_%03d_rec' % d['save'], d['files'][0]['name'])
    ref = oracle_receptor('1h6w.pdb', outdir('oracle_copies_uncheck'), bioMol='1',
                          protToRemove='chid %s resnum %s' % (chid, resnum))
    check(read(web) == read(ref), 'receptor differs from agfr --protToRemove')


def test_altloc_choice_in_every_copy():
    """A:SER:166 defaults to location B (0.58); location A chosen: both copies
    have location A's CA (115.647, in place and moved 80 A)."""
    s = _two_copies('copies_altloc')
    node = find(s.tree(), key='A:SER:166:_')[0]
    s.altloc(node.id, node.altloc['alts'].index('A'))
    atoms, _ = _atoms(s, 'receptor')
    ca = sorted(float(l[30:38]) for l in atoms if l[17:20] == 'SER' and l[22:26] == ' 166'
                and l[12:16].strip() == 'CA')
    check(ca == [115.647, 195.647], ca)


def test_ligand_taken_from_the_original_only():
    """A:SER:96 checked as ligand: one ligand file, and SER 96 in no copy of
    the receptor (the copy kept it before)."""
    s = _two_copies('copies_ligand')
    node = find(s.tree(), key='A:SER:96:_')[0]
    set_state(s, node, LIG, CHECKED)
    rec, _ = _atoms(s, 'receptor')
    check(not [l for l in rec if l[17:20] == 'SER' and l[22:26] == '  96'], 'SER 96 in a copy')
    lig, d = _atoms(s, 'ligand')
    check(len([f for f in d['files'] if f['name'].endswith('.pdbqt')]) == 1, d['files'])


def test_ion_in_every_copy():
    """A:ZN:401 checked into the receptor: both copies' zinc, at +2.00."""
    s = _two_copies('copies_zn')
    set_state(s, find(s.tree(), key='A:ZN:401:_')[0], REC, CHECKED)
    atoms, _ = _atoms(s, 'receptor')
    zn = [float(l[70:76]) for l in atoms if l[17:20] == ' ZN']
    check(zn == [2.0, 2.0], zn)


def test_axis_ions_written_once():
    """4INS Biomol 3, the insulin hexamer: its two zinc ions lie on the 3-fold
    axis, so the copies' zincs fall on the originals; checked into the
    receptor, each is written once, not three times."""
    s = _session_of(os.path.join(HERE, 'data', '4ins.pdb.gz'), '4ins.pdb.gz', 'copies_4ins')
    s.switch('Biomol 3')
    for key in ('B:ZN:101:_', 'D:ZN:101:_'):
        set_state(s, find(s.tree(), key=key)[0], REC, CHECKED)
    atoms, _ = _atoms(s, 'receptor')
    zn = [l for l in atoms if l[17:20] == ' ZN']
    check(len(zn) == 2, len(zn))


def test_chain_of_a_later_block_in_the_tree():
    """Biomol 1 = chain A twice, then chain B once in a second REMARK 350
    block: B's only copy is not in segment A, and the tree, which listed
    segment A, had no chain B."""
    path = _with_assemblies('6ggd_blocks.pdb', [(['A'], [0.0, 80.0]), (['B'], [100.0])])
    s = _session_of(path, '6ggd.pdb', 'copies_blocks')
    s.switch('Biomol 1')
    check([r.label for r in s.tree().roots] == ['A', 'B'], [r.label for r in s.tree().roots])
