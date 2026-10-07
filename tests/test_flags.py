# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""What the checkboxes, altloc radios, charges and menus do to the atoms."""
import os

from common import (CHECKED, LIG, PARTIAL, REC, UNCHECKED, check, find,
                    read, session, set_state)


def _flags(sel):
    return (bool(sel.getFlags('receptor').any()), bool(sel.getFlags('receptor').all()),
            bool(sel.getFlags('ligand').any()), bool(sel.getFlags('ligand').all()))


def test_receptor_and_ligand_are_exclusive():
    s = session('1a3r.pdb', 'excl')
    t = s.tree()
    res = find(t, kind='residue', parent='protein')[5]
    cat = res.parent
    check(res.check_state(REC) == CHECKED and cat.check_state(REC) == CHECKED)
    s.click(res.id, 'lig')
    check(res.check_state(LIG) == CHECKED and res.check_state(REC) == UNCHECKED)
    check(_flags(res.sel) == (False, False, True, True), _flags(res.sel))
    check(cat.check_state(REC) == PARTIAL and cat.check_state(LIG) == PARTIAL)
    s.click(res.id, 'rec')
    check(res.check_state(REC) == CHECKED and res.check_state(LIG) == UNCHECKED)
    check(_flags(res.sel) == (True, True, False, False), _flags(res.sel))
    check(cat.check_state(REC) == CHECKED)


def test_category_click_applies_to_every_residue():
    s = session('1a3r.pdb', 'cat')
    cat = find(s.tree(), label='protein', kind='category')[0]
    s.click(cat.id, 'rec')
    check(cat.check_state(REC) == UNCHECKED)
    check(all(c.check_state(REC) == UNCHECKED for c in cat.children))
    check(cat.sel.select('receptor') is None)
    check(s.buttons()['receptor'])      # the other chains are still receptor
    s.click(cat.id, 'rec')
    check(all(c.check_state(REC) == CHECKED for c in cat.children))
    check(bool(cat.sel.getFlags('receptor').all()))


def test_partial_category_click_checks_all():
    s = session('1a3r.pdb', 'partial')
    cat = find(s.tree(), label='protein', kind='category')[0]
    s.click(cat.children[0].id, 'rec')
    check(cat.check_state(REC) == PARTIAL)
    s.click(cat.id, 'rec')          # Qt: partial -> checked
    check(cat.check_state(REC) == CHECKED)
    check(bool(cat.sel.getFlags('receptor').all()))


def test_altloc_choice_moves_the_flags_and_follows_tabs():
    s = session('6ggd.pdb', 'alt')
    t = s.tree()
    ser = find(t, key='A:SER:166:_')[0]
    B = ser.sel.select('altloc B')
    A = ser.sel.select('altloc A')
    check(bool(B.getFlags('receptor').all()) and not A.getFlags('receptor').any())
    # _altlocsRes lists the letters in an order that varies between runs
    # (set iteration); the default is the most populated one, B (0.58)
    check(ser.altloc['alts'][ser.altloc['current']] == 'B', ser.altloc)
    iA = ser.altloc['alts'].index('A')
    s.altloc(ser.id, iA)                                  # choose A
    check(bool(A.getFlags('receptor').all()) and not B.getFlags('receptor').any())
    check(s.profile._altlocsRes['A:SER:166:_'][0] == iA)
    s.switch('ASU')                                       # same residue, other tab
    ser2 = find(s.tree(), key='A:SER:166:_')[0]
    check(ser2.altloc['current'] == iA, ser2.altloc)
    A2, B2 = ser2.sel.select('altloc A'), ser2.sel.select('altloc B')
    check(bool(A2.getFlags('receptor').all()) and not B2.getFlags('receptor').any())
    d = s.save('receptor')
    path = os.path.join(s.folder, 'save_%03d_rec' % d['save'], d['files'][0]['name'])
    check(b'REMARK 850 AGFR altlocs A:SER:166:_ @A (0.42)' in read(path))


def test_edited_manual_charge_is_written():
    s = session('6ggd.pdb', 'charge')
    zn = [n for n in find(s.tree(), parent='manual charges') if n.key.startswith('A:ZN')][0]
    set_state(s, zn, REC, CHECKED)
    s.charge(zn.id, 1.5)
    d = s.save('receptor')
    path = os.path.join(s.folder, 'save_%03d_rec' % d['save'], d['files'][0]['name'])
    lines = [l for l in read(path).decode().splitlines()
             if l.startswith(('ATOM', 'HETATM')) and l[17:20] == ' ZN' and l[21] == 'A']
    check(len(lines) == 1, lines)
    check(float(lines[0][70:76]) == 1.5, lines[0])
    check(s.saves[-1]['charges'] == {zn.key: 1.5}, s.saves[-1]['charges'])


def test_sidechain_as_is_is_not_completed():
    s = session('6ggd.pdb', 'asis')
    lys = find(s.tree(), key='A:LYS:120:_')[0]
    s.sidechain(lys.id, 'asis')
    check('use as is missing:' in lys.label, lys.label)
    d = s.save('receptor')
    text = read(os.path.join(s.folder, 'save_%03d_rec' % d['save'], d['files'][0]['name'])).decode()
    check('side chain completed A:LYS:120:_' not in text)
    check('side chain completed A:ARG:209:_' in text)
    names = set(l[12:16].strip() for l in text.splitlines()
                if l.startswith('ATOM') and l[21] == 'A' and l[22:26].strip() == '120')
    check('NZ' not in names and 'CB' in names, names)


def test_mutation_after_save_is_honoured():
    """AGFRgui did not invalidate the processed BM on a mutation."""
    s = session('1a3r.pdb', 'mut_after')
    s.save('receptor')
    n = find(s.tree(), key='L:LEU:27:C')[0]
    s.mutate(n.id, 'ALA')
    check(n.label == "LEU'27 -> ALA", n.label)
    d = s.save('receptor')
    text = read(os.path.join(s.folder, 'save_%03d_rec' % d['save'], d['files'][0]['name'])).decode()
    check('REMARK 850 AGFR mutation L:LEU:27:C -> ALA' in text)
    s.mutate(n.id, 'Base Res: LEU')
    check(n.label == "LEU'27C", n.label)
    check('L:LEU:27:C' not in s.profile._resToMutate)


def test_session_record():
    import json
    s = session('1a3r.pdb', 'record')
    n = find(s.tree(), kind='residue', parent='protein')[0]
    s.click(n.id, 'lig')
    s.save('receptor')
    rec = json.load(open(os.path.join(s.folder, 'session.json')))
    check(rec['input']['name'] == '1a3r.pdb' and len(rec['input']['md5']) == 32)
    check(rec['environment']['adfr_dist'].startswith('adfr-'), rec['environment'])
    check(set(rec['environment']['adfr_md5']) >= {'utils/runAGFR.py', 'recFromPDB.py'})
    check([e['what'] for e in rec['events']] == ['check'], rec['events'])
    check(rec['saves'][0]['files'][0]['name'] == '1a3r_rec.pdbqt')


def test_receptor_records_its_biomolecule():
    """REMARK 850 AGFR biomolecule: the tab, its ASU chains, and the number of
    copies when the assembly is built by symmetry."""
    def header(s, tab):
        s.switch(tab)
        d = s.save('receptor')
        path = os.path.join(s.folder, 'save_%03d_rec' % d['save'], d['files'][0]['name'])
        return [l for l in read(path).decode().splitlines() if l.startswith('REMARK 850 AGFR')]
    s = session('6ggd.pdb', 'biomol_6ggd')
    for tab, want in (('ASU', 'ASU: chains A B'), ('Biomol 2', 'Biomol 2: chains B')):
        h = header(s, tab)
        check(h[1].startswith('REMARK 850 AGFR source') and
              h[2] == 'REMARK 850 AGFR biomolecule ' + want, (tab, h[:3]))
    s = session('1h6w.pdb', 'biomol_1h6w')
    h = header(s, 'Biomol 1')            # chains A and B, generated three times
    check(h[2] == 'REMARK 850 AGFR biomolecule Biomol 1: chains A B, 3 copies', h[:3])
    rec = s.saves[-1]
    check(rec['tab'] == 'Biomol 1', rec['tab'])


def test_ligand_records_its_biomolecule():
    """Every ligand file carries the same biomolecule line, after its own
    "REMARK 850 source" line: here two fragments (the EYB of each chain of the
    6GGD ASU), so two numbered files."""
    s = session('6ggd.pdb', 'biomol_lig')
    s.switch('ASU')
    for key in ('A:EYB:402:_', 'B:EYB:402:_'):
        set_state(s, find(s.tree(), key=key)[0], LIG, CHECKED)
    d = s.save('ligand')
    names = [f['name'] for f in d['files']]
    check(names == ['6ggd_lig_001.pdbqt', '6ggd_lig_002.pdbqt'], names)
    for name in names:
        lines = read(os.path.join(s.folder, 'save_%03d_lig' % d['save'], name)).decode().splitlines()
        i = [k for k, l in enumerate(lines) if 'biomolecule' in l]
        check(len(i) == 1 and lines[i[0]] == 'REMARK 850 AGFR biomolecule ASU: chains A B', (name, lines[:4]))
        check(lines[i[0] - 1].startswith('REMARK 850 source'), (name, lines[:4]))
