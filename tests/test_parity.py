# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Byte parity with ADFR's own headless code path.

Each test drives a web Session the way a user would (tab, checkboxes, menus)
and requires the PDBQT it writes to be BYTE-IDENTICAL to what ADFR writes
when given the same choices on its command line (the scenarios of
ADFR/tests/testAGFR.py), or, for ligands, to ADFR's processLigands run on
atoms flagged directly with prody, independently of the tree code.
"""
import os

from common import (CHECKED, LIG, PDBS, REC, UNCHECKED, check, find, outdir,
                    oracle_receptor, origins, quiet, read, session, set_state)


def _same(a, b):
    """ADapt's file `a` equals ADFR's `b`, byte for byte."""
    ba, bb = read(a), read(b)
    check(b'REMARK 850 AGFR biomolecule ' in ba, '%s records its biomolecule' % a)
    if ba != bb:
        la, lb = ba.decode().splitlines(), bb.decode().splitlines()
        for i, (x, y) in enumerate(zip(la, lb)):
            if x != y:
                raise AssertionError('%s vs %s differ at line %d:\n  web   : %s\n  oracle: %s'
                                     % (a, b, i + 1, x, y))
        raise AssertionError('%s (%d lines) vs %s (%d lines)' % (a, len(la), b, len(lb)))


def _save(s, kind):
    d = s.save(kind)
    folder = os.path.join(s.folder, 'save_%03d_%s' % (d['save'], kind[:3]))
    return [os.path.join(folder, f['name']) for f in d['files']]


def test_default_receptor_every_pdb():
    """Every ADFR test PDB: load, change nothing, write the receptor; the
    default biomolecule (Biomol 1, else the first) must match `agfr --toPdbqt`."""
    for pdb in PDBS:
        s = session(pdb, 'default')
        web = _save(s, 'receptor')[0]
        ref = oracle_receptor(pdb, outdir('oracle_default_' + pdb[:-4]))
        _same(web, ref)


def test_4a1x_remove_chain_C():
    """testAGFR.test1_1: protToRemove 'chid C' == unchecking chain C's
    categories in the tree."""
    s = session('4a1x.pdb', 'rm_C')
    t = s.tree()
    root = [r for r in t.roots if r.label == 'C'][0]
    for n in list(root.children) + [c for g in root.children if g.kind == 'group'
                                    for c in g.children]:
        if n.kind in ('category', 'group') and n.live and n.check_state(REC) not in (None, UNCHECKED):
            set_state(s, n, REC, UNCHECKED)
    web = _save(s, 'receptor')[0]
    ref = oracle_receptor('4a1x.pdb', outdir('oracle_rm_C'), protToRemove='chid C')
    _same(web, ref)


def test_6ggd_biomol2_hetero_and_mutations():
    """testAGFR.test2: Biomol 2, keep B:ZN401 and B:EYB402 in the receptor,
    mutate SER96 ILE, PRO98 ALA, GLN100 PHE and the altloc residues SER166 and
    CYS182 to ALA. Also checks the manual charge (ZN) passes through unedited."""
    s = session('6ggd.pdb', '6ggd_b2')
    s.switch('Biomol 2')
    t = s.tree()
    for key in ('B:ZN:401:_', 'B:EYB:402:_'):
        n = find(t, key=key)
        check(len(n) == 1, key)
        set_state(s, n[0], REC, CHECKED)
    for key, new in (('B:SER:96:_', 'ILE'), ('B:PRO:98:_', 'ALA'), ('B:GLN:100:_', 'PHE'),
                     ('B:SER:166:_', 'ALA'), ('B:CYS:182:_', 'ALA')):
        n = find(t, key=key)
        check(len(n) == 1 and n[0].menu == 'mutate', key)
        s.mutate(n[0].id, new)
    web = _save(s, 'receptor')[0]
    mutations = ('B:SER:96::ILE; B:PRO:98::ALA; B:GLN:100::PHE; '
                 'B:SER:166::ALA; B:CYS:182::ALA')
    ref = oracle_receptor('6ggd.pdb', outdir('oracle_6ggd_b2'), bioMol='2',
                          heteroToKeep='chid B resnum 402 401', mutateRes=mutations)
    _same(web, ref)


def test_1a3r_asu_icode_mutation():
    """testAGFR.test3: the only assembly, mutate L:LEU27C (insertion code C)."""
    s = session('1a3r.pdb', '1a3r_asu')
    t = s.tree()
    n = find(t, key='L:LEU:27:C')
    check(len(n) == 1, 'L:LEU:27:C in tree')
    s.mutate(n[0].id, 'ALA')
    web = _save(s, 'receptor')[0]
    ref = oracle_receptor('1a3r.pdb', outdir('oracle_1a3r'), bioMol='ASU',
                          mutateRes='L:LEU:27:C:ALA')
    _same(web, ref)


def test_tab_switch_writes_the_new_biomolecule():
    """AGFRgui kept _BMProcessed across tab switches; here a save after a
    switch must write the new biomolecule (== agfr --bioMol 2)."""
    s = session('6ggd.pdb', 'switch')
    _save(s, 'receptor')                     # Biomol 1 first
    s.switch('Biomol 2')
    web = _save(s, 'receptor')[0]
    ref = oracle_receptor('6ggd.pdb', outdir('oracle_switch'), bioMol='2')
    _same(web, ref)


def _oracle_ligands(pdb, tab, chid, folder):
    """ADFR's ligand path on atoms flagged with prody directly: chain `chid`
    of segment A, no water, no unsupported atom types, the chosen altloc."""
    from adapt.session import gui_options
    from ADFR.utils.runAGFR import runAGFR
    from common import DATA
    with quiet():
        a = runAGFR()
        a.proccessPreLoadCmdLineOptions(gui_options())
        a.loadReceptor(os.path.join(DATA, pdb), processReceptor=False)
        p = a.PDBprofile
        p._currentBM = bm = p._assemblies[tab][1]
        ag = bm._ag
        sel = ag.select('segment A and chid %s and not water and not unsre' % chid)
        sel.setFlags('receptor', False)
        sel.setFlags('ligand', False)
        plain = sel.select('not altlocs')
        if plain is not None:
            plain.setFlags('ligand', True)
        for key, (ind, alts, occ) in p._altlocsRes.items():
            c, resname, resnum, icode = key.split(':')
            if c != chid:
                continue
            res = sel.select('altlocs and chid %s resnum `%s` icode %s' % (c, resnum, icode))
            if res is None:
                continue
            chosen = res.select('altloc _ %s' % alts[ind])
            if chosen is not None:
                chosen.setFlags('ligand', True)
        mol = a.processBM()
        a.processLigands(mol.select('ligand'))
        paths = []
        base = os.path.join(folder, p._pdbid + '_lig')
        many = len(a._adLigs) > 1
        for num, lig in enumerate(a._adLigs):
            path = '%s_%03d.pdbqt' % (base, num + 1) if many else base + '.pdbqt'
            with open(path, 'w') as fh:
                for l in lig.getPDBQTlines():
                    fh.write('%s\n' % l)
            paths.append(path)
    return paths


def test_6bhg_chain_B_as_ligand():
    """testAGFRgui.test_loadPDBSaveLigand: every chain-B category except water
    checked as ligand, then 'lig. ...'."""
    s = session('6bhg.pdb', '6bhg_lig')
    t = s.tree()
    root = [r for r in t.roots if r.label == 'B'][0]
    targets = []
    for n in root.children:
        if n.kind == 'category' and n.label != 'water':
            targets.append(n)
        if n.kind == 'group' and n.label == 'altloc':
            targets.extend(c for c in n.children if c.label != 'water')
    check(targets, 'chain B has categories')
    for n in targets:
        set_state(s, n, LIG, CHECKED)
    web = _save(s, 'ligand')
    ref = _oracle_ligands('6bhg.pdb', s.current, 'B', outdir('oracle_6bhg'))
    check([os.path.basename(p) for p in web] == [os.path.basename(p) for p in ref],
          '%s vs %s' % (web, ref))
    for a, b in zip(web, ref):
        _same(a, b)
    # the peptide's residues, in file order, four keys per line
    raw = read(web[0])
    keys = origins(raw)
    check(keys and all(k.startswith('B:') for k in keys) and len(set(keys)) == len(keys), keys)
    lines = [l for l in raw.split(b'\n') if l.startswith(b'REMARK 850 AGFR ligand from ')]
    check(len(lines) == (len(keys) + 3) // 4, (len(keys), len(lines)))


def test_ion_charge_matches_command_line():
    """An edited ion charge == agfr --ionCharges, header included (the command
    line used to report the default charge, 2.0000, in its custom-charge
    REMARK; ADFR 0.0.23 reports the one written)."""
    s = session('6ggd.pdb', 'ioncharge')
    s.switch('ASU')
    zn = find(s.tree(), key='A:ZN:401:_')[0]
    set_state(s, zn, REC, CHECKED)
    s.charge(zn.id, 1.5)
    web = read(_save(s, 'receptor')[0]).decode().splitlines()
    ref = read(oracle_receptor('6ggd.pdb', outdir('oracle_ioncharge'), bioMol='ASU',
                               heteroToKeep='chid A resnum 401',
                               ionCharges='A:ZN:401::1.5')).decode().splitlines()
    diff = [(a, b) for a, b in zip(web, ref) if a != b]
    check(len(web) == len(ref) and not diff, (len(web), len(ref), diff[:3]))
    check('REMARK 850 AGFR custom charge for A:ZN:401: = 1.5000' in web)
    zn_line = [l for l in web if l[17:20] == ' ZN'][0]
    check(float(zn_line[70:76]) == 1.5, zn_line)


def test_side_chain_as_is_matches_command_line():
    """A side chain kept as is == agfr --sidechainsAsIs, header included (the
    command line had no option for it, so such a receptor could not be rebuilt
    from the command line)."""
    s = session('6ggd.pdb', 'asis')
    s.switch('ASU')
    s.sidechain(find(s.tree(), key='A:LYS:120:_')[0].id, 'asis')
    web = read(_save(s, 'receptor')[0]).decode().splitlines()
    ref = read(oracle_receptor('6ggd.pdb', outdir('oracle_asis'), bioMol='ASU',
                               sidechainsAsIs='A:LYS:120:')).decode().splitlines()
    diff = [(a, b) for a, b in zip(web, ref) if a != b]
    check(len(web) == len(ref) and not diff, (len(web), len(ref), diff[:3]))
    check('REMARK 850 AGFR side chain completed A:LYS:120:_' not in web)
