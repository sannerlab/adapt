# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The tree mirrors ReceptorItemsTree.buildUI."""
from common import (CHECKED, PDBS, REC, UNCHECKED, all_nodes, check, find,
                    session)


def test_tabs_are_the_assemblies():
    for pdb in PDBS:
        s = session(pdb, 'tabs')
        names = list(s.profile._assemblies)
        check(list(s.trees) == names, pdb)
        want = 'Biomol 1' if 'Biomol 1' in names else names[0]
        check(s.current == want, '%s opens on %s, not %s' % (pdb, want, s.current))


def test_every_supported_atom_is_in_exactly_one_leaf():
    """Each segment-A atom of the tab's chains that has no unsupported atom
    type sits in exactly one residue leaf (plain or altloc), for every tab."""
    for pdb in PDBS:
        s = session(pdb, 'cover')
        for name in s.trees:
            s.switch(name)
            t = s.tree()
            ag = t._ag
            want = set()
            for chid in t._chaindInBM:
                sel = ag.select('segment A and chid %s and not unsre' % (chid if chid != ' ' else '_'))
                if sel is not None:
                    want.update(int(i) for i in sel.getIndices())
            seen = {}
            for n in all_nodes(t):
                if n.kind == 'residue' and n.live:
                    for i in n.sel.getIndices():
                        seen[int(i)] = seen.get(int(i), 0) + 1
            dup = [i for i, c in seen.items() if c > 1]
            check(not dup, '%s %s: %d atoms in two leaves' % (pdb, name, len(dup)))
            missing = want - set(seen)
            extra = set(seen) - want
            check(not missing and not extra, '%s %s: %d atoms missing, %d extra' % (
                pdb, name, len(missing), len(extra)))


def test_leaf_states_follow_the_flags():
    """At build, a plain residue shows rec checked iff some atom is receptor,
    and lig checked iff some atom is ligand (which unchecks rec)."""
    for pdb in PDBS:
        s = session(pdb, 'states')
        for n in all_nodes(s.tree()):
            if n.kind != 'residue' or not n.live or n.altloc is not None:
                continue
            rec = n.sel.select('receptor') is not None
            lig = n.sel.select('ligand') is not None
            check(n.state[1] == (CHECKED if lig else UNCHECKED), n.label)
            check(n.state[REC] == (CHECKED if rec and not lig else UNCHECKED), n.label)


def test_categories_labels_and_defaults_6ggd():
    s = session('6ggd.pdb', 'labels')
    t = s.tree()
    A = [r for r in t.roots if r.label == 'A'][0]
    cats = [c.label for c in A.children if c.kind == 'category']
    check(cats[:2] == ['protein', 'sidechains to complete'], cats)
    check([c.label for c in A.children if c.kind == 'group'] == ['altloc'])
    lys = find(t, key='A:LYS:120:_')[0]
    check(lys.label == "LYS'120 build:CG CD CE NZ", lys.label)
    check(lys.menu == 'sc' and lys.rotamer is not None)
    check(0 <= lys.rotamer.best < lys.rotamer.n())
    ser = find(t, key='A:SER:166:_')[0]
    alt = dict(zip(ser.altloc['alts'], ser.altloc['occ']))   # letter order varies
    check(alt == {'B': 0.58, 'A': 0.42}, ser.altloc)
    check(ser.altloc['alts'][ser.altloc['current']] == 'B', ser.altloc)
    check(ser.menu == 'mutate')
    prot = [c for c in A.children if c.label == 'protein'][0]
    check(prot.check_state(REC) == CHECKED)
    water = [c for c in A.children if c.label == 'water'][0]
    check(water.check_state(REC) == UNCHECKED)
    manch = [n for n in all_nodes(t) if n.charge is not None]
    check(manch and all(n.charge == 2.0 for n in manch if 'ZN' in n.label),
          [(n.label, n.charge) for n in manch])
