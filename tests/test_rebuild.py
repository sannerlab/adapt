# Author: Aadhil Haq. Copyright (c) 2026.
"""Missing residues rebuilt on request (adapt/rebuild.py, adapt/loops.py).

Geometry is measured here on the files written, not taken from the code that
built them, and printed: C to N at the joins and inside the loop, CA to CA of
consecutive residues, the closest rebuilt to observed heavy-atom distance, and
the seconds a build took. Files derived from ADFR's test entries (a stretch of
1H6W moved into REMARK 465; 6BHG with its gap pulled apart; 4A1X with MSE, SEP
or a clashing number in its gap) are made by the functions below.
"""
import collections
import hashlib
import json
import os
import subprocess
import sys

import numpy

from common import DATA, HERE, check, find, outdir, quiet, read, session
from adapt import env, loops, prep
from adapt.session import Session

TOL = loops.TOLERANCE


# ---- helpers -----------------------------------------------------------------
def _gap_node(s, key):
    n = [x for x in find(s.tree(), kind='gap') if x.rebuild and x.rebuild['key'] == key]
    check(len(n) == 1, '%s in the tree of %s' % (key, s.current))
    return n[0]


def _save(s):
    d = s.save('receptor')
    folder = os.path.join(s.folder, 'save_%03d_rec' % d['save'])
    return [os.path.join(folder, f['name']) for f in d['files']]


def _atoms(path):
    """[(chain, resnum, icode, resname, name, xyz, line)] of the heavy atoms of a
    PDBQT or PDB file, in file order."""
    out = []
    for l in read(path).decode().splitlines():
        if not l.startswith(('ATOM', 'HETATM')):
            continue
        if path.endswith('.pdbqt'):
            heavy = l[77:79].strip() not in ('H', 'HD', 'HS')
        else:
            heavy = (l[76:78].strip() or l[12:16].strip()[:1]) != 'H'
        if heavy:
            out.append((l[21], int(l[22:26]), l[26].strip(), l[17:20].strip(), l[12:16].strip(),
                        numpy.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l))
    return out


def _residues(atoms, chain):
    """OrderedDict (resnum, icode) -> (resname, {name: xyz}) of one chain, in
    file order (of its first atom: ADFR writes the side chains it completes
    after the rest of the chain)."""
    out = collections.OrderedDict()
    for c, num, icode, resname, name, xyz, line in atoms:
        if c == chain:
            out.setdefault((num, icode), (resname, {}))[1][name] = xyz
    return out


def _blocks(atoms, chain, ids):
    """How many runs of consecutive atoms each residue id of `ids` has in the
    file: 1 when no other residue shares it."""
    runs, last = collections.Counter(), None
    for c, num, icode, resname, name, xyz, line in atoms:
        key = (num, icode) if c == chain else None
        if key != last and key in ids:
            runs[key] += 1
        last = key
    return runs


def _geometry(res, keys):
    """C to N and CA to CA of consecutive residues `keys`, in A."""
    cn = [numpy.linalg.norm(res[a][1]['C'] - res[b][1]['N']) for a, b in zip(keys[:-1], keys[1:])]
    caca = [numpy.linalg.norm(res[a][1]['CA'] - res[b][1]['CA']) for a, b in zip(keys[:-1], keys[1:])]
    return cn, caca


def _contact(atoms, chain, new, flanks):
    """Closest distance between a heavy atom of the new residues and an observed
    heavy atom outside the two residues they are bonded to."""
    a = numpy.array([x[5] for x in atoms if x[0] == chain and (x[1], x[2]) in new])
    b = numpy.array([x[5] for x in atoms if not (x[0] == chain and (x[1], x[2]) in new | flanks)])
    return float(numpy.sqrt(((a[:, None] - b[None]) ** 2).sum(2)).min())


TABLE = []


def _report(entry, chain, gap, n, cn, caca, contact, seconds):
    row = '%-14s %-5s %-6s %2d  C-N %.3f..%.3f  CA-CA %.3f..%.3f  contact %.2f  %s' % (
        entry, chain, gap, n, min(cn), max(cn), min(caca), max(caca), contact, seconds)
    TABLE.append(row)
    print('      ' + row)


def _check_loop(path, chain, before, new, after, entry='', label='', seconds=''):
    """The written loop: residues in order, numbered as asked, geometry in tolerance."""
    atoms = _atoms(path)
    res = _residues(atoms, chain)
    keys = list(res)
    i = keys.index(before)
    check(keys[i:i + len(new) + 2] == [before] + [k for k, name in new] + [after],
          '%s: %s' % (path, keys[i:i + len(new) + 2]))
    check([res[k][0] for k, name in new] == [name for k, name in new], [res[k][0] for k, name in new])
    runs = _blocks(atoms, chain, set(k for k, name in new))
    check(set(runs.values()) == {1}, 'a new residue id is used by another residue: %s' % runs)
    cn, caca = _geometry(res, keys[i:i + len(new) + 2])
    contact = _contact(atoms, chain, set(k for k, n in new), {before, after})
    if entry:
        _report(entry, chain, label, len(new), cn, caca, contact, seconds)
    lo, hi = TOL['C-N']
    check(lo <= min(cn) and max(cn) <= hi, 'C to N %s' % cn)
    lo, hi = TOL['CA-CA']
    check(lo <= min(caca) and max(caca) <= hi, 'CA to CA %s' % caca)
    check(contact >= TOL['contact'], 'contact %.2f' % contact)
    return res


def _seconds(s, segment='A'):
    """'4.8 s, seed 1': how long the build of a rebuilt site took (PDBFixer and
    the refinement, every try; the sites of a group share one build), and the
    seed of the loop kept."""
    r = [r for r in s.agfr._BMProcessed.rebuilt if r['segment'] == segment][0]
    return '%5.1f s, seed %d' % (r['seconds']['total'], r['seed'])


# ---- the buildable gaps of ADFR's test entries ---------------------------------
def test_4a1x_gap_rebuilt_in_chain_order():
    """4A1X Biomol 1, gap SER -3 GKGG PRO 2: the residues REMARK 465 lists, in
    order, numbered -2 to 1; written in the PDBQT and in _rec_H.pdb; the REMARK;
    every observed heavy atom where it was; the unchosen gap of chain B not built."""
    s = session('4a1x.pdb', 'rebuild_4a1x')
    plain_q, plain_h = _save(s)                  # Biomol 1, nothing rebuilt
    node = _gap_node(s, 'A:SER:-3:_ A:PRO:2:_')
    check(node.rebuild == {'key': 'A:SER:-3:_ A:PRO:2:_', 'on': False, 'why': None}, node.rebuild)
    d = s.rebuild(node.id, True)
    check(d['refused'] == [] and node.rebuild['on'] and d['rebuild'] == 1, d['refused'])
    check([r['gap'] for r in d['rebuilt']] == ['A:SER:-3:_ A:PRO:2:_'], d['rebuilt'])
    q, h = _save(s)
    new = [((-2, ''), 'GLY'), ((-1, ''), 'LYS'), ((0, ''), 'GLY'), ((1, ''), 'GLY')]
    _check_loop(q, 'A', (-3, ''), new, (2, ''), '4A1X Biomol 1', 'GKGG', _seconds(s))
    _check_loop(h, 'A', (-3, ''), new, (2, ''))
    text = read(q).decode()
    remark = ('REMARK 850 AGFR rebuilt A:GLY:-2:_ to A:GLY:1:_ (GKGG) using PDBFixer %s and '
              'Amber ff14SB, seed %d' % (loops.versions()[0], s.agfr._BMProcessed.rebuilt[0]['seed']))
    check(remark in text.splitlines(), [l for l in text.splitlines() if 'rebuilt' in l])
    lines = text.splitlines()
    check(lines.index(remark) + 1 == lines.index('REMARK 850 AGFR protonated using reduce'),
          'the rebuilt line comes before the protonation line')
    # every observed heavy atom is where it was, and the new ones are all that was added
    before = dict(((a[0], a[1], a[2], a[4]), a[5]) for a in _atoms(plain_q))
    after = dict(((a[0], a[1], a[2], a[4]), a[5]) for a in _atoms(q))
    moved = [k for k, x in before.items() if k not in after or not numpy.array_equal(x, after[k])]
    check(not moved, 'observed atoms moved or lost: %s' % moved[:5])
    added = set(after) - set(before)
    check(set((k[0], k[1], k[2]) for k in added) == set(('A',) + k for k, n in new), sorted(added)[:5])
    print('      %d observed heavy atoms unchanged (0.000 A), %d added' % (len(before), len(added)))
    check(not any(a[0] == 'B' for a in _atoms(q)), 'Biomol 1 has no chain B')
    # chain B's gap is in another tab and was not chosen
    s.switch('Biomol 2')
    q2, h2 = _save(s)
    res_b = _residues(_atoms(q2), 'B')
    check((-2, '') not in res_b and (-1, '') not in res_b, 'unchosen gap not built')
    check('rebuilt' not in read(q2).decode(), 'no rebuilt line in Biomol 2')


def test_1a3r_numbers_come_from_remark_465():
    """1A3R H: CYS 128, GDTTG, SER 136. The numbers on from 128 would be 129 to
    133; REMARK 465 numbers them 131 to 135, and that is what is written."""
    s = session('1a3r.pdb', 'rebuild_1a3r')
    node = _gap_node(s, 'H:CYS:128:_ H:SER:136:_')
    d = s.rebuild(node.id, True)
    check(d['refused'] == [], d['refused'])
    gap = s.rebuilder.gaps['H:CYS:128:_ H:SER:136:_']
    check(gap.residues == [('GLY', 131, ''), ('ASP', 132, ''), ('THR', 133, ''), ('THR', 134, ''),
                           ('GLY', 135, '')], gap.residues)
    q, h = _save(s)
    new = [((131, ''), 'GLY'), ((132, ''), 'ASP'), ((133, ''), 'THR'), ((134, ''), 'THR'),
           ((135, ''), 'GLY')]
    _check_loop(q, 'H', (128, ''), new, (136, ''), '1A3R ASU', 'GDTTG', _seconds(s))


def test_6bhg_gap_rebuilt():
    """6BHG A: ASP 270, GNQ, VAL 274."""
    s = session('6bhg.pdb', 'rebuild_6bhg')
    node = _gap_node(s, 'A:ASP:270:_ A:VAL:274:_')
    d = s.rebuild(node.id, True)
    check(d['refused'] == [], d['refused'])
    q, h = _save(s)
    new = [((271, ''), 'GLY'), ((272, ''), 'ASN'), ((273, ''), 'GLN')]
    _check_loop(q, 'A', (270, ''), new, (274, ''), '6BHG ASU', 'GNQ', _seconds(s))


def test_unclosable_gap_is_refused_with_its_reason():
    """4A1X B: two residues (GK) cannot join SER -3 to GLY 0 with trans peptide
    bonds and sane bonds and angles (an exhaustive search over their four
    torsions, ideal geometry, gets no closer than 1.4 A RMSD on the residue
    after the gap; observed two-residue stretches close to 0.1 A). The loop is
    built, fails the checks, and is refused: the box goes back, the reason is in
    the tree and in session.json, nothing is written."""
    s = session('4a1x.pdb', 'rebuild_refused')
    s.switch('Biomol 2')
    node = _gap_node(s, 'B:SER:-3:_ B:GLY:0:_')
    d = s.rebuild(node.id, True)
    check(len(d['refused']) == 1 and 'passed the checks' in d['refused'][0], d['refused'])
    why = node.rebuild['why']
    check(node.rebuild['on'] is False and why.startswith(
        'none of %d loops built for it passed the checks; the first had ' % len(loops.SEEDS)), why)
    check('LYS -1' in why or 'GLY 0' in why or 'SER -3' in why, 'faults name the real residues')
    print('      4A1X B GK refused: %s' % node.rebuild['why'])
    check('B:SER:-3:_ B:GLY:0:_' not in s.rebuilder.chosen and d['rebuilt'] == [])
    rec = json.load(open(os.path.join(s.folder, 'session.json')))
    ev = [e for e in rec['events'] if e['what'] == 'rebuild']
    check(len(ev) == 1 and ev[0]['rebuild'] is False and ev[0]['asked'] and ev[0]['refused'], ev)
    q, h = _save(s)
    check('rebuilt' not in read(q).decode() and (-2, '') not in _residues(_atoms(q), 'B'))
    tree = s.view(full=True)['tree']
    boxes = [n for r in tree['roots'] for n in _walk_json(r) if n.get('rebuild')]
    check(boxes and boxes[0]['rebuild']['why'] == node.rebuild['why'], boxes)


def _walk_json(n):
    yield n
    for c in n.get('children', []):
        for x in _walk_json(c):
            yield x


# ---- what is offered, and what is not -----------------------------------------
def _derived(name, edit, folder):
    """A copy of ADFR's test file `name` with every line passed through
    edit(line) (None drops it, a list replaces it), written into this run's
    outputs (derived/<folder>); its path."""
    out = []
    for l in read(os.path.join(DATA, name)).decode().splitlines():
        x = edit(l)
        if x is None:
            continue
        out.extend(x if isinstance(x, list) else [x])
    path = os.path.join(outdir(os.path.join('derived', folder)), name)
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return path


_N = [0]


def _open(path, name):
    _N[0] += 1
    with quiet():
        return Session('d%015d' % _N[0], read(path), os.path.basename(path), outdir(name),
                       env.provenance())


def _4a1x_with(resname):
    """4A1X with chain A's LYS -1 (missing, in REMARK 465) written as `resname`
    in its SEQRES and REMARK 465 records."""
    seen = []

    def edit(l):
        if l.startswith('SEQRES') and l[11] == 'A' and ' LYS ' in l and not seen:
            seen.append(1)
            return l.replace(' LYS ', ' %s ' % resname, 1)
        if l.startswith('REMARK 465     LYS A    -1'):
            return l.replace('LYS', resname, 1)
        return l
    return _derived('4a1x.pdb', edit, resname)


def test_residue_names_come_from_remark_465():
    """SEQRES's one-letter codes cannot tell SEP from SER (ProDy writes both S),
    so names come from REMARK 465: a selenomethionine is built as MET, as ADFR
    writes MSE; a phosphoserine is not one of the 20, and the gap is not offered."""
    s = _open(_4a1x_with('MSE'), 'names_mse')
    gap = s.rebuilder.gaps['A:SER:-3:_ A:PRO:2:_']
    check(gap.why is None and gap.residues == [('GLY', -2, ''), ('MET', -1, ''), ('GLY', 0, ''),
                                               ('GLY', 1, '')], (gap.why, gap.residues))
    s = _open(_4a1x_with('SEP'), 'names_sep')
    gap = s.rebuilder.gaps['A:SER:-3:_ A:PRO:2:_']
    check(gap.why == 'SEP -1 is not one of the 20 amino acids', gap.why)
    node = _gap_node(s, 'A:SER:-3:_ A:PRO:2:_')
    try:
        s.rebuild(node.id, True)
        raise AssertionError('a gap with SEP was rebuilt')
    except Exception as e:
        check('not one of the 20' in str(e), e)


def test_numbers_already_used_are_refused():
    """A water of chain A numbered -1, as the missing LYS: the gap is not offered."""
    def edit(l):
        if l.startswith('HETATM') and l[17:20] == 'HOH' and l[21] == 'A' and not edit.done:
            edit.done = True
            return l[:22] + '%4d' % -1 + l[26:]
        return l
    edit.done = False
    s = _open(_derived('4a1x.pdb', edit, 'water'), 'numbers_used')
    why = s.rebuilder.gaps['A:SER:-3:_ A:PRO:2:_'].why
    check(why == 'REMARK 465 numbers LYS -1, which the chain already uses', why)


def test_unspannable_gap_is_refused():
    """6BHG with chain A from VAL 274 on moved 20 A: three residues cannot span
    the CA to CA distance of ASP 270 to VAL 274; the gap is not offered."""
    def edit(l):
        if l.startswith(('ATOM', 'HETATM')) and l[21] == 'A' and l[17:20] != 'HOH' \
                and int(l[22:26]) >= 274:
            return l[:30] + '%8.3f' % (float(l[30:38]) + 20.0) + l[38:]
        return l
    s = _open(_derived('6bhg.pdb', edit, 'apart'), 'unspannable')
    why = s.rebuilder.gaps['A:ASP:270:_ A:VAL:274:_'].why
    check(why is not None and why.startswith('its ends are') and
          why.endswith('3 residues span at most 15.2 A'), why)
    print('      6BHG, ends moved apart: %s' % why)
    node = _gap_node(s, 'A:ASP:270:_ A:VAL:274:_')
    try:
        s.rebuild(node.id, True)
        raise AssertionError('an unspannable gap was rebuilt')
    except Exception as e:
        check('its ends are' in str(e), e)


def test_termini_are_not_offered():
    s = session('4a1x.pdb', 'termini')
    s.switch('ASU')
    for n in find(s.tree(), kind='gap'):
        check((n.rebuild is not None) == n.label.startswith('gap:'), (n.label, n.rebuild))


def test_nothing_chosen_writes_master_files():
    """With nothing chosen (or chosen and unchosen again) the processed molecule
    is processBM's own and the files are byte for byte those of a session that
    never saw the rebuild box (the parity tests compare those with ADFR)."""
    s = session('4a1x.pdb', 'parity_none')
    q0, h0 = _save(s)
    with quiet():
        mol = prep.processed_bm(s)
    check(not hasattr(mol, 'rebuilt') and mol is s.agfr._BMProcessed, 'processBM molecule itself')
    node = _gap_node(s, 'A:SER:-3:_ A:PRO:2:_')
    s.rebuild(node.id, True)
    s.rebuild(node.id, False)
    q1, h1 = _save(s)
    t = session('4a1x.pdb', 'parity_fresh')
    q2, h2 = _save(t)
    for a, b in ((q0, q2), (q1, q2), (h0, h2), (h1, h2)):
        check(read(a) == read(b), '%s differs from %s' % (a, b))


# ---- provenance, reopening, reduce ----------------------------------------------
def test_record_and_reopened_session():
    """session.json records the choice as an edit and the build with the save;
    a page that opens the session again (GET api/session) shows the box ticked,
    and api/rebuilt draws the same coordinates."""
    import threading
    import urllib.request
    from adapt.server import build
    httpd = build('127.0.0.1', 0, outdir('rebuild_server'))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base, t = 'http://127.0.0.1:%d' % httpd.server_address[1], '?t=' + httpd.app.token

    def get(path):
        with urllib.request.urlopen(base + path) as r:
            return json.loads(r.read().decode())

    def post(path, obj=None, data=None):
        req = urllib.request.Request(base + path, method='POST',
                                     data=data if data is not None else json.dumps(obj).encode())
        with urllib.request.urlopen(req) as r:
            return json.loads(r.read().decode())
    try:
        up = post('/api/upload' + t + '&name=4a1x.pdb', data=read(os.path.join(DATA, '4a1x.pdb')))
        sid = up['session']
        gap = [n for r in up['tree']['roots'] for n in _walk_json(r) if n.get('rebuild')][0]
        check(gap['rebuild'] == {'key': 'A:SER:-3:_ A:PRO:2:_', 'on': False, 'why': None}, gap)
        d = post('/api/rebuild' + t, {'s': sid, 'node': gap['id'], 'on': True})
        check(d['refused'] == [] and len(d['rebuilt']) == 1, d.get('error'))
        first = d['rebuilt'][0]['pdb']
        sv = post('/api/save' + t, {'s': sid, 'kind': 'receptor'})
        again = get('/api/session' + t + '&s=' + sid)
        box = [n for r in again['tree']['roots'] for n in _walk_json(r) if n.get('rebuild')][0]
        check(box['rebuild']['on'] is True, box)
        check(get('/api/rebuilt' + t + '&s=' + sid)['rebuilt'][0]['pdb'] == first,
              'the reopened session draws the same residues')
        s = httpd.app.session(sid)
        rec = json.load(open(os.path.join(s.folder, 'session.json')))
        ev = [e for e in rec['events'] if e['what'] == 'rebuild']
        check(ev == [{'node': 'gap: SER-3--GKGG--PRO2', 'key': 'A:SER:-3:_ A:PRO:2:_',
                      'rebuild': True, 'what': 'rebuild', 'tab': 'Biomol 1', 'at': ev[0]['at']}], ev)
        r = rec['saves'][-1]['rebuilt']
        check(len(r) == 1 and r[0]['residues'] == ['A:GLY:-2:_', 'A:LYS:-1:_', 'A:GLY:0:_', 'A:GLY:1:_']
              and r[0]['seed'] in loops.SEEDS and r[0]['platform'] == 'Reference'
              and r[0]['method'] == 'PDBFixer %s and Amber ff14SB' % loops.versions()[0]
              and r[0]['openmm'] == loops.versions()[1], r)
        check(sv['save'] == 1, sv)
    finally:
        httpd.shutdown()


def test_rec_h_pdb_matches_the_pdbqt():
    """_rec_H.pdb with rebuilt residues: the same header, the same heavy atoms,
    every polar H of the PDBQT; reduce protonated the new residues and saw no
    chain end beside the gap (its residues are in chain order)."""
    s = session('6bhg.pdb', 'rebuild_reduce')
    s.rebuild(_gap_node(s, 'A:ASP:270:_ A:VAL:274:_').id, True)
    q, h = _save(s)
    qt, ht = read(q).decode().splitlines(), read(h).decode().splitlines()
    head = [l for l in qt if l.startswith('REMARK 850 AGFR')]
    check(ht[:len(head)] == head and any('rebuilt' in l for l in head), 'same header')
    xyz = lambda l: tuple(round(float(l[i:i + 8]), 3) for i in (30, 38, 46))
    hq = [l for l in qt if l.startswith(('ATOM', 'HETATM'))]
    hh = [l for l in ht if l.startswith(('ATOM', 'HETATM'))]
    heavy_q = sorted(xyz(l) for l in hq if l[77:79].strip() not in ('H', 'HD', 'HS'))
    heavy_h = sorted(xyz(l) for l in hh if l[76:78].strip() != 'H')
    check(heavy_q == heavy_h, (len(heavy_q), len(heavy_h)))
    polar = set(xyz(l) for l in hq if l[77:79].strip() in ('HD', 'HS'))
    check(polar <= set(xyz(l) for l in hh if l[76:78].strip() == 'H'), 'polar H of the PDBQT in _rec_H')
    for num in (271, 272, 273):
        names = [l[12:16].strip() for l in hh if l[21] == 'A' and int(l[22:26]) == num]
        check('H' in names and len([x for x in names if x.startswith('H')]) >= 3, (num, names))
    # reduce adds H1 H2 H3 to a chain's first residue; none beside the gap
    near = [l[12:16].strip() for l in hh if l[21] == 'A' and 270 <= int(l[22:26]) <= 274]
    check(not set(near) & {'H1', 'H2', 'H3', 'OXT'}, 'no chain end beside the gap')


# ---- reproducible -----------------------------------------------------------
_TWICE = r'''
import hashlib, sys
sys.path.insert(0, %r)
from common import find, quiet, session
s = session('4a1x.pdb', 'twice')
n = [x for x in find(s.tree(), kind='gap') if x.rebuild][0]
s.rebuild(n.id, True)
mol = s.agfr._BMProcessed
a, b = mol.rebuilt[0]['atoms']
print('SHA ' + hashlib.sha1(mol._ag.getCoords()[a:b].astype('<f8').tobytes()).hexdigest())
'''


def test_same_seed_same_coordinates_in_another_process():
    """The same input and seed give bit-identical coordinates in two separate
    processes, and in this one; PDBFixer's dynamics run on the Reference
    platform, and OPENMM_DEFAULT_PLATFORM is left as it was."""
    import openmm as mm
    made, real = [], mm.Context

    def recording(*a):
        c = real(*a)
        made.append(c.getPlatform().getName())
        return c
    s = session('4a1x.pdb', 'twice')
    n = [x for x in find(s.tree(), kind='gap') if x.rebuild][0]
    before = os.environ.get('OPENMM_DEFAULT_PLATFORM')
    mm.Context = recording
    try:
        s.rebuild(n.id, True)
    finally:
        mm.Context = real
    check(made and set(made) == {'Reference'}, made)
    check(os.environ.get('OPENMM_DEFAULT_PLATFORM') == before, 'environment left as it was')
    mol = s.agfr._BMProcessed
    a, b = mol.rebuilt[0]['atoms']
    here = hashlib.sha1(mol._ag.getCoords()[a:b].astype('<f8').tobytes()).hexdigest()
    shas = []
    for _ in range(2):
        out = subprocess.run([sys.executable, '-c', _TWICE % HERE], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, universal_newlines=True).stdout
        line = [l for l in out.splitlines() if l.startswith('SHA ')]
        check(line, out[-2000:])
        shas.append(line[0][4:])
    print('      coordinates sha1: this process %s, two others %s %s' % (here[:12], shas[0][:12],
                                                                         shas[1][:12]))
    check(shas == [here, here], (here, shas))


# ---- assemblies -------------------------------------------------------------
GAP_1H6W = (385, 387)


def _1h6w_gap(first=GAP_1H6W[0], last=GAP_1H6W[1]):
    """1H6W with chain A's residues first..last taken out of the coordinates and
    listed in REMARK 465 (with the file's other missing residues), as a gap the
    depositors would have reported; and their coordinates, to compare with."""
    taken = collections.OrderedDict()

    def edit(l):
        if l.startswith(('ATOM', 'ANISOU')) and l[21] == 'A' and first <= int(l[22:26]) <= last:
            if l.startswith('ATOM'):
                taken.setdefault(int(l[22:26]), [l[17:20], {}])[1][l[12:16].strip()] = \
                    numpy.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            return None
        if l.startswith('REMARK 465     VAL A   245'):         # the last of chain A's
            return [l] + [('REMARK 465     %s A %5d' % (taken_names[num], num)).ljust(80)
                          for num in range(first, last + 1)]
        return l
    taken_names = {}
    for l in read(os.path.join(DATA, '1h6w.pdb')).decode().splitlines():
        if l.startswith('ATOM') and l[21] == 'A' and first <= int(l[22:26]) <= last:
            taken_names[int(l[22:26])] = l[17:20]
    return _derived('1h6w.pdb', edit, 'gap'), taken


def test_copies_follow_their_chain():
    """1H6W Biomol 1 is three copies of chains A and B. A gap of chain A, chosen
    once, is rebuilt in every copy, each in its own place: one REMARK line, one
    record per copy, the residues in chain order in each copy. In the ASU tab
    the same choice gives one copy."""
    path, crystal = _1h6w_gap()
    first, last = GAP_1H6W
    s = _open(path, 'assembly')
    key = [k for k in s.rebuilder.gaps if k.startswith('A:') and k.split()[0].endswith(':%d:_' % (first - 1))]
    check(len(key) == 1, list(s.rebuilder.gaps))
    key = key[0]
    s.switch('Biomol 1')
    d = s.rebuild(_gap_node(s, key).id, True)
    check(d['refused'] == [], d['refused'])
    check(sorted(r['segment'] for r in d['rebuilt']) == ['A', 'B', 'C'], [r['segment'] for r in d['rebuilt']])
    q, h = _save(s)
    text = read(q).decode()
    check(len([l for l in text.splitlines() if l.startswith('REMARK 850 AGFR rebuilt')]) == 1, 'one line')
    rec = json.load(open(os.path.join(s.folder, 'session.json')))['saves'][-1]
    check(sorted(r['segment'] for r in rec['rebuilt']) == ['A', 'B', 'C'], rec['rebuilt'])
    renames = rec['reduce']['chain_renames']               # reduce's ids for the copies
    chains = ['A'] + sorted(c for c, v in renames.items() if v['chain'] == 'A')
    check(len(chains) == 3, renames)
    new = [((num, ''), crystal[num][0]) for num in range(first, last + 1)]
    segment = dict([('A', 'A')] + [(c, v['segment']) for c, v in renames.items()])
    for chain in chains:
        res = _check_loop(q, chain, (first - 1, ''), new, (last + 1, ''), '1H6W Biomol 1',
                          'copy %s' % segment[chain], _seconds(s, segment[chain]))
        if chain == 'A':
            names = sorted(set(crystal[first][1]) & set(res[(first, '')][1]))
            rmsd = numpy.sqrt(numpy.mean([((crystal[num][1][n] - res[(num, '')][1][n]) ** 2).sum()
                                          for num in range(first, last + 1) for n in ('N', 'CA', 'C', 'O')]))
            print('      1H6W A %d to %d: backbone RMSD to the deposited coordinates %.2f A'
                  % (first, last, rmsd))
            check(names, names)
    s.switch('ASU')
    d = s.rebuilt()
    check([r['segment'] for r in d['rebuilt']] == ['A'], d['rebuilt'])


# ---- the table ---------------------------------------------------------------
def test_zz_print_the_measured_table():
    """Every loop measured above, in one table (this test runs last)."""
    print('      entry          chain gap  n  C to N (A)        CA to CA (A)      contact (A)  time')
    for row in TABLE:
        print('      ' + row)
    check(TABLE or True)
