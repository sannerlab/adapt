# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The receptor as reduce protonated it: <pdbid>_rec_H.pdb next to the PDBQT,
and reduce's decisions in the session record (adapt/reduce_log.py)."""
import json
import os

from common import check, read, session


def _save_receptor(s):
    d = s.save('receptor')
    folder = os.path.join(s.folder, 'save_%03d_rec' % d['save'])
    names = [f['name'] for f in d['files']]
    return d, [os.path.join(folder, n) for n in names]


def _atoms(lines):
    for l in lines:
        if l.startswith(('ATOM', 'HETATM')):
            yield l


def _xyz(l):
    return tuple(round(float(l[i:i + 8]), 3) for i in (30, 38, 46))


def _is_h_pdbqt(l):
    return l[77:79].strip() in ('H', 'HD', 'HS')


def _is_h_pdb(l):
    el = l[76:78].strip()
    return el == 'H' if el else l[12:16].strip().lstrip('0123456789').startswith('H')


def test_protonated_receptor_pdb():
    s = session('6ggd.pdb', 'reduce_pdb')
    d, (pdbqt, pdb) = _save_receptor(s)
    check(os.path.basename(pdb) == '6ggd_rec_H.pdb', pdb)
    q = read(pdbqt).decode().splitlines()
    p = read(pdb).decode().splitlines()
    # the same provenance header, then reduce's own output
    head = [l for l in q if l.startswith('REMARK 850 AGFR')]
    check(p[:len(head)] == head, 'same REMARK 850 header')
    check(head[-1] == 'REMARK 850 AGFR end' and p[len(head)].startswith('USER  MOD reduce'),
          p[len(head) - 1:len(head) + 1])
    # the same heavy atoms; every polar H of the PDBQT, and the non-polar ones
    heavy_q = [_xyz(l) for l in _atoms(q) if not _is_h_pdbqt(l)]
    heavy_p = [_xyz(l) for l in _atoms(p) if not _is_h_pdb(l)]
    check(sorted(heavy_q) == sorted(heavy_p), (len(heavy_q), len(heavy_p)))
    h_q = set(_xyz(l) for l in _atoms(q) if _is_h_pdbqt(l))
    h_p = set(_xyz(l) for l in _atoms(p) if _is_h_pdb(l))
    check(h_q and h_q <= h_p, 'every polar H of the PDBQT is in the PDB (%d of %d)'
          % (len(h_q & h_p), len(h_q)))
    check(len(h_p) > len(h_q), 'non-polar hydrogens too (%d > %d)' % (len(h_p), len(h_q)))
    # reduce's decisions, recorded with the save
    rec = json.load(open(os.path.join(s.folder, 'session.json')))
    red = rec['saves'][-1]['reduce']
    check(red['version'].startswith('reduce.'), red['version'])
    check(red['chain_renames'] == {}, 'one copy: no chain renamed')
    dec = [x for x in red['decisions'] if 'residue' in x]
    check(dec, 'His/Asn/Gln decisions recorded')
    for x in dec:
        chain, name, num, icode = x['residue'].split(':')
        check(chain == 'A' and name in ('HIS', 'ASN', 'GLN'), x)      # Biomol 1 is chain A
        check(x['category'] in ('K', 'C', 'X', 'F', None) and x['decision'], x)
    his = [x for x in dec if ':HIS:' in x['residue']]
    check(his and all(x['decision'].endswith(('HD1', 'HE2', 'HN')) or 'FLIP' in x['decision']
                      for x in his), his)
    # every recorded decision is one of reduce's USER MOD lines in the PDB
    user = [l for l in p if l.startswith('USER  MOD')]
    for x in dec:
        chain, name, num, icode = x['residue'].split(':')
        check(any(' %s %4s%s%s ' % (chain, num, '' if icode == '_' else icode, name) in l or
                  ' %s%4s%s%s ' % (chain, num, ' ' if icode == '_' else icode, name) in l
                  for l in user), x)


def test_symmetry_copies_are_renamed_for_reduce():
    """1H6W Biomol 1: chains A B, three copies. reduce needs unique chain ids,
    so the copies are renamed; the record says how."""
    s = session('1h6w.pdb', 'reduce_sym')
    s.switch('Biomol 1')
    d, (pdbqt, pdb) = _save_receptor(s)
    rec = json.load(open(os.path.join(s.folder, 'session.json')))
    ren = rec['saves'][-1]['reduce']['chain_renames']
    check(len(ren) == 4, ren)       # 3 copies x 2 chains, the first copy keeps A and B
    check(sorted(set(v['chain'] for v in ren.values())) == ['A', 'B'], ren)
    check(sorted(set(v['segment'] for v in ren.values())) == ['B', 'C'], ren)
    chains_pdb = set(l[21] for l in _atoms(read(pdb).decode().splitlines()))
    check(set(ren) <= chains_pdb, (ren, chains_pdb))


def test_protonate_copy_is_verbatim():
    """adapt/reduce_log.py:protonate is ADFR's, plus only the lines marked ADAPT."""
    import inspect
    from adapt import reduce_log
    ours = inspect.getsource(reduce_log.protonate).splitlines()
    theirs = inspect.getsource(reduce_log._ORIGINAL).splitlines()
    i = next(k for k, l in enumerate(ours) if '# ADAPT:' in l)
    added = ours[i:i + 4]
    check(added[2].strip() == 'LAST.clear()' and added[3].strip().startswith('LAST.update('), added)
    check(ours[:i] + ours[i + 4:] == theirs, 'only the marked lines differ')
