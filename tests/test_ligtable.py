# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The ligand table: the one shipped, a rebuilt one taking over, sessions keeping
the table they were made with, and an entry released after 2020 (10HY) whose
ligand is found."""
import gzip
import json
import os

from common import check, outdir, DATA
from adapt import ligtable


def _state(name):
    """ADAPT_STATE_DIR under tests/_out; returns a function that restores it."""
    old = os.environ.get('ADAPT_STATE_DIR')
    os.environ['ADAPT_STATE_DIR'] = outdir(name)

    def restore():
        if old is None:
            os.environ.pop('ADAPT_STATE_DIR', None)
        else:
            os.environ['ADAPT_STATE_DIR'] = old
    return restore


def test_shipped_table():
    t = ligtable.Table(ligtable.SHIPPED)
    d = t.describe()
    check(d['entries'] > 200000 and d['codes'] > 50000, d)
    check(t.built.startswith('2026-'), t.built)
    check('ALY' in t.ligands('6BHG') and 'UNX' in t.ligands('6bhg'), t.ligands('6bhg'))
    check(t.ligands('10hy') == ['A1C5F'], t.ligands('10hy'))
    check(t.name('A1C5F')['formula'] == 'C17 H23 Cl N6 O2', t.name('A1C5F'))
    check(t.ligands('1fv1') == ['GOL', 'SO4'] and not t.has('9zzz'))
    check(ligtable.looks_like_pdb_id('10HY') and ligtable.looks_like_pdb_id('pdb_00001abc'))
    check(not ligtable.looks_like_pdb_id('myprotein') and not ligtable.looks_like_pdb_id('abcd'))


def test_rebuilt_table_takes_over_and_old_ones_stay():
    restore = _state('lig_state')
    try:
        d = ligtable.user_dir()
        for f in os.listdir(d) if os.path.isdir(d) else []:
            os.remove(os.path.join(d, f))
        shipped = ligtable.current()
        check(shipped.describe()['shipped'], 'without a rebuilt table: the shipped one')
        t = ligtable.write({'ZZZ': {'9zzz'}, 'GOL': {'9zzz', '1abc'}},
                           {'ZZZ': ('C1', 'TEST'), 'GOL': ('C3 H8 O3', 'GLYCEROL')}, 2, d)
        check(ligtable.current().id == t.id, 'the most recent table is used')
        check(ligtable.by_id(t.id).id == t.id and ligtable.by_id(shipped.id).id == shipped.id)
        check(ligtable.by_id(ligtable.ADFR_2020) == ligtable.ADFR_2020)
        check(ligtable.by_id('000000000000') is None)
        check(t.ligands('9ZZZ') == ['GOL', 'ZZZ'] and t.describe()['entries'] == 2)
    finally:
        restore()


def test_build_asks_only_for_new_names():
    """build() with the RCSB answered from here: entries, components, and names
    asked only for components the current table does not name."""
    restore = _state('lig_build')
    calls = []
    real_get, real_post = ligtable._get, ligtable._post

    def fake_get(url, data=None, tries=5):
        calls.append(url)
        return ['1FV1', '10HY', '9XYZ']

    def fake_post(query, ids, tries=5):
        calls.append(('post', tuple(ids)))
        if 'chem_comps' in query:
            return {'data': {'chem_comps': [{'chem_comp': {'id': c, 'name': 'NEW ' + c, 'formula': 'C9'}}
                                            for c in ids]}}
        ent = {'1FV1': {'nonpolymer_entities': [{'rcsb_nonpolymer_entity_container_identifiers': {'nonpolymer_comp_id': 'GOL'}}]},
               '10HY': {'nonpolymer_entities': [{'rcsb_nonpolymer_entity_container_identifiers': {'nonpolymer_comp_id': 'A1C5F'}}]},
               '9XYZ': {'branched_entities': [{'rcsb_branched_entity_container_identifiers': {'chem_comp_monomers': ['NAG', 'QQQQQ']}}],
                        'polymer_entities': [{'rcsb_polymer_entity_container_identifiers': {'chem_comp_nstd_monomers': ['MSE']}}]}}
        return {'data': {'entries': [dict(ent[i], rcsb_id=i) for i in ids]}}
    ligtable._get, ligtable._post = fake_get, fake_post
    try:
        t = ligtable.build()
        check(t.ligands('9xyz') == ['MSE', 'NAG', 'QQQQQ'] and t.ligands('10hy') == ['A1C5F'])
        asked = [c for c in calls if isinstance(c, tuple) and set(c[1]) <= {'QQQQQ', 'GOL', 'A1C5F', 'NAG', 'MSE'}
                 and 'QQQQQ' in c[1]]
        check(asked == [('post', ('QQQQQ',))], 'only the unknown component is named: %s' % asked)
        check(t.name('GOL')['name'] == 'GLYCEROL' and t.name('QQQQQ')['name'] == 'NEW QQQQQ')
        check(t.counts == {'entries': 3, 'with_ligands': 3, 'codes': 5}, t.counts)
    finally:
        ligtable._get, ligtable._post = real_get, real_post
        restore()


def _categories(session, tab=None):
    """{category label: [residue labels]} over every chain of a tab."""
    out = {}
    tree = session.tree(tab)

    def walk(n, cat=None):
        if n.kind == 'category':
            cat = n.label
        elif cat:
            out.setdefault(cat, []).append(n.label)
        for c in n.children:
            walk(c, cat)
    for r in tree.roots:
        walk(r)
    return out


def test_10hy_ligand_is_found_and_sessions_keep_their_table():
    """10HY (2026, mmCIF only): its ligand A1C5F, renamed L01 on conversion, is a
    ligand with the current table; with ADFR's 2020 table it was not (and the
    dictionary's own L01 made it look incomplete); each session records its
    table and is restored with it."""
    from adapt import fetch
    from adapt.session import Session
    from adapt import env
    try:
        path, _, _ = fetch.fetch('10HY', outdir('rcsb_10hy'))
    except fetch.FetchError as e:
        print('      skip: no route to the RCSB (%s)' % str(e)[:80])
        return
    data = open(path, 'rb').read()
    name = os.path.basename(path)
    new = Session('lig10hy_new0001', data, name, outdir('lig_sessions'), env.provenance())
    old = Session('lig10hy_old0001', data, name, outdir('lig_sessions'), env.provenance(),
                  ligand_table=ligtable.ADFR_2020)
    cn, co = _categories(new), _categories(old)
    check(any('L01' in r for r in cn.get('ligands', [])), cn)
    # with ADFR's table and checks, L01 is taken for the dictionary's own L01 and
    # listed as missing atoms
    check(not any('L01' in r for r in co.get('ligands', [])) and any('L01' in r for r in co.get('missing atoms', [])), co)
    check(new.record['ligand_table']['entry_listed'] is True, new.record['ligand_table'])
    check(old.record['ligand_table']['id'] == ligtable.ADFR_2020)
    for s in (new, old):
        s.flush()
        r = Session.restore(s.folder)
        check(r.record['ligand_table']['id'] == s.record['ligand_table']['id'])
        check(_categories(r) == _categories(s), 'restored with its own table')
    # a session recorded before ADapt had a table restores with ADFR's
    rec_path = os.path.join(old.folder, 'session.json')
    rec = json.load(open(rec_path))
    rec.pop('ligand_table')
    json.dump(rec, open(rec_path, 'w'))
    check(Session.restore(old.folder).ligand_table == ligtable.ADFR_2020)
