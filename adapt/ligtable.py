# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The ligand table: the chemical components each PDB entry holds, which tells
an entry's ligands apart from the other hetero groups (ADFR's PDBProfiler.getLigands).

ADFR carries a table compiled in March 2020 from Ligand Expo (cc-to-pdb.tdd:
126,335 entries, 29,761 codes); Ligand Expo was retired in February 2026.
ADapt ships a table built from the RCSB PDB Data API with the same content:
for every current entry, its non-polymer components, the monomers of its
branched entities (carbohydrates) and the non-standard residues of its
polymers, with each component's name and formula. It can be rebuilt on request
(`adapt ligands --update`, or from the page when an entry is not in it); a
rebuilt table is kept under its checksum in ~/.adapt/ligands/ ($ADAPT_STATE_DIR),
so a session is restored with the table it was made with.

File: gzipped text, '#' lines first (the build date and counts), then one line
per component: code TAB formula TAB name TAB the PDB ids (lower case, space
separated).
"""
import contextlib
import datetime
import gzip
import hashlib
import io
import json
import os
import re
import threading
import time
import urllib.request

SHIPPED = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ligands.tsv.gz')
ADFR_2020 = 'adfr-2020'          # ADFR's own table (sessions made before this one existed)
DATA_API = 'https://data.rcsb.org'
PDB_ID = re.compile(r'^([0-9][a-z0-9]{3}|pdb_[0-9]{4}[0-9a-z]{4})$')

_entries_query = """query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id
  nonpolymer_entities { rcsb_nonpolymer_entity_container_identifiers { nonpolymer_comp_id } }
  branched_entities { rcsb_branched_entity_container_identifiers { chem_comp_monomers } }
  polymer_entities { rcsb_polymer_entity_container_identifiers { chem_comp_nstd_monomers } } } }"""
_names_query = """query($ids: [String!]!) { chem_comps(comp_ids: $ids) { chem_comp { id name formula } } }"""


def user_dir():
    return os.path.join(os.path.expanduser(os.environ.get('ADAPT_STATE_DIR', '~/.adapt')), 'ligands')


def looks_like_pdb_id(name):
    return bool(PDB_ID.match((name or '').lower()))


class Table(object):
    """One ligand table file."""

    def __init__(self, path):
        self.path = path
        with open(path, 'rb') as fh:
            raw = fh.read()
        self.id = hashlib.sha256(raw).hexdigest()[:12]
        self.built = ''
        self.counts = {}
        self._cc = {}            # code -> (formula, name)
        self._pdb = {}           # pdb id -> [codes]
        for line in gzip.decompress(raw).decode('utf-8').splitlines():
            if line.startswith('#'):
                m = re.match(r'# built (\S+) from (\S+) \((.*)\)', line)
                if m:
                    self.built, self.source = m.group(1), m.group(2)
                    self.counts = dict((k, int(v)) for k, v in re.findall(r'(\w+) (\d+)', m.group(3)))
                continue
            code, formula, name, ids = line.split('\t')
            self._cc[code] = (formula, name)
            for p in ids.split():
                self._pdb.setdefault(p, []).append(code)

    def describe(self):
        return {'id': self.id, 'built': self.built, 'date': self.built[:10],
                'entries': len(self._pdb), 'codes': len(self._cc),
                'shipped': os.path.abspath(self.path) == os.path.abspath(SHIPPED)}

    def has(self, pdbid):
        return (pdbid or '').lower() in self._pdb

    def ligands(self, pdbid):
        return list(self._pdb.get((pdbid or '').lower(), []))

    def name(self, code):
        """{'name', 'formula'} as ADFR's chemComp holds them, or None."""
        f = self._cc.get(code)
        return {'name': f[1], 'formula': f[0]} if f else None


_tables = {}
_tables_lock = threading.Lock()


def _load(path):
    with _tables_lock:
        key = os.path.abspath(path)
        if key not in _tables or _tables[key][0] != os.path.getmtime(path):
            _tables[key] = (os.path.getmtime(path), Table(path))
        return _tables[key][1]


def available():
    """Every table on this machine: the shipped one and those rebuilt here."""
    out = [_load(SHIPPED)] if os.path.exists(SHIPPED) else []
    d = user_dir()
    if os.path.isdir(d):
        out += [_load(os.path.join(d, f)) for f in sorted(os.listdir(d)) if f.endswith('.tsv.gz')]
    return out


def current():
    """The most recently built table."""
    tabs = available()
    return max(tabs, key=lambda t: t.built) if tabs else None


def by_id(tid):
    """The table with this id, ADFR_2020 for ADFR's own, or None if it is not here."""
    if tid == ADFR_2020:
        return ADFR_2020
    return next((t for t in available() if t.id == tid), None)


# ---------------------------------------------------------------- in ADFR
_active = threading.local()
_installed = []


def install():
    """Route ADFR's residue classification through the table chosen with
    using(): during PDBProfiler._classify, the entry's ligand codes (renamed as
    the mmCIF conversion renamed them) and the components' names, which
    getLigands and getSites look up in ADFR's module-level pdb2lig and chemComp."""
    if _installed:
        return
    import ADFR.recFromPDB as R
    original = R.PDBProfiler._classify
    adfr_pdb2lig, adfr_chem, adfr_res = R.pdb2lig, R.chemComp, R.resData

    def _classify(self, pdbid, *args, **kw):
        table = getattr(_active, 'table', None)
        if table is None or table == ADFR_2020:
            return original(self, pdbid, *args, **kw)
        renamed = (getattr(self._receptor, 'conversion', None) or {}).get('residues') or {}
        codes = table.ligands(pdbid)
        names = dict(adfr_chem)
        for c in codes:
            n = table.name(c)
            if n:
                names[c] = names[renamed.get(c, c)] = n
        # a residue renamed by the mmCIF conversion (A1C5F -> L01) is checked for
        # missing atoms against its own component, not the one that owns the
        # new name in the Chemical Component Dictionary (L01 is a real code);
        # a component the dictionary ADFR carries does not know is not checked
        res = adfr_res
        if renamed:
            res = dict(adfr_res)
            for orig, new in renamed.items():
                if orig in adfr_res:
                    res[new] = adfr_res[orig]
                else:
                    res.pop(new, None)
        R.pdb2lig = {(pdbid or '').lower(): [renamed.get(c, c) for c in codes]}
        R.chemComp, R.resData = names, res
        try:
            return original(self, pdbid, *args, **kw)
        finally:
            R.pdb2lig, R.chemComp, R.resData = adfr_pdb2lig, adfr_chem, adfr_res

    R.PDBProfiler._classify = _classify
    _installed.append(original)


@contextlib.contextmanager
def using(table):
    """ADFR's ligand lookup uses `table` (a Table, or ADFR_2020) in this thread."""
    install()
    old = getattr(_active, 'table', None)
    _active.table = table
    try:
        yield
    finally:
        _active.table = old


# ---------------------------------------------------------------- building
def _post(query, ids, tries=5):
    body = json.dumps({'query': query, 'variables': {'ids': ids}}).encode()
    return _get(DATA_API + '/graphql', body, tries)


def _get(url, data=None, tries=5):
    # the RCSB is reached directly or through the proxy the environment names
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={
                'User-Agent': 'ADapt ligand table', 'Content-Type': 'application/json'})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


_building = threading.Lock()


def build(progress=None, threads=4, folder=None):
    """Rebuild the table from the RCSB PDB Data API (a few minutes; some 600
    requests) into `folder` (default user_dir()); the new Table."""
    from concurrent.futures import ThreadPoolExecutor
    if not _building.acquire(False):
        raise RuntimeError('the ligand table is already being rebuilt')
    try:
        say = progress or (lambda msg: None)
        ids = _get(DATA_API + '/rest/v1/holdings/current/entry_ids')
        say('%d entries in the PDB' % len(ids))
        batches = [ids[i:i + 500] for i in range(0, len(ids), 500)]
        cc2pdb = {}
        with ThreadPoolExecutor(threads) as pool:
            for k, d in enumerate(pool.map(lambda b: _post(_entries_query, b), batches)):
                for e in (d.get('data') or {}).get('entries') or []:
                    pid = e['rcsb_id'].lower()
                    for n in e.get('nonpolymer_entities') or []:
                        c = (n.get('rcsb_nonpolymer_entity_container_identifiers') or {}).get('nonpolymer_comp_id')
                        if c:
                            cc2pdb.setdefault(c, set()).add(pid)
                    for b in e.get('branched_entities') or []:
                        for c in (b.get('rcsb_branched_entity_container_identifiers') or {}).get('chem_comp_monomers') or []:
                            cc2pdb.setdefault(c, set()).add(pid)
                    for p in e.get('polymer_entities') or []:
                        for c in (p.get('rcsb_polymer_entity_container_identifiers') or {}).get('chem_comp_nstd_monomers') or []:
                            cc2pdb.setdefault(c, set()).add(pid)
                if (k + 1) % 100 == 0:
                    say('%d of %d batches of entries' % (k + 1, len(batches)))
        codes = sorted(cc2pdb)
        # names and formulas: those of the current table are kept (a component's
        # name rarely changes); only components new to it are asked for
        known = current()
        names = dict((c, known._cc[c]) for c in codes if known is not None and c in known._cc)
        missing = [c for c in codes if c not in names]
        cbatches = [missing[i:i + 1000] for i in range(0, len(missing), 1000)]
        with ThreadPoolExecutor(threads) as pool:
            for d in pool.map(lambda b: _post(_names_query, b), cbatches):
                for c in (d.get('data') or {}).get('chem_comps') or []:
                    cc = (c or {}).get('chem_comp') or {}
                    if cc.get('id'):
                        names[cc['id']] = (cc.get('formula') or '', cc.get('name') or '')
        say('%d components named' % len(names))
        return write(cc2pdb, names, len(ids), folder or user_dir())
    finally:
        _building.release()


def write(cc2pdb, names, n_entries, folder):
    """Write a table file into `folder` (named by its build date and checksum)."""
    clean = lambda s: ' '.join(str(s).split())
    built = datetime.datetime.now().astimezone().isoformat(timespec='seconds')
    with_lig = len({p for v in cc2pdb.values() for p in v})
    lines = ['# ADapt ligand table: the chemical components of each PDB entry',
             '# built %s from %s (entries %d with_ligands %d codes %d)'
             % (built, DATA_API, n_entries, with_lig, len(cc2pdb)),
             '# code\tformula\tname\tPDB ids']
    for c in sorted(cc2pdb):
        f, n = names.get(c, ('', ''))
        lines.append('%s\t%s\t%s\t%s' % (c, clean(f), clean(n), ' '.join(sorted(cc2pdb[c]))))
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode='wb', mtime=0) as gz:     # same text, same bytes
        gz.write(('\n'.join(lines) + '\n').encode('utf-8'))
    raw = buf.getvalue()
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, 'ligands-%s-%s.tsv.gz' % (built[:10], hashlib.sha256(raw).hexdigest()[:12]))
    tmp = path + '.part'
    with open(tmp, 'wb') as fh:
        fh.write(raw)
    os.replace(tmp, path)
    return _load(path)
