# Author: Aadhil Haq. Copyright (c) 2026.
"""Missing residues, rebuilt on request, gap by gap.

ADapt perceives missing residues (gaps.py) and by default builds nothing.
For a gap between two residues of a chain (an internal gap) the user can ask
for the missing residues to be built. A missing terminus is not offered: only
one of its ends is held, so where it goes would be a guess.

What is built:
  * the residues the file's REMARK 465 records list between the two residues
    either side of the gap, with REMARK 465's names and numbers (SEQRES gives
    one letter per residue, and ProDy writes some modified residues with the
    letter of their parent); a selenomethionine (MSE) is built as the
    methionine ADFR writes for every MSE, any other residue that is not one of
    the 20 amino acids means the gap is not offered;
  * heavy atoms only, by adapt/loops.py: placed by PDBFixer, then refined in
    Amber ff14SB with every observed atom held fixed, and refused when the
    loop they make is not sane. reduce protonates them with everything else.

A gap is not offered, and the tree says why, when REMARK 465 does not list
exactly its residues, a number it gives is used by the chain already, a
residue either side is not one of the 20 amino acids or lacks its backbone
atoms, the residue before it ends with an OXT, or its ends are further apart
than its residues can span (STEP per peptide bond).

Where it happens: prep.processed_bm() hands the molecule runAGFR.processBM()
made (side chains built, mutations made, alternate locations cleared) to
Rebuilder.apply(), which inserts the new residues in chain order. Then
processReceptor() protonates, charges and types them with everything else.
The new atoms carry the data and flags of the molecule's other atoms, are in
the receptor when the residues on both sides of the gap are, and are never in
the ligand.

What a build sees: the observed heavy atoms within reach of the gap, one
location each (_observed, _near), whether or not they will be written. Gaps
whose new residues could touch are built together, so they avoid each other;
a gap far from the others is built alone, so choosing another gap does not
change it. A copy of a chain in an assembly gets the gap rebuilt like the
chain itself, as ADFR completes side chains in every copy, each copy's
residues in their own place.

With nothing chosen nothing here runs, and the files are ADFR's.
"""
import hashlib
from collections import OrderedDict

import numpy

from . import gaps as gapmod
from . import loops

# How far new residues can reach, in A: a missing residue's CA is at most one
# trans peptide (3.8 A) per residue from the end of the gap it hangs from, so
# every new atom lies within STEP * (n + 1) / 2 + SIDE of one of the two ends
# (SIDE: CA to the far end of an Arg or Lys side chain). Observed atoms a
# MARGIN beyond that are kept for the build: PDBFixer's repulsion has fallen
# to 0.2% of its contact value at 10 A, and loops.refine() cuts its
# Lennard-Jones terms off at 9 A.
STEP = 3.8
SIDE = 7.5
MARGIN = 10.0
CONTACT = 4.0           # A: new residues of two gaps this close are built together

# the residues PDBFixer and ff14SB have templates for, and their SEQRES letters
STANDARD = OrderedDict([
    ('ALA', 'A'), ('ARG', 'R'), ('ASN', 'N'), ('ASP', 'D'), ('CYS', 'C'),
    ('GLN', 'Q'), ('GLU', 'E'), ('GLY', 'G'), ('HIS', 'H'), ('ILE', 'I'),
    ('LEU', 'L'), ('LYS', 'K'), ('MET', 'M'), ('PHE', 'F'), ('PRO', 'P'),
    ('SER', 'S'), ('THR', 'T'), ('TRP', 'W'), ('TYR', 'Y'), ('VAL', 'V')])
# a missing selenomethionine is built as the methionine ADFR writes for one
# (runAGFR.processBM turns every MSE into MET)
AS_WRITTEN = {'MSE': 'MET'}


class NotBuilt(Exception):
    """A chosen gap could not be built; .key names it, str() says why."""

    def __init__(self, key, why):
        Exception.__init__(self, why)
        self.key = key


def _order(r):
    """(resnum, icode): residues of a chain sort this way ('' before 'A')."""
    return (int(r[1]), r[2] or '')


def _key(chid, r):
    return '%s:%s:%s:%s' % (chid, r[0], r[1], r[2] or '_')


def _sel(chid, r):
    """ProDy selection of residue r = (resname, resnum, icode) of chain chid."""
    return 'chid %s and resnum `%d` and icode %s' % (chid if chid != ' ' else '_', int(r[1]),
                                                     r[2] or '_')


class Gap(object):
    """An internal gap of a chain, as gaps.py perceives it, and the residues
    that rebuilding it adds.

    before, after: (resname, resnum, icode) of the residues on either side
    seq          : the missing residues, one SEQRES letter each
    residues     : [(resname, resnum, icode)] to build, in chain order
    why          : why it cannot be rebuilt, or None
    """
    __slots__ = ('chid', 'before', 'after', 'seq', 'residues', 'why')

    def __init__(self, chid, before, after, seq):
        self.chid = chid
        self.before = tuple(before)
        self.after = tuple(after)
        self.seq = seq
        self.residues = None
        self.why = None

    @property
    def key(self):
        return gapmod.key(self.chid, self.before, self.after)

    def keys(self):
        """The residue keys of the new residues, as the REMARKs write them."""
        return [_key(self.chid, r) for r in self.residues]


def perceive(profile):
    """OrderedDict {key: Gap} of the internal gaps gaps.py found, chain by chain."""
    ag = profile._receptor._ag
    listed = _remark465(profile)
    out = OrderedDict()
    for chid, chgaps in profile._gapsPerChain.items():
        used = _numbers_in_use(ag, chid)
        for g in chgaps[1]:
            gap = Gap(chid, g[0], g[-1], ''.join(g[1:-1]))
            try:
                gap.residues, gap.why = _residues(gap, listed.get(chid, []), used)
                gap.why = gap.why or _ends(ag, gap)
            except Exception as e:        # a file ADapt opens must not fail here
                gap.residues = None
                gap.why = 'could not be examined (%s: %s)' % (type(e).__name__, e)
            out[gap.key] = gap
    return out


def _remark465(profile):
    """{chain: [(resname, resnum, icode)]}: the file's REMARK 465 records, in
    file order (those of the first model, which is the one ADapt reads)."""
    out = {}
    for chid, resname, resnum, icode, model in \
            (profile._receptor.pdbHeader or {}).get('missing_residues', []):
        if model in ('', '1', None):
            out.setdefault(chid, []).append((resname, int(resnum), icode or ''))
    return out


def _numbers_in_use(ag, chid):
    """(resnum, icode) of every residue the chain id has in the file, waters and
    hetero groups included: a new residue must not share one."""
    sel = ag.select('chid %s' % (chid if chid != ' ' else '_'))
    if sel is None:
        return set()
    return set(zip(sel.getResnums().tolist(), sel.getIcodes().tolist()))


def _residues(gap, listed, used):
    """([(resname, resnum, icode)], None) for the residues of the gap, or
    (None, why): REMARK 465's records between the two ends of the gap, which
    must be as many as the gap's SEQRES letters and agree with them."""
    lo, hi = _order(gap.before), _order(gap.after)
    run = [r for r in listed if lo < _order(r) < hi]
    if len(run) != len(gap.seq):
        return None, ('REMARK 465 does not list its %d residues, so their names and numbers '
                      'are not known' % len(gap.seq))
    for r, letter in zip(run, gap.seq):
        if AS_WRITTEN.get(r[0], r[0]) not in STANDARD:
            return None, '%s %d%s is not one of the 20 amino acids' % r
        if r[0] in STANDARD and STANDARD[r[0]] != letter:
            return None, 'REMARK 465 lists %s %d%s where SEQRES has %s' % (r + (letter,))
        if (r[1], r[2]) in used:
            return None, 'REMARK 465 numbers %s %d%s, which the chain already uses' % r
    return [(AS_WRITTEN.get(r[0], r[0]), r[1], r[2]) for r in run], None


def _ends(ag, gap):
    """Why the residues on either side cannot hold new ones, or None: they must
    be amino acids with their backbone, the residue before must not end with an
    OXT, and their CA must be close enough for the gap's residues to span."""
    at = {}
    for res, names in ((gap.before, ('CA', 'C')), (gap.after, ('N', 'CA'))):
        if AS_WRITTEN.get(res[0], res[0]) not in STANDARD:
            return '%s %d%s, beside it, is not one of the 20 amino acids' % res
        for name in names:
            atom = ag.select('%s and name %s' % (_sel(gap.chid, res), name))
            if atom is None:
                return '%s %d%s has no %s atom' % (res + (name,))
            at[res, name] = atom.getCoords()[0]
    if ag.select('%s and name OXT' % _sel(gap.chid, gap.before)) is not None:
        return '%s %d%s ends with an OXT atom' % gap.before
    d = numpy.linalg.norm(at[gap.before, 'CA'] - at[gap.after, 'CA'])
    span = STEP * (len(gap.seq) + 1)
    if d > span:
        return ('its ends are %.1f A apart (CA to CA); %d residues span at most %.1f A'
                % (d, len(gap.seq), span))
    return None


class Site(object):
    """A chosen gap in one copy (segment) of its chain, in a processed molecule:
    the atoms the new residues join (index of the C before, of the N after),
    and whether those two are in the receptor."""
    __slots__ = ('segment', 'gap', 'c', 'n', 'ends', 'receptor')

    @property
    def names(self):
        return [r[0] for r in self.gap.residues]

    @property
    def labels(self):
        """'SER -3', ...: the residue before the gap, the new ones, the one after."""
        return ['%s %d%s' % r for r in [self.gap.before] + self.gap.residues + [self.gap.after]]

    def reach(self):
        return STEP * (len(self.gap.seq) + 1) / 2. + SIDE

    def spec(self):
        return (self.segment, self.gap.chid, self.gap.before, self.gap.after, tuple(self.names))

    def where(self, chains):
        """(chain, residue) indices of the residue after the gap in _pdb's text."""
        gap = self.gap
        after = (self.segment, gap.chid, int(gap.after[1]), gap.after[2] or '')
        before = (self.segment, gap.chid, int(gap.before[1]), gap.before[2] or '')
        for k, residues in enumerate(chains):
            if after in residues:
                j = residues.index(after)
                if j > 0 and residues[j - 1] == before:
                    return k, j
        raise RuntimeError('%s does not follow %s in the molecule' % (
            _key(gap.chid, gap.after), _key(gap.chid, gap.before)))


class Rebuilder(object):
    """The internal gaps of one PDB file, the ones the user chose to rebuild,
    and the builds made so far (by everything that goes into one, so a save
    writes the coordinates the preview showed)."""

    def __init__(self, profile):
        self.profile = profile
        self.gaps = perceive(profile)
        self.chosen = set()
        self._builds = {}
        if not loops.available():
            for gap in self.gaps.values():
                gap.why = gap.why or 'PDBFixer and OpenMM are not installed in this Python'

    def choose(self, key, on):
        gap = self.gaps[key]
        if on and gap.why:
            raise ValueError(gap.why)
        if on:
            self.chosen.add(key)
        else:
            self.chosen.discard(key)

    def refuse(self, key, why):
        """A chosen gap whose loop was not sane: no longer chosen, and why."""
        self.chosen.discard(key)
        self.gaps[key].why = why

    def chosen_in(self, chids):
        """The chosen gaps of these chains, in chain order."""
        return [g for k, g in self.gaps.items() if k in self.chosen and g.chid in chids]

    def apply(self, mol, bm):
        """processBM's molecule `mol` (made from the tab's molecule `bm`) with
        the chosen gaps of its chains rebuilt in every copy of them: a new
        molecule, or `mol` itself when no chosen gap is in it."""
        sites = self._sites(mol)
        if not sites:
            return mol
        keep = _observed(mol, bm, sites)
        done = []
        for group in _groups(sites):
            done.extend(zip(group, self._build(mol._ag, keep, group)))
        return insert(mol, done)

    def _sites(self, mol):
        ag = mol._ag
        segments = list(OrderedDict.fromkeys(ag.getSegnames()))
        out = []
        for gap in self.chosen_in(set(ag.getChids())):
            for seg in segments:
                c = _join_atom(mol, seg, gap.chid, gap.before, 'C')
                n = _join_atom(mol, seg, gap.chid, gap.after, 'N')
                if c is None or n is None:
                    continue        # this copy lacks the chain or an end of the gap
                s = Site()
                s.segment, s.gap, s.c, s.n = seg, gap, c, n
                s.ends = ag.getCoords()[[c, n]]
                s.receptor = bool(ag.getFlags('receptor')[c] and ag.getFlags('receptor')[n])
                out.append(s)
        return out

    def _build(self, ag, keep, group):
        """[loops.Loop] for the sites of `group`, from one build (or from the
        build that had the same input)."""
        text, chains = _pdb(ag, keep & _near(ag, group))
        spec = [s.spec() for s in group]
        key = '%r %r %s %s %r\n' % (spec, loops.SEEDS, loops.PLATFORM, loops.METHOD,
                                    loops.versions())
        digest = hashlib.sha1((key + text).encode()).hexdigest()
        if digest not in self._builds:
            try:
                built = loops.build(text, chains, group)
            except loops.NotBuilt as e:
                raise NotBuilt(e.site.gap.key, str(e))
            fixer = loops.versions()[0]
            for s, b in zip(group, built):
                print('REBUILD gap /%s/%s %s as %s to %s (PDBFixer %s and %s, seed %d, %.1f s)' % (
                    s.segment, s.gap.key, s.gap.seq, s.gap.keys()[0], s.gap.keys()[-1],
                    fixer, loops.METHOD, b.seed, b.seconds['total']))
            self._builds[digest] = (built, digest)
        return [(b, self._builds[digest][1]) for b in self._builds[digest][0]]


def _join_atom(mol, seg, chid, res, name):
    """Index of atom `name` of residue `res` in copy `seg`; of the location
    written when the atom has two (processBM keeps both, unflagged)."""
    sel = mol.select('segment %s and %s and name %s' % (seg, _sel(chid, res), name))
    if sel is None:
        return None
    ind = sel.getIndices()
    written = sel.getFlags('receptor') | sel.getFlags('ligand')
    return int(ind[written][0] if written.any() else ind[0])


def _observed(mol, bm, sites):
    """Mask of the processed molecule's atoms a build sees: the heavy atoms that
    are not deleted, at every location the file gives them. A loop then avoids
    whichever location is written, and does not change when another is chosen
    (ADFR breaks a tie of occupancies by set order, which differs from run to
    run). The residues either side of a gap are seen as written, since the new
    residues join them. processBM cleared the altloc letters; they are read
    from the tab's molecule, whose atoms come first and in the same order."""
    ag = mol._ag
    keep = ~ag.getFlags('deleted') & ~_hydrogen(ag)
    alt = numpy.full(len(ag), ' ')
    alt[:len(bm._ag)] = bm._ag.getAltlocs()
    written = ag.getFlags('receptor') | ag.getFlags('ligand')
    for s in sites:
        for res in (s.gap.before, s.gap.after):
            ind = mol.select('segment %s and %s' % (s.segment, _sel(s.gap.chid, res))).getIndices()
            loc = alt[ind]
            if (loc == ' ').all():
                continue
            mine = written[ind] & (loc != ' ')
            letter = loc[mine][0] if mine.any() else min(set(loc) - set(' '))
            keep[ind[(loc != ' ') & (loc != letter)]] = False
    return keep


def _hydrogen(ag):
    h = ag.getFlags('hydrogen')
    return h if h is not None else numpy.zeros(len(ag), bool)


def _groups(sites):
    """The sites, in groups of those whose new residues could touch."""
    group = list(range(len(sites)))

    def top(i):
        while group[i] != i:
            i = group[i]
        return i
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            d = numpy.linalg.norm(sites[i].ends[:, None] - sites[j].ends[None], axis=2).min()
            if d < sites[i].reach() + sites[j].reach() + CONTACT:
                group[top(j)] = top(i)
    out = OrderedDict()
    for i, s in enumerate(sites):
        out.setdefault(top(i), []).append(s)
    return list(out.values())


def _residue_ids(ag):
    """One integer per atom, the same for the atoms of one residue (by segment,
    chain, number and insertion code), counted in order of appearance."""
    seen = {}
    return numpy.array([seen.setdefault(k, len(seen)) for k in zip(
        ag.getSegnames(), ag.getChids(), ag.getResnums(), ag.getIcodes())])


def _near(ag, group):
    """Mask of the atoms of the residues that come within reach + MARGIN of an
    end of a gap of the group."""
    xyz = ag.getCoords()
    near = numpy.zeros(len(ag), bool)
    for s in group:
        d = numpy.linalg.norm(xyz[:, None] - s.ends[None], axis=2).min(1)
        near |= d <= s.reach() + MARGIN
    rid = _residue_ids(ag)
    return numpy.isin(rid, rid[near])


def _pdb(ag, mask):
    """PDB text of the atoms in `mask` for PDBFixer, and chains[k] = [(segment,
    chid, resnum, icode)], the residues of its chain k in order. One chain per
    segment and chain id, each residue's atoms together (ADFR appends the side
    chains it builds), residues numbered 1, 2, ... (PDBFixer goes by their
    order, and a file's own numbers may not fit a PDB line). An atom's second
    location goes into a last chain of one-atom residues, which PDBFixer keeps
    as obstacles (it has no template for them)."""
    rid = _residue_ids(ag)
    seg, chid = ag.getSegnames(), ag.getChids()
    res = list(zip(seg, chid, ag.getResnums().tolist(), ag.getIcodes()))
    chains = OrderedDict()                       # (seg, chid) -> OrderedDict(rid -> [indices])
    for i in numpy.flatnonzero(mask):
        chains.setdefault((seg[i], chid[i]), OrderedDict()).setdefault(rid[i], []).append(i)
    names, resnames, elements = ag.getNames(), ag.getResnames(), ag.getElements()
    hetero, xyz = ag.getFlags('hetatm'), ag.getCoords()
    out, order, extra = [], [], []

    def line(i, hetatm, name, resname, c, num):
        name = name if len(name) == 4 or len(elements[i]) == 2 else ' ' + name
        out.append('%-6s%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s' % (
            'HETATM' if hetatm else 'ATOM', (len(out) + 1) % 100000, name, resname, c[:1], num,
            xyz[i][0], xyz[i][1], xyz[i][2], elements[i]))
    for (s, c), residues in chains.items():
        order.append([])
        for k, (r, atoms) in enumerate(residues.items()):
            order[-1].append(tuple(res[atoms[0]]))
            seen = set()
            for i in atoms:
                if names[i] in seen:
                    extra.append(i)
                    continue
                seen.add(names[i])
                line(i, hetero[i], names[i], resnames[i], c, k + 1)
        out.append('TER')
    if extra:
        order.append([])
        for k, i in enumerate(extra):
            order[-1].append(('', '', k + 1, ''))
            line(i, True, elements[i], 'ALT', 'Z', k + 1)
        out.append('TER')
    out.append('END')
    return '\n'.join(out) + '\n', order


# ---- the new atoms in ADFR's molecule ------------------------------------
def insert(mol, done):
    """A new Molecule: `mol`'s atoms that are not deleted, in order, with the
    new residues of each (site, (loop, digest)) of `done` placed before the
    residue after their gap, carrying the data and flags of the molecule's
    atoms. Its .rebuilt lists what was inserted (record())."""
    from MolKit2.molecule import Molecule
    ag = mol._ag
    pieces = []
    for site, (loop, digest) in done:
        gap = site.gap
        after = mol.select('segment %s and %s' % (site.segment, _sel(gap.chid, gap.after)))
        pieces.append((int(after.getIndices().min()), site, loop, _piece(ag, site, loop), digest))
    pieces.sort(key=lambda p: p[0])
    new, spans = _merged(ag, pieces)
    out = Molecule(mol.name, new)
    out.rebuilt = [record(site, loop, span, digest)
                   for (at, site, loop, piece, digest), span in zip(pieces, spans)]
    return out


def _piece(ag, site, loop):
    """An AtomGroup of the new atoms of one site, with every data label and
    flag of `ag` set as for the molecule's own protein atoms."""
    from prody import AtomGroup
    from prody.atomic.functions import SAVE_SKIP_POINTER as SKIP
    gap, n = site.gap, len(loop.atoms)
    p = AtomGroup('rebuilt')
    p.setCoords(loop.coords)
    values = {'name': [name for k, name, el in loop.atoms],
              'resname': [gap.residues[k][0] for k, name, el in loop.atoms],
              'resnum': [gap.residues[k][1] for k, name, el in loop.atoms],
              'icode': [gap.residues[k][2] for k, name, el in loop.atoms],
              'element': [el for k, name, el in loop.atoms],
              'chain': [gap.chid] * n, 'segment': [site.segment] * n,
              'altloc': [' '] * n, 'occupancy': [1.0] * n, 'beta': [0.0] * n,
              'serial': [0] * n, 'charge': [0.0] * n, 'category': ['prote'] * n,
              'secondary': ['C'] * n}
    for label in ag.getDataLabels():
        if label in SKIP:
            continue
        old = ag.getData(label)
        v = values.get(label)
        # other labels (atomicNumber, radius) are made again by Molecule()
        p.setData(label, numpy.array(v, old.dtype) if v is not None
                  else numpy.zeros((n,) + old.shape[1:], old.dtype))
    for label in ag.getFlagLabels():
        computed = p.getFlags(label)        # protein, bb, sc, ...: from the names
        if computed is None:
            computed = numpy.zeros(n, bool)
            if label == 'receptor':
                computed[:] = site.receptor
        p._setFlags(label, numpy.array(computed, bool))
    p._setFlags('rebuilt', numpy.ones(n, bool))
    return p


def _merged(ag, pieces):
    """The AtomGroup of `ag`'s atoms that are not deleted with the pieces
    inserted, and the (start, end) of each piece in it."""
    from prody import AtomGroup
    from prody.atomic.functions import SAVE_SKIP_POINTER as SKIP
    keep = numpy.flatnonzero(~ag.getFlags('deleted'))
    parts, where, spans = [], numpy.full(len(ag), -1), []
    start = n = 0
    for at, site, loop, piece, digest in pieces:
        cut = int(numpy.searchsorted(keep, at))
        parts.append(keep[start:cut])
        where[keep[start:cut]] = numpy.arange(n, n + cut - start)
        n += cut - start
        parts.append(piece)
        spans.append((n, n + len(piece)))
        n += len(piece)
        start = cut
    parts.append(keep[start:])
    where[keep[start:]] = numpy.arange(n, n + len(keep) - start)

    def column(get):
        return numpy.concatenate([get(x) if isinstance(x, AtomGroup) else get(ag)[x]
                                  for x in parts])
    new = AtomGroup(ag.getTitle())
    new.setCoords(column(lambda a: a.getCoords()))
    for label in ag.getDataLabels():
        if label not in SKIP:
            new.setData(label, column(lambda a: a.getData(label)))
    for label in ag.getFlagLabels() + ['rebuilt']:
        flags = ag.getFlags(label)
        new._setFlags(label, numpy.concatenate([
            x.getFlags(label) if isinstance(x, AtomGroup)
            else (flags[x] if flags is not None else numpy.zeros(len(x), bool)) for x in parts]))
    # bonds: the molecule's own among the atoms kept, each piece's, and the
    # peptide bonds that join a piece to the residues on either side
    bonds = []
    if ag._bonds is not None:
        mapped = where[ag._bonds]
        bonds.extend(mapped[(mapped >= 0).all(1)].tolist())
    for (at, site, loop, piece, digest), (s, e) in zip(pieces, spans):
        bonds.extend([s + i, s + j] for i, j in loop.bonds)
        bonds.append([int(where[site.c]), s + _atom(loop, 0, 'N')])
        bonds.append([s + _atom(loop, len(site.gap.residues) - 1, 'C'), int(where[site.n])])
    new.setBonds(numpy.array(bonds, int))
    return new, spans


def _atom(loop, k, name):
    return next(i for i, (r, n, el) in enumerate(loop.atoms) if r == k and n == name)


def record(site, loop, span, digest):
    """What the session records about one rebuilt site."""
    fixer, openmm = loops.versions()
    gap = site.gap
    return {'gap': gap.key, 'segment': site.segment, 'sequence': gap.seq,
            'residues': gap.keys(), 'written': site.receptor, 'atoms': list(span),
            'method': 'PDBFixer %s and %s' % (fixer, loops.METHOD), 'openmm': openmm,
            'platform': loops.PLATFORM, 'seed': loop.seed, 'geometry': loop.geometry,
            'seconds': loop.seconds, 'build': digest}


def remarks(mol):
    """REMARK 850 lines for the gaps rebuilt in the receptor, once per gap
    (its copies in an assembly are rebuilt with it, as mutations are)."""
    out = []
    for r in getattr(mol, 'rebuilt', ()):
        if not r['written']:
            continue
        line = 'REMARK 850 AGFR rebuilt %s to %s (%s) using %s, seed %d' % (
            r['residues'][0], r['residues'][-1], r['sequence'], r['method'], r['seed'])
        if line not in out:
            out.append(line)
    return out


def preview(mol):
    """What the page draws for the rebuilt residues of a processed molecule:
    one entry per site, with its PDB text and the two joins."""
    import io
    import prody
    out = []
    xyz = mol._ag.getCoords()
    for r in getattr(mol, 'rebuilt', ()):
        s, e = r['atoms']
        text = io.StringIO()
        prody.writePDBStream(text, mol._ag[s:e])
        joins = [[xyz[i].tolist(), xyz[j].tolist()] for i, j in mol._ag._bonds.tolist()
                 if (s <= i < e) != (s <= j < e)]
        out.append({'gap': r['gap'], 'segment': r['segment'], 'written': r['written'],
                    'residues': r['residues'], 'sequence': r['sequence'],
                    'pdb': text.getvalue(), 'joins': joins, 'geometry': r['geometry']})
    return out
