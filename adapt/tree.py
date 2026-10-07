# Author: Michel Sanner, PhD. Copyright (c) 2026. Derived from ADFR / AGFRgui code; see THIRD_PARTY_NOTICES.md.
"""The residue tree of one biomolecule tab, as data instead of Qt widgets.

A port of ADFR/GUI/pdbAnalyser.py:ReceptorItemsTree.buildUI (adcpsuite ADFR,
adfr-0.0.22). The same prody selections are made in the same order, so every
node carries exactly the atoms its QTreeWidgetItem._pmvObj carried in AGFRgui,
and the receptor/ligand flags set through it (flags.py) land on the same atoms.

Check states follow Qt5's QTreeWidgetItem rules, because the GUI relies on
them rather than setting parent states itself:
  * a node built with ItemIsTristate (Qt5: ItemIsAutoTristate) that has
    children DISPLAYS the state computed from its children
    (childrenCheckState): all checked -> checked, all unchecked -> unchecked,
    otherwise partial; and NO checkbox at all if any child has none;
  * setting the state of such a node pushes it down to every descendant that
    has a state in that column (QTreeWidgetItem::setData).
"""
from collections import OrderedDict

from ADFR.recFromPDB import getResidueKey, getResidueKeyAlt

CHECKED, PARTIAL, UNCHECKED = 'checked', 'partial', 'unchecked'
REC, LIG = 0, 1          # our two check columns (Qt columns 0 and 2)

# pdbAnalyser.py:residueNames, the targets offered by the mutate menu
RESIDUE_NAMES = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL']

# buildUI's category loop: tree label -> PDBProfiler category code
CATEGORIES = [
    ('protein', 'prote'), ('nucleic', 'nucle'),
    ('sidechains to complete', 'cmpsc'), ('missing atoms', 'missa'),
    ('cofactors', 'cofac'), ('modified', 'modif'),
    ('manual charges', 'manch'), ('ligands', 'ligan'),
    ('additives', 'addit'), ('water', '_wate'), ('other', 'other')]

DEFAULT_CHECKED = ('protein', 'nucleic', 'sidechains to complete')


class Node(object):
    """One QTreeWidgetItem.

    sel      : the atoms of _pmvObj (None when the item had none)
    gapsel   : _pmvGapObj, atoms a "missing residues" entry focuses on
    state    : explicit check state per column, None = no checkbox
    tristate : built with ItemIsTristate
    checkable: ItemIsUserCheckable (buildUI removes it for unsupported types)
    live     : clicking it changes flags (onClick needs _pmvObj + checkable;
               the whole "unknown atom types" subtree is inert here, see
               flags.py)
    """
    __slots__ = ('id', 'label', 'parent', 'children', 'sel', 'gapsel',
                 'state', 'tristate', 'checkable', 'live', 'kind', 'key',
                 'altloc', 'charge', 'rotamer', 'menu', 'expanded', 'resname',
                 'resnum', 'icode')

    def __init__(self, tree, label, parent=None, kind='item'):
        self.id = tree._newid()
        tree.nodes[self.id] = self
        self.label = label
        self.parent = parent
        self.children = []
        if parent is not None:
            parent.children.append(self)
        self.sel = self.gapsel = None
        self.state = [None, None]
        self.tristate = False
        self.checkable = True        # QTreeWidgetItem default flags
        self.live = True
        self.kind = kind
        self.key = None
        self.altloc = None           # {'alts','occ','current'} for altloc residues
        self.charge = None           # float for manual-charge residues
        self.rotamer = None          # rotamer.Preview for side chains to build
        self.menu = None             # 'sc' | 'mutate' | None
        self.expanded = False
        self.resname = self.resnum = self.icode = None

    # ---- Qt check-state semantics -------------------------------------
    def check_state(self, col):
        if self.children and self.tristate:
            seen = set()
            for c in self.children:
                s = c.check_state(col)
                if s is None:
                    return None
                if s == PARTIAL:
                    return PARTIAL
                seen.add(s)
                if len(seen) > 1:
                    return PARTIAL
            return seen.pop() if seen else None
        return self.state[col]

    def set_state(self, col, value):
        """QTreeWidgetItem::setData(col, CheckStateRole, value)."""
        self.state[col] = value
        if self.tristate and value != PARTIAL:
            for c in self.children:
                if c.check_state(col) is not None:
                    c.set_state(col, value)

    def leaves(self):
        if not self.children:
            yield self
        else:
            for c in self.children:
                for x in c.leaves():
                    yield x

    def to_json(self):
        d = {'id': self.id, 'label': self.label, 'kind': self.kind}
        if self.expanded:
            d['open'] = True
        if self.menu:
            d['menu'] = self.menu
        if self.altloc is not None:
            d['altloc'] = self.altloc
        if self.charge is not None:
            d['charge'] = self.charge
        if self.rotamer is not None:
            d['rotamer'] = self.rotamer.describe()
        if self.children:
            d['children'] = [c.to_json() for c in self.children]
        return d


class Tree(object):
    """ReceptorItemsTree for one entry of PDBProfiler._assemblies."""

    def __init__(self, session, name):
        self.session = session
        self.profile = session.profile
        self.name = name
        self._chaindInBM, self._mol = self.profile._assemblies[name]
        self._ag = self._mol._ag
        self.nodes = session.nodes
        self.roots = []
        self._typItems = {}
        self._altlocItems = {}       # {chid: {residue key: node}}
        self._hasUI = False

    def _newid(self):
        return self.session._newid()

    def header(self):
        return self._mol.name + ' receptor'

    # ------------------------------------------------------------------
    def build(self, rotamers=True):
        """buildUI (pdbAnalyser.py:485), minus the widgets."""
        if self._hasUI:
            return
        self._hasUI = True
        from . import rotamer as rotmod

        for chain in self._ag.select('segment A').getHierView().iterChains():
            chidSel = chid = chain.getChid()
            if chidSel == ' ':
                chidSel = '_'
            # only consider chains from ASU that appear in this BM
            if chid not in self._chaindInBM:
                continue
            root = chItem = Node(self, chid, kind='chain')
            root.expanded = True
            self.roots.append(root)
            self._typItems[chid] = typItems = {}

            altlocNO = chain.select('not altlocs and not unsre')
            altloc = chain.select('altlocs and not unsre')
            unsup = chain.select('unsre')

            for _iter, chainAtoms in enumerate([altlocNO, altloc, unsup]):
                if _iter == 1:
                    self._altlocItems[chid] = altlocItems = {}
                    if chainAtoms is None or len(chainAtoms) == 0:
                        continue
                    chItem = Node(self, 'altloc', root, kind='group')
                    chItem.expanded = True
                    typItems['altloc'] = chItem
                    chItem.tristate = True
                    chItem.sel = chainAtoms
                if _iter == 2:
                    if chainAtoms is None or len(chainAtoms) == 0:
                        continue
                    chItem = Node(self, 'unknown atom types', root, kind='group')
                    chItem.expanded = True
                    typItems['unsre'] = chItem
                    chItem.tristate = True
                if chainAtoms is None:
                    # buildUI would raise on None.select below; an empty set
                    # of a chain simply has no categories
                    continue

                for typ, cat in CATEGORIES:
                    atoms = chainAtoms.select('category %s' % cat)
                    if not atoms:
                        continue
                    parentItem = typItems[typ] = Node(self, typ, chItem,
                                                      kind='category')
                    parentItem.sel = atoms
                    parentItem.tristate = True
                    if typ in DEFAULT_CHECKED:
                        parentItem.state[REC] = CHECKED
                    elif not _iter == 2:
                        parentItem.state[REC] = UNCHECKED
                    # a category directly under the chain or under "altloc"
                    # is level 1 in openResMenu
                    if typ == 'sidechains to complete':
                        parentItem.menu = 'sc'
                    if _iter == 2:
                        parentItem.live = False

                    if _iter == 0 or _iter == 2:
                        self._plain_residues(atoms, parentItem, chidSel, typ,
                                             _iter, rotamers, rotmod)
                    else:
                        self._altloc_residues(atoms, parentItem, altlocItems)

            # entries for gaps, under the chain root (chItem at _iter 0)
            chgaps = self.profile._gapsPerChain.get(chid, None)
            if chgaps:
                self._gaps(root, chid, chgaps)
        return self

    def _plain_residues(self, atoms, parentItem, chidSel, typ, _iter,
                        rotamers, rotmod):
        profile = self.profile
        for res in atoms.getHierView().iterResidues():
            allres = self._ag.select('segment A chid %s resindex %d' % (
                chidSel, res.getResindices()[0]))
            if len(allres) > len(res):
                continue  # some atoms in this residue have altlocs
            resnum = res.getResnum()
            resname = res.getResname()
            icode = icodeSel = res.getIcode()
            if not icode:
                icodeSel = '_'
            item = Node(self, '', parentItem, kind='residue')
            item.sel = res
            item.key = key = getResidueKey(res)
            item.resname, item.resnum, item.icode = resname, resnum, icode
            if res.select('category modif'):
                old, new = profile._resToMutate.get(key, (None, None))
                if old is not None:
                    txt = "%s'%s -> %s" % (old, resnum, new)
                else:
                    txt = "%s'%s" % (resname, resnum)
            elif res.select('category cmpsc'):
                txt = "%s'%s%s build:%s" % (resname, resnum, icode,
                                            ' '.join(profile._missingAtoms[key]))
                if rotamers:
                    item.rotamer = rotmod.Preview.make(self, res, chidSel,
                                                       resname, resnum, icodeSel)
            elif res.select('category missa'):
                txt = "%s'%s%s missing:%s" % (resname, resnum, icode,
                                              ' '.join(profile._missingAtoms[key]))
            else:
                txt = "%s'%s%s" % (resname, resnum, icode)

            if res.select('receptor'):
                item.state[REC] = CHECKED
            elif not _iter == 2:
                item.state[REC] = UNCHECKED

            if res.select('category manch'):
                item.charge = float(res.getCharges()[0])

            item.label = txt
            if res.select('ligand'):
                item.state[LIG] = CHECKED
                item.state[REC] = UNCHECKED
            elif _iter == 2:
                item.checkable = False
            else:
                item.state[LIG] = UNCHECKED

            if _iter == 2:
                item.live = False
                item.expanded = True
                unsat = res.select('unsat')
                if unsat:
                    for atom in unsat:
                        a = Node(self, '%s (%s)' % (atom.getName(), atom.getElement()),
                                 item, kind='atom')
                        a.sel = _one_atom(self._ag, atom.getIndex())
                        a.checkable = False
                        a.live = False
            elif typ in ('sidechains to complete',):
                item.menu = 'sc'
            elif typ in ('modified', 'protein'):
                item.menu = 'mutate'

    def _altloc_residues(self, atoms, parentItem, altlocItems):
        altres = self.profile._altlocsRes
        typ = parentItem.label
        for res in atoms.getHierView().iterResidues():
            key = getResidueKey(res)
            item = Node(self, "%s'%s%s" % (res.getResname(), res.getResnum(),
                                           res.getIcode()),
                        parentItem, kind='residue')
            altlocItems[key] = item
            item.sel = res
            item.key = key
            item.resname, item.resnum, item.icode = (res.getResname(),
                                                     res.getResnum(), res.getIcode())
            item.state[LIG] = UNCHECKED
            item.state[REC] = CHECKED if res.select('receptor') else UNCHECKED
            ind, alts, occ = altres[key]
            item.altloc = {'alts': list(alts), 'occ': [float(o) for o in occ],
                           'current': int(ind)}
            if typ in ('modified', 'protein'):
                item.menu = 'mutate'
            elif typ == 'sidechains to complete':
                item.menu = 'sc'

    def _gaps(self, chItem, chid, chgaps):
        ag = self._ag

        def sel(num, icode):
            if icode:
                return ag.select('segment A chid %s resnum `%d` icode %s' % (chid, num, icode))
            return ag.select('segment A chid %s resnum `%d`' % (chid, num))

        misResItem = Node(self, 'missing residues', chItem, kind='gaps')
        if len(chgaps[0]):
            n = Node(self, 'N term: %s' % chgaps[0][1], misResItem, kind='gap')
            rName, rNum, rIcode = chgaps[0][0]
            n.gapsel = sel(rNum, rIcode)
        for gap in chgaps[1]:
            r1Name, r1Num, r1Icode = gap[0]
            r2Name, r2Num, r2Icode = gap[-1]
            rr = tuple(gap[1:-1])
            g = Node(self, 'gap: %s%d--' % (r1Name, r1Num) + '%s' * len(rr) % rr +
                     '--%s%d' % (r2Name, r2Num), misResItem, kind='gap')
            a, b = sel(r1Num, r1Icode), sel(r2Num, r2Icode)
            g.gapsel = a + b if (a is not None and b is not None) else (a or b)
        if len(chgaps[2]):
            c = Node(self, 'C term: %s' % chgaps[2][1], misResItem, kind='gap')
            cName, cNum, cIcode = chgaps[2][0]
            c.gapsel = sel(cNum, cIcode)
        for n in misResItem.children:
            n.checkable = False
            n.live = False
        misResItem.checkable = False
        misResItem.live = False

    # ------------------------------------------------------------------
    def to_json(self):
        return {'name': self.name, 'header': self.header(),
                'roots': [r.to_json() for r in self.roots]}

    def states(self):
        """{node id: [rec, lig]} for every node, computed as Qt displays them."""
        out = {}
        for r in self.roots:
            self._collect(r, out)
        return out

    def _collect(self, node, out):
        rec, lig = node.check_state(REC), node.check_state(LIG)
        out[node.id] = [rec, lig, node.live and node.checkable]
        for c in node.children:
            self._collect(c, out)

    def altloc_nodes(self):
        for chid, items in self._altlocItems.items():
            for key, node in items.items():
                yield key, node

    def manual_charge_nodes(self):
        for n in self.nodes.values():
            if n.charge is not None and _in_tree(n, self):
                yield n


def _one_atom(ag, index):
    from prody.atomic.selection import Selection
    return Selection(ag, [index], '')


def _in_tree(node, tree):
    while node.parent is not None:
        node = node.parent
    return node in tree.roots


def modres_parent(profile, res):
    """The standard residue a modified residue derives from, as its MODRES
    record gives it (resnum written with the insertion code, as prody keeps it)."""
    num = '%d%s' % (res.getResnum(), res.getIcode().strip())
    for poly in (profile._receptor.pdbHeader or {}).get('polymers', []):
        if poly.chid == res.getChid():
            for resname, resnum, base, _ in poly.modified or ():
                if resname == res.getResname() and resnum.strip() == num:
                    return base
    raise ValueError('no MODRES record for %s' % getResidueKeyAlt(res))


def mutate_options(tree, node):
    """The entries of openResMenu's mutate menu for a residue node."""
    profile = tree.profile
    res = node.sel
    resname = res.getResname()
    mutated = ' -> ' in node.label
    parent = node.parent.label
    if parent == 'modified':
        # AGFRgui reads profile._modifiedRes here, which only ADFR's commented-out
        # _classify2 ever built, so the menu raised AttributeError. The parent
        # residue is the one the file's MODRES record names for this residue.
        baseResName = modres_parent(profile, res)
        allKnown = not res.getFlags('unsre').any()
        if allKnown and mutated:
            names = [resname, 'Base Res: %s' % baseResName]
        else:
            names = ['Base Res: %s' % baseResName]
        names.extend(sorted(set(RESIDUE_NAMES) - set([baseResName])))
    else:
        if mutated:
            names = ['Base Res: %s' % resname]
            names.extend(sorted(set(RESIDUE_NAMES) - set([resname])))
        else:
            names = sorted(set(RESIDUE_NAMES) - set([resname]))
    current = profile._resToMutate.get(node.key, (None, None))[1]
    return {'names': names, 'current': current}
