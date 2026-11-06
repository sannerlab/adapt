# Author: Michel Sanner, PhD. Copyright (c) 2026. Derived from ADFR / AGFRgui code; see THIRD_PARTY_NOTICES.md.
"""What a click in the tree does to the receptor/ligand flags.

Ports of ReceptorItemsTree.setRecFlag, setLigFlag, _setItemFlags, onClick,
AltlocsChooser.onclick, _setSCcompletion and _selected (pdbAnalyser.py,
adcpsuite ADFR). They act on the same prody atoms with the same calls, so the
atoms runAGFR.processReceptor / processLigands later select as 'receptor' /
'ligand' are the ones AGFRgui would have selected.

Deliberate differences from the Qt code, each a GUI defect rather than a
behaviour anyone depends on:
  * An altloc chosen in one tab is applied to the atoms of a tab opened for
    the first time afterwards. AGFRgui showed the choice there (and wrote it
    in the REMARK 850 altlocs line) but left the atoms on the file's default.
  * A save right after a tab switch writes the new tab's biomolecule;
    AGFRgui kept the previous tab's processed molecule (session.switch).
  * Clicking a node's LABEL only focuses the viewer. In AGFRgui a label click
    re-applies the checkbox state, and on a PARTIAL category that state (1) is
    truthy, so every residue under it silently became receptor while the
    checkbox still showed partial.
  * Mutations and side-chain completion choices invalidate the processed
    biomolecule. In AGFRgui they do not, so a choice made after a receptor had
    been written was ignored until some checkbox changed.
  * An edited manual charge is used. AGFRgui's _getCharges reads
    checkState(1) of the charge items, a column that never gets a check state,
    so the spin box value was never passed on and the default charge
    (PDBProfiler._manualCharges) was always written. With unedited charges
    both give the same file (tests/test_parity.py checks this).
  * The "unknown atom types" subtree is display-only: its residues are not
    user-checkable in AGFRgui either, and their unsupported atoms cannot be
    written to a PDBQT.
"""
from .tree import CHECKED, LIG, REC, UNCHECKED


class Refused(Exception):
    """The request names something AGFRgui would not let the user click."""


def invalidate(agfr, receptor=True, ligand=True):
    if receptor:
        agfr._processedReceptor = None
    if ligand:
        agfr._adLigs = []
        agfr._processedLigands = []
    agfr._BMProcessed = None


# ---- setRecFlag / setLigFlag -------------------------------------------
def _chosen(node):
    return node.altloc['alts'][node.altloc['current']]


def _set_leaf(node, flag, other, other_col, on):
    sel = node.sel
    if sel is None:
        return
    if on:
        # none of the atoms of this residue can be the other thing
        sel.setFlags(other, False)
        if node.checkable:
            node.set_state(other_col, UNCHECKED)
        if node.altloc is not None:
            # only the atoms with no altloc or the chosen one
            sel.setFlags(flag, False)
            atoms = sel.select('altloc _ %s' % _chosen(node))
            if atoms is not None:
                atoms.setFlags(flag, True)
        else:
            sel.setFlags(flag, True)
    else:
        sel.setFlags(flag, False)


def set_rec_flag(node, on):
    if node.children:
        for c in node.children:
            set_rec_flag(c, on)
    else:
        _set_leaf(node, 'receptor', 'ligand', LIG, on)


def set_lig_flag(node, on):
    if node.children:
        for c in node.children:
            set_lig_flag(c, on)
    else:
        _set_leaf(node, 'ligand', 'receptor', REC, on)


def set_item_flags(agfr, node, flag_name, on):
    """_setItemFlags."""
    assert flag_name in ('receptor', 'ligand')
    if flag_name == 'receptor':
        set_rec_flag(node, on)
        invalidate(agfr, receptor=True, ligand=False)
    else:
        set_lig_flag(node, on)
        invalidate(agfr, receptor=False, ligand=True)


# ---- the user actions ----------------------------------------------------
def click(agfr, node, col):
    """The checkbox of `node` in column `col` (REC or LIG) was clicked."""
    if not (node.live and node.checkable) or node.sel is None:
        raise Refused('%s cannot be checked' % node.label)
    cur = node.check_state(col)
    if cur is None:
        raise Refused('%s has no %s checkbox' % (
            node.label, 'receptor' if col == REC else 'ligand'))
    # Qt: a non-user-tristate box goes checked -> unchecked, else -> checked
    new = UNCHECKED if cur == CHECKED else CHECKED
    node.set_state(col, new)         # Qt pushes it down the subtree first
    set_item_flags(agfr, node, 'receptor' if col == REC else 'ligand',
                   new == CHECKED)


def choose_altloc(agfr, profile, node, index):
    """AltlocsChooser.onclick: another alternate location was chosen."""
    alts = node.altloc['alts']
    index = int(index)
    if not 0 <= index < len(alts):
        raise Refused('no alternate location %d for %s' % (index, node.label))
    node.altloc['current'] = index
    set_item_flags(agfr, node, 'receptor', node.check_state(REC) == CHECKED)
    set_item_flags(agfr, node, 'ligand', node.check_state(LIG) == CHECKED)
    profile._altlocsRes[node.key][0] = index


def sync_altlocs(agfr, profile, tree):
    """switchBM: altloc choices are shared by all tabs through
    PDBProfiler._altlocsRes; a choice made in another tab is re-applied here
    (setChecked on the radio button fires onclick when it changes)."""
    for key, node in tree.altloc_nodes():
        ind = profile._altlocsRes[key][0]
        if node.altloc['current'] != ind:
            choose_altloc(agfr, profile, node, ind)


def apply_altlocs(agfr, profile, tree, initial):
    """A tab built for the first time: its atoms carry the altloc flagged
    when the file was read (`initial`), its tree shows _altlocsRes. Re-apply
    the residues where those differ."""
    for key, node in tree.altloc_nodes():
        if node.altloc['current'] != initial.get(key, node.altloc['current']):
            choose_altloc(agfr, profile, node, node.altloc['current'])


def set_sc_completion(agfr, profile, node, act):
    """_setSCcompletion: 'complete' or 'asis' for one residue or a category."""
    if act not in ('complete', 'asis'):
        raise Refused('unknown side chain action %r' % act)
    if node.menu != 'sc':
        raise Refused('%s has no side chain menu' % node.label)
    items = node.children if node.kind == 'category' else [node]
    for item in items:
        res = item.sel
        key = item.key
        atoms = ' '.join(profile._missingAtoms.get(key, []))
        try:
            # ADFR's definition, the one agfr --sidechainsAsIs uses
            profile.completeSidechain(key, act == 'complete')
        except ValueError as e:
            raise Refused(str(e))
        if act == 'complete':
            item.label = "%s'%s%s build:%s" % (res.getResname(), res.getResnum(),
                                              res.getIcode(), atoms)
        else:
            item.label = "%s'%s%s use as is missing:%s" % (
                res.getResname(), res.getResnum(), res.getIcode(), atoms)
    invalidate(agfr)


def mutate(agfr, profile, node, choice):
    """_selected: an entry of the mutate menu was picked."""
    if node.menu != 'mutate':
        raise Refused('%s has no mutate menu' % node.label)
    res = node.sel
    resnum = res.getResnum()
    resname = res.getResname()
    parts = choice.split('Base Res:')
    mutResName = (parts[1] if len(parts) > 1 else parts[0]).strip()
    if mutResName == resname:
        node.label = "%s'%s%s" % (resname, resnum, res.getIcode())
        profile._resToMutate.pop(node.key, None)
    else:
        node.label = "%s'%s -> %s" % (resname, resnum, mutResName)
        profile._resToMutate[node.key] = [resname, mutResName]
    invalidate(agfr)


def set_charge(agfr, node, value):
    """The manual-charge spin box (range -10..10, step 0.05, 2 decimals)."""
    if node.charge is None:
        raise Refused('%s takes no manual charge' % node.label)
    value = round(max(-10.0, min(10.0, float(value))), 2)
    node.charge = value
    invalidate(agfr, receptor=True, ligand=False)
    return value


def charges(tree):
    """{residue key: charge} for manual-charge residues in the receptor."""
    return dict((n.key, n.charge) for n in tree.manual_charge_nodes()
                if n.check_state(REC) == CHECKED)


def focus_indices(node):
    """The atoms onClick hands to visuallyIdentify."""
    obj = None
    if node.sel is not None:
        obj = node.sel
        if node.altloc is not None:
            obj = obj.select('receptor')
    elif node.gapsel is not None:
        obj = node.gapsel
    return [] if obj is None else [int(i) for i in obj.getIndices()]
