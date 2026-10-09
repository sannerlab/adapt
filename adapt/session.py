# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""One uploaded PDB file and everything the panel does to it.

A Session is what GridGUI.loadReceptor + PDBAnalyserWidget hold in AGFRgui:
a runAGFR instance configured with agfrgui's own option defaults, its
PDBProfiler, and one residue tree per biomolecule. Its folder (on NFS, never
/tmp) keeps the input file, every PDBQT written, and session.json, which
records the ADFR code that ran and every edit in order, so any output can be
traced back and reproduced.
"""
import contextlib
import datetime
import hashlib
import io
import json
import os
import re
import threading
import traceback
from collections import OrderedDict

from . import env  # noqa: F401  (sets ADS_ROOT before ADFR is used)
from . import flags, gaps, prep, rebuild, reduce_log
from .tree import CHECKED, LIG, REC, Tree

# prody, OpenBabel, reduce's pipes and redirect_stdout are all process-wide.
# The page is one user's, so one lock costs nothing and removes every race.
LOCK = threading.RLock()

# ADFR's gap detection never fires (adapt/gaps.py); use the fixed copy
gaps.install()
# keep what reduce writes when a receptor is protonated (adapt/reduce_log.py)
reduce_log.install()

PDB_SUFFIXES = ('.pdb', '.ent', '.pdb.gz', '.ent.gz')


def gui_options():
    """The options AGFRgui starts with: agfrgui's parser, no arguments
    (ADFR/GUI/__init__.py:makeGridGUI -> GridGUI.configure)."""
    from ADFR.utils.optParser import AGFRParser
    opts = vars(AGFRParser(prog='agfrgui').parse_args([]))
    opts.pop('debug', None)
    opts.pop('inputFolder', None)
    return opts


@contextlib.contextmanager
def captured(log):
    """Run ADFR with its prints collected into `log` (a list of lines)."""
    buf = io.StringIO()
    with LOCK, contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            yield
        finally:
            log.extend(l for l in buf.getvalue().splitlines()
                       if l.strip() and 'BiopythonDeprecationWarning' not in l
                       and 'pairwise2' not in l)


def safe_name(name):
    base = os.path.basename(name or 'input.pdb')
    base = re.sub(r'[^A-Za-z0-9._-]', '_', base).lstrip('.') or 'input.pdb'
    if not base.lower().endswith(PDB_SUFFIXES):
        raise ValueError('%s: ADapt reads .pdb, .ent, .pdb.gz or .ent.gz files'
                         % base)
    return base


class Session(object):

    def __init__(self, sid, data, filename, outroot, provenance=None, source=None):
        self.id = sid
        self.nodes = {}
        self._next = 0
        self.lock = LOCK
        self.log = []
        self._logsent = 0
        name = safe_name(filename)
        stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        stem = name.split('.')[0]
        self.folder = os.path.join(outroot, '%s_%s_%s' % (stem, stamp, sid[:6]))
        os.makedirs(self.folder)
        # the name is kept: PDBProfiler takes the PDB id from it, and the id
        # selects the known ligands of that entry (ADFR.PDBLigDicts.pdb2lig)
        self.input = os.path.join(self.folder, name)
        with open(self.input, 'wb') as fh:
            fh.write(data)
        self.saves = []
        self.events = []
        from . import APP, __version__
        self.record = {
            'tool': {'name': APP, 'version': __version__},
            'created': datetime.datetime.now().isoformat(timespec='seconds'),
            'input': {'name': name, 'bytes': len(data),
                      'md5': hashlib.md5(data).hexdigest(),
                      'source': source or {'upload': filename}},
            'environment': provenance,
            'events': self.events,
            'saves': self.saves,
        }

        from ADFR.utils.runAGFR import runAGFR
        with captured(self.log):
            self.options = gui_options()
            self.agfr = runAGFR()
            self.agfr.proccessPreLoadCmdLineOptions(self.options)
            self.agfr.loadReceptor(self.input, processReceptor=False)
        if not self.agfr._pdbInput:
            raise ValueError('%s was not read as a PDB file' % name)
        self.profile = self.agfr.PDBprofile
        # the internal gaps, and the ones the user asks to rebuild (rebuild.py)
        self.rebuilder = rebuild.Rebuilder(self.profile)
        # the altloc each residue was flagged with when the file was read
        self._altloc0 = dict((k, v[0]) for k, v in self.profile._altlocsRes.items())
        self.record['options'] = _jsonable(self.options)
        self.record['pdbid'] = self.profile._pdbid
        self.trees = OrderedDict((n, Tree(self, n))
                                 for n in self.profile._assemblies)
        self.current = None
        first = 'Biomol 1' if 'Biomol 1' in self.trees else list(self.trees)[0]
        self.switch(first)
        self.record['opened_tab'] = first
        self.flush()        # the record exists from the moment a file is opened

    def _newid(self):
        self._next += 1
        return self._next

    def tree(self, name=None):
        return self.trees[name or self.current]

    def node(self, nid):
        try:
            n = self.nodes[int(nid)]
        except (KeyError, TypeError, ValueError):
            raise flags.Refused('no such tree item %r' % (nid,))
        if n.label is None:
            raise flags.Refused('no such tree item %r' % (nid,))
        return n

    def _event(self, what, **kw):
        kw['what'] = what
        kw['tab'] = self.current
        kw['at'] = datetime.datetime.now().isoformat(timespec='seconds')
        self.events.append(kw)
        self.flush()

    def flush(self):
        path = os.path.join(self.folder, 'session.json')
        tmp = path + '.part'
        with open(tmp, 'w') as fh:
            json.dump(self.record, fh, indent=1)
        os.replace(tmp, path)

    # ---- header (PDBAnalyserWidget.buildUI label) ----------------------
    def label(self):
        hdr = self.profile._receptor.pdbHeader or {}
        exp = hdr.get('experiment', None)
        txt = '%s.pdb ' % self.profile._pdbid
        if exp:
            if exp.find('NMR') >= 0:
                txt += '%s, n models: %d' % (exp, hdr.get('n_models', 1))
            elif exp.find('CRYSTAL') >= 0 or exp.find('X-RAY') >= 0:
                res = hdr.get('resolution', None)
                txt += ('%s, resolution: %0.2f' % (exp, res)) if res is not None else exp
            else:
                txt += exp
        return txt

    # ---- tabs (switchBM) ------------------------------------------------
    def switch(self, name):
        if name not in self.trees:
            raise flags.Refused('no biomolecule %r' % name)
        tree = self.trees[name]
        with captured(self.log):
            if not tree._hasUI:
                tree.build()
                # AGFRgui builds a new tab's radio buttons from _altlocsRes but
                # leaves that molecule's atoms flagged with the file's default
                # altloc, so a choice made in another tab showed here, went
                # into the REMARK, and was not in the coordinates. Apply it.
                flags.apply_altlocs(self.agfr, self.profile, tree, self._altloc0)
            else:
                flags.sync_altlocs(self.agfr, self.profile, tree)
            self.profile._currentBM = tree._mol
            # AGFRgui keeps _BMProcessed across a tab switch, so a save right
            # after switching wrote the PREVIOUS biomolecule. Drop it here.
            flags.invalidate(self.agfr)
        changed = self.current is not None and self.current != name
        self.current = name
        if changed:
            self._event('tab', name=name)
        return tree

    # ---- what the page shows -------------------------------------------
    def buttons(self):
        """PDBAnalyserWidget._configureButtons."""
        bm = self.profile._currentBM
        return {'receptor': bool(bm.select('receptor')),
                'ligand': bool(bm.select('ligand')),
                'other': bool(bm.select('not receptor and not ligand'))}

    def atom_states(self, name=None):
        """One character per atom of the tab's molecule, in atom-index order:
        r receptor, l ligand, o neither; upper case for segment A (the
        others are the biomolecule's symmetry copies)."""
        ag = self.tree(name)._mol._ag
        rec = ag.getFlags('receptor')
        lig = ag.getFlags('ligand')
        seg = ag.getSegnames()
        out = []
        for r, l, s in zip(rec, lig, seg):
            c = 'r' if r else ('l' if l else 'o')
            out.append(c.upper() if s == 'A' else c)
        return ''.join(out)

    def view(self, full=False, name=None):
        tree = self.tree(name)
        d = {'tab': tree.name, 'states': tree.states(), 'buttons': self.buttons(),
             'atoms': self.atom_states(name),
             # the page redraws the rebuilt residues after an edit while any are chosen
             'rebuild': len(self.rebuilder.chosen_in(tree._chaindInBM))}
        if full:
            d['tree'] = tree.to_json()
        d['log'] = self.new_log()
        return d

    def new_log(self):
        """The log lines not yet sent to the page."""
        lines = self.log[self._logsent:]
        self._logsent = len(self.log)
        return lines

    def info(self):
        return {'session': self.id, 'pdbid': self.profile._pdbid,
                'label': self.label(), 'folder': self.folder,
                'folder_display': _home_tilde(self.folder),
                'tabs': list(self.trees), 'current': self.current}

    def structure(self, name=None):
        import prody
        tree = self.tree(name)
        out = io.StringIO()
        prody.writePDBStream(out, tree._mol._ag)
        with LOCK:
            geo = gaps.geometry(self.profile._gapsPerChain, tree._mol)
        return {'tab': tree.name, 'pdb': out.getvalue(),
                'n': len(tree._mol._ag), 'atoms': self.atom_states(name),
                'gaps': geo}

    # ---- edits ------------------------------------------------------------
    def click(self, nid, column):
        node = self.node(nid)
        col = {'rec': REC, 'lig': LIG}.get(column)
        if col is None:
            raise flags.Refused('column must be rec or lig')
        with captured(self.log):
            flags.click(self.agfr, node, col)
        self._event('check', node=node.label, key=node.key, column=column,
                    state=node.check_state(col))
        return self.view()

    def altloc(self, nid, index):
        node = self.node(nid)
        if node.altloc is None:
            raise flags.Refused('%s has no alternate locations' % node.label)
        with captured(self.log):
            flags.choose_altloc(self.agfr, self.profile, node, index)
        self._event('altloc', node=node.label, key=node.key,
                    altloc=node.altloc['alts'][node.altloc['current']])
        return self.view()

    def charge(self, nid, value):
        node = self.node(nid)
        with captured(self.log):
            v = flags.set_charge(self.agfr, node, value)
            # the atoms carry it too: processBM copies them, and
            # getPDBQTHeaderLines reports the charge of the ASU's atoms
            node.sel.setCharges([v] * len(node.sel))
            chid, resname, resnum, icode = node.key.split(':')
            asu = self.profile._receptor._ag.select(
                'chid %s resnum `%s` icode %s' % (chid, resnum, icode))
            if asu is not None:
                asu.setCharges([v] * len(asu))
        self._event('charge', node=node.label, key=node.key, charge=v)
        d = self.view()
        d['value'] = v
        return d

    def sidechain(self, nid, act):
        node = self.node(nid)
        with captured(self.log):
            flags.set_sc_completion(self.agfr, self.profile, node, act)
        self._event('sidechain', node=node.label, key=node.key, action=act)
        return self.view(full=True)

    def mutate(self, nid, choice):
        node = self.node(nid)
        with captured(self.log):
            flags.mutate(self.agfr, self.profile, node, choice)
        self._event('mutate', node=node.label, key=node.key, choice=choice)
        return self.view(full=True)

    def mutate_options(self, nid):
        from .tree import mutate_options
        node = self.node(nid)
        if node.menu != 'mutate':
            raise flags.Refused('%s has no mutate menu' % node.label)
        return mutate_options(self.tree(), node)

    def rotamer(self, nid, choice):
        node = self.node(nid)
        if node.rotamer is None:
            raise flags.Refused('%s has no rotamer preview' % node.label)
        with LOCK:
            pdb, index = node.rotamer.pdb(int(choice))
            node.rotamer.current = int(choice)
        return {'pdb': pdb, 'index': index, 'message': node.rotamer.message(index),
                'auto': node.rotamer.best}

    def rebuild(self, nid, on):
        """The rebuild box of an internal gap: choose it (or not), and build it
        at once, so the page can show the residues before they are written. A
        loop that is not sane is refused: the box goes back, with the reason."""
        node = self.node(nid)
        if node.rebuild is None:
            raise flags.Refused('%s is not a gap that can be rebuilt' % node.label)
        key, on = node.rebuild['key'], bool(on)
        try:
            self.rebuilder.choose(key, on)
        except ValueError as e:
            raise flags.Refused('%s cannot be rebuilt: %s' % (node.label, e))
        with captured(self.log):
            flags.invalidate(self.agfr)
        try:
            built, refused = self._preview(key)
        except Exception as e:
            # leave nothing chosen that cannot be built
            (self.rebuilder.chosen.discard if on else self.rebuilder.chosen.add)(key)
            with captured(self.log):
                flags.invalidate(self.agfr)
            self.log.append(traceback.format_exc())
            raise flags.Refused('%s could not be rebuilt: %s: %s'
                                % (node.label, type(e).__name__, e))
        if key not in refused:
            self._boxes(key)
            self._event('rebuild', node=node.label, key=key, rebuild=on)
        d = self.view(full=True)
        d['rebuilt'] = built
        d['refused'] = list(refused.values())
        return d

    def rebuilt(self, name=None):
        """The rebuilt residues of the current tab, as the receptor would get
        them now (rebuild.preview); the same molecule a save writes. An edit
        can leave a chosen gap's loop not sane: it is refused on the way, and
        the page gets the tree again."""
        if name and name != self.current:
            raise flags.Refused('rebuilt residues are drawn for the current tab only')
        built, refused = self._preview()
        d = {'rebuilt': built}
        if refused:
            d.update(self.view(full=True))
            d['refused'] = list(refused.values())
        return d

    def _preview(self, asked=None):
        """(rebuild.preview of the processed molecule, {key: message} of the
        chosen gaps refused while making it)."""
        refused = OrderedDict()
        while True:
            try:
                with captured(self.log):
                    mol = prep.processed_bm(self)
                return rebuild.preview(mol), refused
            except rebuild.NotBuilt as e:
                refused[e.key] = self._refuse(e, asked)

    def _refuse(self, e, asked=None):
        """A chosen gap whose loop was not sane (rebuild.NotBuilt e): no longer
        chosen, its box says why, and the record says so."""
        self.rebuilder.refuse(e.key, str(e))
        with captured(self.log):
            flags.invalidate(self.agfr)
        self._boxes(e.key)
        label = next((n.label for n in self.nodes.values()
                      if n.rebuild is not None and n.rebuild['key'] == e.key), 'gap ' + e.key)
        self.log.append('NOT REBUILT %s: %s' % (label, e))
        self._event('rebuild', node=label, key=e.key, rebuild=False, refused=str(e),
                    asked=e.key == asked)
        return '%s was not rebuilt: %s' % (label, self.rebuilder.gaps[e.key].why)

    def _boxes(self, key):
        """The rebuild box of a gap, in every tab built so far."""
        for n in self.nodes.values():
            if n.rebuild is not None and n.rebuild['key'] == key:
                n.rebuild['on'] = key in self.rebuilder.chosen
                n.rebuild['why'] = self.rebuilder.gaps[key].why

    def focus(self, nid):
        return {'indices': flags.focus_indices(self.node(nid))}

    def select(self, what):
        """The atoms showCurrentReceptor / showCurrentLigand / showOther show."""
        sel = {'receptor': 'receptor', 'ligand': 'ligand',
               'other': 'not ligand and not receptor'}.get(what)
        if sel is None:
            raise flags.Refused('show receptor, ligand or other')
        atoms = self.profile._currentBM.select(sel)
        return {'indices': [] if atoms is None else [int(i) for i in atoms.getIndices()]}

    # ---- the save buttons -------------------------------------------------
    def save(self, kind):
        if kind not in ('receptor', 'ligand'):
            raise flags.Refused('save receptor or ligand')
        n = len(self.saves) + 1
        folder = os.path.join(self.folder, 'save_%03d_%s' % (n, kind[:3]))
        os.makedirs(folder)
        log = []
        try:
            with captured(log):
                paths = (prep.receptor if kind == 'receptor' else prep.ligands)(self, folder)
        except prep.Nothing as e:
            os.rmdir(folder)     # empty: nothing was written into it
            raise flags.Refused(str(e))
        except rebuild.NotBuilt as e:
            os.rmdir(folder)     # the molecule was made before any file
            self.log.extend(log)
            raise flags.Refused('%s. Nothing was written; save again to write without it.'
                                % self._refuse(e))
        except Exception:
            log.append(traceback.format_exc())
            self.log.extend(log)
            self.saves.append({'n': n, 'kind': kind, 'tab': self.current,
                               'error': log[-1], 'log': log})
            self.flush()
            raise
        self.log.extend(log)
        files = []
        for p in paths:
            with open(p, 'rb') as fh:
                files.append({'name': os.path.basename(p),
                              'md5': hashlib.md5(fh.read()).hexdigest()})
        rebuilt = None
        if kind == 'receptor':
            # what was rebuilt and written, as the REMARK lines say it
            rebuilt = [dict((k, v) for k, v in r.items() if k not in ('atoms', 'written'))
                       for r in getattr(self.agfr._BMProcessed, 'rebuilt', ()) if r['written']]
        reduce = None
        if kind == 'receptor' and reduce_log.LAST.get('stdout'):
            out = reduce_log.LAST['stdout']
            reduce = {'version': reduce_log.version(out),
                      'chain_renames': reduce_log.chain_renames(reduce_log.LAST.get('taken')),
                      'decisions': reduce_log.decisions(out)}
        entry = {'n': n, 'kind': kind, 'tab': self.current,
                 'at': datetime.datetime.now().isoformat(timespec='seconds'),
                 'folder': os.path.basename(folder), 'files': files,
                 'mutations': dict(self.profile._resToMutate),
                 'altlocs': dict((k, v[1][v[0]]) for k, v in
                                 self.profile._altlocsRes.items()),
                 'charges': flags.charges(self.tree()),
                 'log': log}
        if rebuilt:
            entry['rebuilt'] = rebuilt
        if reduce is not None:
            entry['reduce'] = reduce
        self.saves.append(entry)
        self.flush()
        return {'save': n, 'files': files, 'log': log}

    def save_paths(self, n):
        for s in self.saves:
            if s['n'] == int(n) and 'files' in s:
                folder = os.path.join(self.folder, s['folder'])
                return s, [os.path.join(folder, f['name']) for f in s['files']]
        raise flags.Refused('no save %r' % (n,))


def _jsonable(d):
    out = {}
    for k, v in d.items():
        try:
            json.dumps(v)
            out[k] = v
        except TypeError:
            out[k] = repr(v)
    return out


def _home_tilde(path):
    """A path as shown to the user: the home folder written as ~."""
    home = os.path.expanduser('~')
    if home and home != '/' and (path == home or path.startswith(home + os.sep)):
        return '~' + path[len(home):]
    return path
