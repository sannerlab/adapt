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

import numpy

from . import env  # noqa: F401  (sets ADS_ROOT before ADFR is used)
from . import copies, flags, gaps, prep, reduce_log
from .tree import CHECKED, LIG, REC, Tree

# prody, OpenBabel, reduce's pipes and redirect_stdout are all process-wide.
# The page is one user's, so one lock costs nothing and removes every race.
LOCK = threading.RLock()

PDB_SUFFIXES = ('.pdb', '.ent', '.pdb.gz', '.ent.gz', '.cif', '.cif.gz', '.mmcif')


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


class RestoreError(Exception):
    """A session folder cannot be restored exactly."""


def safe_name(name):
    base = os.path.basename(name or 'input.pdb')
    base = re.sub(r'[^A-Za-z0-9._-]', '_', base).lstrip('.') or 'input.pdb'
    if not base.lower().endswith(PDB_SUFFIXES):
        raise ValueError('%s: ADapt reads PDB files (.pdb, .ent) and mmCIF files '
                         '(.cif, .mmcif), gzipped or not'
                         % base)
    return base


class Session(object):

    def __init__(self, sid, data, filename, outroot, provenance=None, source=None,
                 folder=None):
        """A new session of `data` (the file `filename`) in a new folder under
        `outroot`; or, with `folder`, the session of that folder (restore():
        nothing in it is written until its edits are replayed)."""
        self.id = sid
        self.nodes = {}
        self._next = 0
        self.lock = LOCK
        self.log = []
        self._logsent = 0
        self._replaying = folder is not None
        name = safe_name(filename)
        if folder is None:
            stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
            stem = name.split('.')[0]
            folder = os.path.join(outroot, '%s_%s_%s' % (stem, stamp, sid[:6]))
            os.makedirs(folder)
            # the name is kept: PDBProfiler takes the PDB id from it, and the id
            # selects the known ligands of that entry (ADFR.PDBLigDicts.pdb2lig)
            with open(os.path.join(folder, name), 'wb') as fh:
                fh.write(data)
        self.folder = folder
        self.input = os.path.join(self.folder, name)
        self.saves = []
        self.events = []
        from . import APP, __version__
        self.record = {
            'schema': 2,
            'session': sid,
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

    # ---- restoring a session from its folder ----------------------------
    @classmethod
    def restore(cls, folder, sid=None, provenance=None):
        """The session recorded in `folder`, its edits replayed: after a
        server restart, or for a session another server made. Refused (a
        RestoreError naming the event) when the input file changed or an edit
        cannot be replayed exactly; nothing in the folder is changed then."""
        path = os.path.join(folder, 'session.json')
        with open(path) as fh:
            rec = json.load(fh)
        name = rec['input']['name']
        with open(os.path.join(folder, name), 'rb') as fh:
            data = fh.read()
        if hashlib.md5(data).hexdigest() != rec['input']['md5']:
            raise RestoreError('%s changed since the session was made' % name)
        sid = sid or rec.get('session')
        if not sid:
            raise RestoreError('the session in %s has no id (made before restore '
                               'existed); restore it with a new one' % folder)
        s = cls(sid, data, name, None, provenance, folder=folder)
        s._replay(rec)
        return s

    def _replay(self, rec):
        saves = dict((e.get('after_event'), e) for e in rec['saves'] if 'after_event' in e)
        try:
            for k, e in enumerate(rec['events']):
                if k in saves:
                    self._check_state(saves[k])
                self._apply(k, e)
            if len(rec['events']) in saves:
                self._check_state(saves[len(rec['events'])])
        finally:
            self._replaying = False
        # from now on the session continues its own record
        self.record = rec
        self.events = rec['events']
        self.saves = rec['saves']
        self.record['session'] = self.id
        self._event('restored', events=len(self.events))

    def _check_state(self, save):
        if save.get('state_md5') and save['state_md5'] != self._state_md5():
            raise RestoreError('the replayed edits do not give the state of save %d'
                               % save['n'])

    def _apply(self, k, e):
        what = e.get('what')
        if what in ('restored', None):
            return
        if e.get('tab') and e['tab'] != self.current and what != 'tab':
            self.switch(e['tab'])
        if what == 'tab':
            self.switch(e['name'])
            return
        node = self._find(k, e)
        if what == 'check':
            col = {'rec': REC, 'lig': LIG}[e['column']]
            for _ in range(3):
                if node.check_state(col) == e['state']:
                    break
                self.click(node.id, e['column'])
            if node.check_state(col) != e['state']:
                raise RestoreError('event %d (%s %s) does not give the recorded state'
                                   % (k + 1, what, e.get('node')))
        elif what == 'altloc':
            self.altloc(node.id, node.altloc['alts'].index(e['altloc']))
        elif what == 'charge':
            self.charge(node.id, e['charge'])
        elif what == 'sidechain':
            self.sidechain(node.id, e['action'])
        elif what == 'mutate':
            self.mutate(node.id, e['choice'])
        elif what == 'rotamer':
            self.rotamer(node.id, -1 if e['rotamer'] == 'auto' else int(e['rotamer']) - 1)
        elif what == 'substitute':
            self.substitute(node.id, e.get('element'), e.get('charge'))
        else:
            raise RestoreError('event %d: unknown edit %r' % (k + 1, what))

    def _find(self, k, e):
        """The node an event was made on: by its ref; for sessions recorded
        before refs, by its residue key, else by its label if no other node of
        the tab has it."""
        tree = self.tree()
        if e.get('ref'):
            node = tree.resolve(e['ref'])
        elif e.get('key'):
            node = tree.resolve({'key': e['key']})
        else:
            same = [n for n in tree.walk() if n.label == e.get('node')]
            if len(same) != 1:
                raise RestoreError('event %d (%s %s) names %d items of tab %s; it '
                                   'cannot be replayed' % (k + 1, e.get('what'), e.get('node'),
                                                           len(same), self.current))
            node = same[0]
        if node is None:
            raise RestoreError('event %d (%s %s): no such item in tab %s'
                               % (k + 1, e.get('what'), e.get('node'), self.current))
        return node

    def _state_md5(self):
        """A digest of what a save of the current tab writes from: the atoms'
        receptor and ligand flags and charges, and the choices."""
        ag = self.tree()._ag
        h = hashlib.md5()
        h.update(self.current.encode())
        h.update(ag.getFlags('receptor').tobytes())
        h.update(ag.getFlags('ligand').tobytes())
        q = ag.getCharges()
        if q is not None:
            h.update(numpy.round(q, 4).tobytes())
        p = self.profile
        h.update(json.dumps([sorted(p._resToMutate.items()),
                             sorted((k, v[1][v[0]]) for k, v in p._altlocsRes.items()),
                             sorted(p._rotamerChoice.items()),
                             sorted(p._substitutions.items())], default=str).encode())
        return h.hexdigest()

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

    def _mirror(self):
        """The symmetry copies of the current tab follow its original chains
        (adapt/copies.py)."""
        tree = self.tree()
        copies.mirror(tree._ag, tree._copies)

    def _event(self, what, **kw):
        if self._replaying:
            return
        node = kw.pop('ref_node', None)
        if node is not None:
            kw['ref'] = node.ref()
        kw['what'] = what
        kw['tab'] = self.current
        kw['at'] = datetime.datetime.now().isoformat(timespec='seconds')
        self.events.append(kw)
        self.flush()

    def flush(self):
        if self._replaying:
            return
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
            copies.mirror(tree._ag, tree._copies)
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
        r receptor, l ligand, o neither; upper case for the original chains
        (the others are the biomolecule's symmetry copies)."""
        ag = self.tree(name)._mol._ag
        rec = ag.getFlags('receptor')
        lig = ag.getFlags('ligand')
        orig = ag.getFlags(copies.ORIGINAL)
        out = []
        for r, l, o in zip(rec, lig, orig):
            c = 'r' if r else ('l' if l else 'o')
            out.append(c.upper() if o else c)
        return ''.join(out)

    def view(self, full=False, name=None):
        tree = self.tree(name)
        d = {'tab': tree.name, 'states': tree.states(), 'buttons': self.buttons(),
             'atoms': self.atom_states(name)}
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
            self._mirror()
        self._event('check', ref_node=node, node=node.label, key=node.key, column=column,
                    state=node.check_state(col))
        return self.view()

    def altloc(self, nid, index):
        node = self.node(nid)
        if node.altloc is None:
            raise flags.Refused('%s has no alternate locations' % node.label)
        with captured(self.log):
            flags.choose_altloc(self.agfr, self.profile, node, index)
            self._mirror()
        self._event('altloc', ref_node=node, node=node.label, key=node.key,
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
            self._mirror()
        self._event('charge', ref_node=node, node=node.label, key=node.key, charge=v)
        d = self.view()
        d['value'] = v
        return d

    def sidechain(self, nid, act):
        node = self.node(nid)
        with captured(self.log):
            flags.set_sc_completion(self.agfr, self.profile, node, act)
            self._prune_rotamers()
        self._event('sidechain', ref_node=node, node=node.label, key=node.key, action=act)
        return self.view(full=True)

    def mutate(self, nid, choice):
        node = self.node(nid)
        with captured(self.log):
            flags.mutate(self.agfr, self.profile, node, choice)
            self._prune_rotamers()
        self._event('mutate', ref_node=node, node=node.label, key=node.key, choice=choice)
        return self.view(full=True)

    def mutate_options(self, nid):
        from .tree import mutate_options
        node = self.node(nid)
        if node.menu != 'mutate':
            raise flags.Refused('%s has no mutate menu' % node.label)
        return mutate_options(self.tree(), node)

    def rotamer(self, nid, choice):
        """Show rotamer `choice` (-1: auto, the best-scoring one) of a side
        chain to complete, and build that one when the receptor is written."""
        node = self.node(nid)
        if node.rotamer is None:
            raise flags.Refused('%s has no rotamers' % node.label)
        choice = int(choice)
        with LOCK:
            try:
                self.profile.chooseRotamer(node.key, None if choice < 0 else choice)
            except ValueError as e:
                raise flags.Refused(str(e))
            pdb, index = node.rotamer.pdb(choice)
            node.rotamer.current = choice
            flags.invalidate(self.agfr)
        rot = 'auto' if choice < 0 else choice + 1
        last = self.events[-1] if self.events else {}
        if last.get('what') == 'rotamer' and last.get('key') == node.key and not self._replaying:
            last['rotamer'] = rot          # one event per residue, not per wheel step
            self.flush()
        else:
            self._event('rotamer', ref_node=node, node=node.label, key=node.key, rotamer=rot)
        return {'pdb': pdb, 'index': index, 'message': node.rotamer.message(index),
                'auto': node.rotamer.best}

    def substitute(self, nid, element, charge=None):
        """Write the one-atom residue of node `nid`, an ion AD4.1 cannot type
        (copper, cobalt, ...), as an ion of `element` (Zn, Mg, Ca, Fe, Mn;
        None undoes it). The tree is rebuilt: the ion moves to manual charges."""
        node = self.node(nid)
        if node.menu != 'substitute':
            raise flags.Refused('%s cannot be replaced' % node.label)
        q = None if charge in (None, '') else float(charge)
        with captured(self.log):
            try:
                self.profile.substitute(node.key, element or None, q)
            except ValueError as e:
                raise flags.Refused(str(e))
            flags.invalidate(self.agfr)
            self._rebuild()
        sub = self.profile._substitutions.get(node.key)
        self._event('substitute', ref_node=node, node=node.label, key=node.key,
                    element=sub[0] if sub else None, charge=sub[1] if sub else None)
        return self.view(full=True)

    def _rebuild(self):
        """New trees for every tab (the current one built now): after a
        substitution a residue changes category."""
        self.nodes.clear()
        current, self.current = self.current, None
        self.trees = OrderedDict((n, Tree(self, n)) for n in self.profile._assemblies)
        self.switch(current)

    def _prune_rotamers(self):
        # a chosen rotamer belongs to a residue completed or mutated as that type
        for key, (restype, index) in list(self.profile._rotamerChoice.items()):
            if self.profile._resToMutate.get(key, (None, None))[1] != restype:
                self.profile.chooseRotamer(key, None)

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
        reduce = None
        out = getattr(self.agfr, 'reduceOutput', None)
        if kind == 'receptor' and out:
            reduce = {'version': reduce_log.version(out),
                      'chain_renames': reduce_log.chain_renames(self.agfr.reduceChainRenames),
                      'decisions': reduce_log.decisions(out)}
        entry = {'n': n, 'kind': kind, 'tab': self.current,
                 'at': datetime.datetime.now().isoformat(timespec='seconds'),
                 'folder': os.path.basename(folder), 'files': files,
                 'mutations': dict(self.profile._resToMutate),
                 'altlocs': dict((k, v[1][v[0]]) for k, v in
                                 self.profile._altlocsRes.items()),
                 'charges': flags.charges(self.tree()),
                 'rotamers': dict((k, i + 1) for k, (t, i) in
                                  self.profile._rotamerChoice.items()),
                 'substitutions': dict((k, list(v)) for k, v in
                                       self.profile._substitutions.items()),
                 'log': log}
        entry['after_event'] = len(self.events)       # for restore: replay up to here
        entry['state_md5'] = self._state_md5()        # ... and check the state
        if kind == 'receptor':
            # residues dropped for an atom type AD4.1 does not have
            entry['left_out'] = dict(self.profile._leftOut)
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
