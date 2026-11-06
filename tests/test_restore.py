# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""A session survives a server restart: Session.restore replays the edits
its folder records and the session continues where it stopped. It used to
live in the server's memory only ("no such session; load the PDB file
again")."""
import json
import os
import shutil

import numpy

from common import CHECKED, DATA, LIG, REC, UNCHECKED, check, find, outdir, quiet, read, set_state
from adapt import env
from adapt.server import App
from adapt.session import RestoreError, Session

_n = [0]


def _new(label):
    _n[0] += 1
    with open(os.path.join(DATA, '6ggd.pdb'), 'rb') as fh:
        data = fh.read()
    with quiet():
        return Session('%016x' % (0xabc000 + _n[0]), data, '6ggd.pdb', outdir(label),
                       env.provenance())


def _edit_everything(s):
    """One edit of every kind, in Biomol 2 (chain B)."""
    s.switch('Biomol 2')
    t = s.tree()
    water = [n for n in t.walk() if n.label == 'water' and n.kind == 'category'][0]
    set_state(s, water, REC, CHECKED)                       # a category: no residue key
    zn = find(t, key='B:ZN:401:_')[0]
    set_state(s, zn, REC, CHECKED)
    s.charge(zn.id, 1.5)
    alt = [n for n in t.walk() if n.altloc is not None and len(n.altloc['alts']) > 1][0]
    s.altloc(alt.id, (alt.altloc['current'] + 1) % len(alt.altloc['alts']))
    sc = [n for n in t.walk() if n.rotamer is not None]
    s.sidechain(sc[0].id, 'asis')
    s.mutate(find(t, key='B:SER:96:_')[0].id, 'ILE')
    s.rotamer(sc[1].id, 2)
    return s


def _states(s):
    return dict((json.dumps(n.ref(), sort_keys=True), (n.check_state(REC), n.check_state(LIG)))
                for n in s.tree().walk())


def _file(s, d, kind='receptor'):
    return read(os.path.join(s.folder, 'save_%03d_%s' % (d['save'], kind[:3]),
                             d['files'][0]['name']))


def test_restart_restores_everything():
    """Edits of every kind, a save, then one more edit and a save: the
    restored session has the same tab, tree states and atom flags, and a new
    receptor save is byte-identical to the one before the restart."""
    a = _edit_everything(_new('restore_all'))
    rec_a = a.save('receptor')
    set_state(a, find(a.tree(), key='B:EYB:402:_')[0], LIG, CHECKED)
    a.save('ligand')
    with quiet():
        b = Session.restore(a.folder, provenance=env.provenance())
    check(b.id == a.id and b.current == a.current == 'Biomol 2', (b.current, a.current))
    check(_states(b) == _states(a), 'tree states differ')
    for flag in ('receptor', 'ligand'):
        check((b.tree()._ag.getFlags(flag) == a.tree()._ag.getFlags(flag)).all(), flag)
    check(numpy.allclose(b.tree()._ag.getCharges(), a.tree()._ag.getCharges()), 'charges')
    rec_b = b.save('receptor')
    check(_file(b, rec_b) == _file(a, rec_a), 'receptor differs after the restore')
    events = json.load(open(os.path.join(a.folder, 'session.json')))['events']
    check(events[-1]['what'] == 'restored' and len(b.saves) == 3, (events[-1:], len(b.saves)))


def test_changed_input_refused():
    a = _new('restore_changed')
    with open(a.input, 'ab') as fh:
        fh.write(b'REMARK changed\n')
    try:
        Session.restore(a.folder)
    except RestoreError as e:
        check('changed' in str(e), str(e))
    else:
        raise AssertionError('a changed input was restored')


def _old_style(a, label):
    """A copy of a's folder as recorded before refs, ids and state digests."""
    folder = os.path.join(outdir(label), os.path.basename(a.folder))
    shutil.copytree(a.folder, folder)
    path = os.path.join(folder, 'session.json')
    rec = json.load(open(path))
    rec.pop('session')
    for e in rec['events']:
        e.pop('ref', None)
    for e in rec['saves']:
        e.pop('after_event', None)
        e.pop('state_md5', None)
    json.dump(rec, open(path, 'w'))
    return folder


def test_old_sessions_replay_what_is_unambiguous():
    """Before refs, a category click recorded only its label: 'water' in
    Biomol 2 (one chain) replays; 'protein' in the ASU (under chains A and B,
    and under their altloc groups) is refused, naming the event."""
    a = _new('restore_old_ok')
    a.switch('Biomol 2')
    set_state(a, [n for n in a.tree().walk() if n.label == 'water'][0], REC, CHECKED)
    with quiet():
        b = Session.restore(_old_style(a, 'restore_old_ok_copy'), sid='0123456789abcdef')
    check(_states(b) == _states(a), 'old-style session not replayed')
    c = _new('restore_old_bad')
    c.switch('ASU')
    set_state(c, [n for n in c.tree().walk() if n.label == 'protein'][0], REC, UNCHECKED)
    try:
        Session.restore(_old_style(c, 'restore_old_bad_copy'), sid='0123456789abcdee')
    except RestoreError as e:
        check('protein' in str(e) and 'cannot be replayed' in str(e), str(e))
    else:
        raise AssertionError('an ambiguous event was replayed')


def test_server_restart():
    """A new server on the same output folder lists the session and opens it
    by its id (api/session?s=, the page's ?session=)."""
    folder = outdir('restore_server')
    app = App(folder)
    with open(os.path.join(DATA, '6ggd.pdb'), 'rb') as fh:
        with quiet():
            sid = app.upload('6ggd.pdb', fh.read())['session']
    _edit_everything(app.session(sid))
    states = _states(app.session(sid))
    app2 = App(folder)                                   # the restart
    listed = [d for d in app2.list_sessions()['sessions'] if d['session'] == sid]
    check(listed and not listed[0]['live'], listed)
    with quiet():
        s = app2.session(sid)
    check(_states(s) == states, 'restored through the server: states differ')
    try:
        App(folder).session('feedfacefeedface')
    except Exception as e:
        check('no such session' in str(e), str(e))
