# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""What reduce wrote when ADFR protonated a receptor, read for the session.

ADFR's protonate runs `reduce -build -DB ... -` and keeps its output, a PDB
text with reduce's USER MOD records (its choice for every His, Asn and Gln,
among others) and every atom with its hydrogens, and the chain renames it made
for reduce (copies of a chain in an assembly get new chain ids): after
processReceptor they are runAGFR.reduceOutput and runAGFR.reduceChainRenames.
This module reads them: version(), decisions() and chain_renames() for the
session record, renumbered() for the receptor PDB file.
"""
import re


# ---- reading reduce's output -------------------------------------------------
_DECISION = re.compile(
    r'^USER  MOD (?P<set>Set\s*\d+\.\d+|Single)\s*:\s?(?P<chain>.)(?P<resnum>[ \d-]{4})'
    r'(?P<icode>.)(?P<resname>HIS|ASN|GLN)\s[^:]*:(?P<decision>[^:]*):sc=\s*(?P<score>\S+)'
    r'(?:\s+(?P<category>[KCXF])\((?P<detail>[^)]*)\))?')


def renumbered(stdout):
    """reduce's output with its atoms numbered 1, 2, ... in file order.
    reduce writes the hydrogens it adds with serial number 0; nothing in the
    text refers to a serial (no CONECT), so only columns 7-11 of the ATOM and
    HETATM lines change. From the 100000th atom on, the serial is written in
    hexadecimal, as ProDy writes it (and reads it back)."""
    out, n = [], 0
    for l in stdout.splitlines(True):
        if l.startswith(('ATOM  ', 'HETATM')):
            n += 1
            l = l[:6] + ('%5d' % n if n < 100000 else '%5x' % n) + l[11:]
        out.append(l)
    return ''.join(out)


def version(stdout):
    """'reduce.3.23.130521' from the first USER MOD line, or None."""
    m = re.match(r'USER  MOD (reduce\S*)', stdout or '')
    return m.group(1) if m else None


def decisions(stdout):
    """reduce's choice for each His, Asn and Gln, from its USER MOD lines:
    [{residue: 'B:HIS:214:_', set: 'Set 1.1' | 'Single', decision: 'no HE2',
      score: 0.828, category: 'K' | 'C' | 'X' | 'F' | None, detail: 'o=2.4,f=-1'}]
    plus {raw: line} for the NoAdj-H lines of His, Asn or Gln that reduce did
    not adjust (e.g. a His bound to a metal). Chains are as reduce saw them."""
    out = []
    for line in (stdout or '').split('\n'):
        if not line.startswith('USER  MOD'):
            continue
        m = _DECISION.match(line)
        if m:
            icode = m.group('icode').strip() or '_'
            try:
                score = float(m.group('score'))
            except ValueError:
                score = None
            out.append({'residue': '%s:%s:%s:%s' % (m.group('chain'), m.group('resname'),
                                                     m.group('resnum').strip(), icode),
                        'set': ' '.join(m.group('set').replace('Set', 'Set ').split()),
                        'decision': m.group('decision').strip(),
                        'score': score,
                        'category': m.group('category'),
                        'detail': m.group('detail')})
        elif line.startswith('USER  MOD NoAdj-H') and re.search(r' (HIS|ASN|GLN) ', line):
            out.append({'raw': line})
    return out


def chain_renames(taken):
    """{new chain id: {'segment': s, 'chain': c}} for the chains reduce saw
    under another id (copies of a chain in an assembly)."""
    out = {}
    for new, old in (taken or {}).items():
        seg, chid = old.split(' ', 1) if ' ' in old else ('', old)
        if new != chid:
            out[new] = {'segment': seg, 'chain': chid}
    return out
