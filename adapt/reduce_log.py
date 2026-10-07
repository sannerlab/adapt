# Author: Michel Sanner, PhD. Copyright (c) 2026. Derived from ADFR / AGFRgui code; see THIRD_PARTY_NOTICES.md.
"""Keeping what reduce wrote when ADFR protonates a receptor.

ADFR/recFromPDB.py:protonate runs `reduce -build -DB ... -` and parses its
output, a PDB file with reduce's USER MOD lines (its choice for every His,
Asn and Gln, among others) and every atom with its hydrogens, into a molecule;
the text itself is dropped. This module holds a verbatim copy of the installed
function (adcpsuite ADFR, adfr-0.0.22), extracted from the file, with ONE
marked change: the output and the chain renames made for reduce (copies of a
chain in an assembly get new chain ids) are kept in LAST. Nothing else in the
function changes, so the receptor PDBQT is the same (tests/test_parity.py).

install() points ADFR's two references to protonate (recFromPDB and runAGFR,
which imported the name) at this copy, inside the ADapt process only; the
installed ADFR is not modified.
"""
import os           # noqa: F401  (names the copied code uses)
import platform     # noqa: F401
import re
import subprocess   # noqa: F401

import numpy        # noqa: F401
import prody        # noqa: F401
from mglutil.util.io import StringIO   # noqa: F401
from MolKit2 import Molecule           # noqa: F401

from ADFR import recFromPDB

_ORIGINAL = recFromPDB.protonate

# what the last call left: {'stdout': reduce's output, 'taken': {chain: "seg chain"}}
LAST = {}


def protonate(atoms, molName):
    """execute reduce to protonate a set of Prody Atoms
    handle segments by temporarily renaming chains to have unique chain names
    The original chain bnames are restored after reduce added hydrogen atoms
    """
    heavy = atoms.select("not hydrogen")

    # since reduce does not like chains with the same name in different segments
    # and creates unbound hydrogen atoms we make sure each chain has its own name
    #allsegs = numpy.unique(added[comp]['_contactingSegments']+[v['pepSeg']])
    allsegs = numpy.unique(heavy.getSegnames())
    taken = {}
    if len(allsegs)>1:
        chids = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
                 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
                 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

        for chid in numpy.unique(atoms.getChids()):
            chids.remove(chid)

        for chain in heavy.getHierView().iterChains():
            chid = chain.getChid()
            seg = chain.getSegname()
            if chid in taken:
                #we need to rename this chain
                if len(chids)==0:
                    # we ran out of chain ids to rename e.e. 5vf1.pdb 6 chains and 24 BIOMAT for 144-mer
                    print('ERROR: protonate, too many chains')
                    return None
                newchid = chids.pop(0)
                chain.setChids(newchid)
                taken[newchid] = "%c %c"%(seg, chid)
            else:
                taken[chid] = "%c %c"%(seg, chid)

    ## run reduce as a sbprocess to avoid having to write a file
    #reducePath = os.path.join(os.path.split(sys.executable)[0], 'reduce')
    from ADFR.utils.MakeGrids import findBinary
    reducePath = findBinary("reduce")
    #print "REDUCE PATH:", reducePath
    system_info = platform.uname()
    if system_info[0] == 'Windows':
        shell=False
    else:
        shell=True

    try:
        #import pdb; pdb.set_trace()
        proc = subprocess.Popen('%s -build -DB %s_wwPDB_het_dict.txt -'%(reducePath,reducePath),
                                stdin=subprocess.PIPE ,
                                stdout=subprocess.PIPE ,
                                stderr=subprocess.PIPE,
                                bufsize=1, shell=shell,
                                encoding="ascii"
                                )
        prody.writePDBStream(proc.stdin, atoms.select('not hydrogen'))
    except IOError:
        print('ERROR: cmd %s failed'%(reducePath+' -build -'))
        return None

    (stdout, stderr) = proc.communicate()
    # ADAPT: keep what reduce wrote (a PDB: USER MOD lines + all atoms with H)
    # and the chain renames made for it; ADFR parses it and drops the text
    LAST.clear()
    LAST.update(stdout=stdout, taken=dict(taken))

    # reduce command failed
    if stdout[:16]!="USER  MOD reduce":
        print("ERROR: output of reduce is wrong")
        return None

    #out = stdout.split('\n')
    #for l in out: print l
    agH = prody.parsePDBStream(StringIO(stdout.split('\n')))
    molh = Molecule(molName, agH)
    heavyH = molh._ag.select("not hydrogen")
    assert len(heavy)==len(heavyH)
    hat = molh._ag.select("hydrogen")
    if hat is None: # reduce failed to protonate e.g. 6ann.pdb
        print("ERROR: no hydrogen atoms after executing reduce")
        return None
    numBonds = [x.numBonds() for x in hat]
    #import pdb;pdb.set_trace()
    #if min(numBonds)==0:
    #    raise RuntimeError('%d H atoms with no bond after reduce'%syst)

    # restore CHIDs
    ## v = taken
    ## atomsets = {}
    ## if len(v):
    ##     # build atoms sets for each chain before resetting any chid to
    ##     # avoid affect selection of atoms later no
    ##     for newchid, old in v.items():
    ##         atomsets[newchid] = molh._ag.select('chid %s'%newchid)
    ##     for newchid, old in v.items():
    ##         oldSeg, oldChid = old.split()
    ##         _atoms = atomsets[newchid]
    ##         _atoms.setSegnames([oldSeg]*len(atoms))
    ##         _atoms.setChids([oldChid]*len(atoms))

    return molh


def install():
    """Use this copy for every receptor protonated from now on."""
    from ADFR.utils import runAGFR
    recFromPDB.protonate = protonate
    runAGFR.protonate = protonate


def uninstall():
    from ADFR.utils import runAGFR
    recFromPDB.protonate = _ORIGINAL
    runAGFR.protonate = _ORIGINAL


# ---- reading reduce's output -------------------------------------------------
_DECISION = re.compile(
    r'^USER  MOD (?P<set>Set\s*\d+\.\d+|Single)\s*:\s?(?P<chain>.)(?P<resnum>[ \d-]{4})'
    r'(?P<icode>.)(?P<resname>HIS|ASN|GLN)\s[^:]*:(?P<decision>[^:]*):sc=\s*(?P<score>\S+)'
    r'(?:\s+(?P<category>[KCXF])\((?P<detail>[^)]*)\))?')


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
