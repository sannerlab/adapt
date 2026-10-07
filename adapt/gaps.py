# Author: Michel Sanner, PhD. Copyright (c) 2026. Derived from ADFR / AGFRgui code; see THIRD_PARTY_NOTICES.md.
"""Gap detection for AGFR's PDB panel: PDBProfiler._gaps with its defects fixed.

A verbatim copy of ADFR/recFromPDB.py:PDBProfiler._gaps (adcpsuite ADFR,
adfr-0.0.22), extracted from the installed file, with six marked fixes:

  FIX 1  the test for missing residues never fired, so _gapsPerChain was always
         empty and neither AGFRgui nor ADapt ever showed "missing residues";
  FIX 2  once it fires, a chain with no missing C-terminus reported its whole
         sequence as one (aln[0][-0:] is the whole string);
  FIX 3  a protein chain without a SEQRES polymer raised KeyError;
  FIX 4  modified residues and caps written as HETATM (MLZ, ALY, NH2, ACE, ...)
         were left out of the structure sequence and reported as missing;
  FIX 5  a missing terminal residue of a homopolymer run (a His tag) could be
         aligned inside the run, and was then dropped as "no physical gap".
  FIX 6  a residue letter the substitution matrix lacks made the aligner raise
         SystemError, so the file could not be opened at all: '?', which
         AAnames gives MSE (selenomethionine), HID/HIE/HIP, CYX, ACE, ...; and
         U, O, J, which SEQRES gives SEC, PYL, XLE. They are aligned as 'X'.

Checked against the REMARK 465 records of ADFR's ten test PDBs
(tests/test_gaps.py). KEPT AS ADFR HAS IT: a stretch listed as missing whose
flanking residues' C and N are within 3.5 A is not reported (ADFR's own
comment cites 1jff GLY34-HIS61: the chain is continuous there, 2.68 A).

install() replaces the method on the class inside the ADapt process only. The
installed ADFR is not modified: it is shared with agfr and the ADCP runs.
Only the GUI reads _gapsPerChain (the tree's "missing residues" node and the
3D gap markers); nothing that writes a PDBQT does, so the files ADapt writes
are unchanged (tests/test_parity.py).
"""
from Bio import pairwise2 as pw          # noqa: F401  (used by the copied code)
import numpy                             # noqa: F401

from ADFR import recFromPDB
from ADFR.recFromPDB import AAnames, blosum62   # noqa: F401

_ORIGINAL = recFromPDB.PDBProfiler._gaps

# A residue is part of the chain's backbone when one of its atoms is within this
# distance of the neighbouring residue's C or N. A peptide C-N bond is 1.33 A;
# non-bonded contacts are longer than 2.5 A, so nothing else passes.
LINK = 2.0


def _near(res, atom):
    if atom is None:
        return False
    d = res.getCoords() - atom.getCoords()[0]
    return bool(((d * d).sum(1)).min() <= LINK * LINK)


def _polymer_residues(chain):
    """[(residue, one-letter code)] of the chain's polymer, in file order.

    Standard ('protein') residues as ADFR codes them; also any other non-water
    residue bonded into the backbone (a HETATM modified residue or cap), coded
    'X' as SEQRES codes it. A ligand or ion is never bonded to a backbone C or
    N, so it is not taken.
    """
    allres = list(chain.iterResidues())
    prot = [r.select('protein') is not None for r in allres]
    out = []
    for k, res in enumerate(allres):
        if prot[k]:
            out.append((res, AAnames.get(res.getResname(), '*')))  # for blosum62
            continue
        if res.select('water') is not None:
            continue
        prev = out[-1][0] if out else None
        nxt = next((allres[m] for m in range(k + 1, len(allres)) if prot[m]), None)
        if (prev is not None and _near(res, prev.select('name C'))) or \
                (nxt is not None and _near(res, nxt.select('name N'))):
            out.append((res, 'X'))
    return out


_ALPHABET = frozenset(blosum62.alphabet)


def _alignable(seq):
    """FIX 6: `seq` with every letter blosum62 lacks written as 'X'."""
    return ''.join(c if c in _ALPHABET else 'X' for c in seq)


def _with_letters(aln, seqA, seqB):
    """The alignment `aln` of _alignable(seqA) and _alignable(seqB), spelled
    with the original letters (the gaps, '-', are where they were)."""
    def spell(aligned, original):
        it = iter(original)
        return ''.join(c if c == '-' else next(it) for c in aligned)
    a, b = spell(aln[0], seqA), spell(aln[1], seqB)
    if hasattr(aln, '_replace'):
        return aln._replace(seqA=a, seqB=b)
    return (a, b) + tuple(aln[2:])


def _terminal_runs(aln):
    """Move structure gaps ('-' in aln[1]) to the chain-end side of a run of
    identical residues that reaches the start or the end of the alignment.

    Inside such a run every position of the gap scores the same, so the choice
    is the aligner's; at the run's end it is a missing terminus, which is what
    it is. Letters keep their order, so residues still map to the structure.
    """
    a, b = aln[0], list(aln[1])

    def settle(lo, hi, gaps_last):
        seg = b[lo:hi + 1]
        c = a[hi] if gaps_last else a[lo]
        if '-' in seg and all(x in (c, '-') for x in seg):
            k = sum(1 for x in seg if x == c)
            b[lo:hi + 1] = ([c] * k + ['-'] * (len(seg) - k)) if gaps_last \
                else (['-'] * (len(seg) - k) + [c] * k)

    n = len(a)
    if n and a[-1] != '-':
        i = n - 1
        while i > 0 and a[i - 1] == a[-1]:
            i -= 1
        settle(i, n - 1, True)
    if n and a[0] != '-':
        j = 0
        while j < n - 1 and a[j + 1] == a[0]:
            j += 1
        settle(0, j, False)
    seqB = ''.join(b)
    if hasattr(aln, '_replace'):
        return aln._replace(seqB=seqB)
    return (aln[0], seqB) + tuple(aln[2:])


def _gaps(self, polymers):
    # find gaps: for each gap we will have the residue before and after the gap
    # and the list of missing residues inbetween.
    # missing N- and C-termini are dealt with separately
    #from time import time
    #t0 = time()
    ag = self._receptor._ag
    self._gapsPerChain = {}
    if polymers:
        for chain in ag.getHierView().iterChains():
            chid = chain.getChid()
            if chid not in polymers:
                continue  # FIX 3: a protein chain without SEQRES has nothing to align to (KeyError)
            fasta = ''
            reslist = []
            sel = chain.select('protein')
            if sel is None: continue # chain can be only water
            # FIX 4: the structure sequence was built from prody's 'protein' residues
            # only, so a modified residue or cap written as HETATM (MLZ, ALY, NH2,
            # ACE, ...) looked absent and was reported as a missing residue although
            # its atoms are in the file. Residues bonded into the backbone are kept,
            # as 'X' (what SEQRES uses for them).
            gaps = [] # list of [[pdbRes1, restyp, ... restyp, pdbRes2], ...]
            for res, letter in _polymer_residues(chain):
                fasta += letter
                reslist.append(res)

            #aln = pairwise2.align.globalms(
            #    polymers[chid].sequence, fasta, one_alignment_only=True)[0]
            #if aln[1].count('-'): # there are missing residues

            ## biopython alignment changed see issue #3303
            #aln = pw.align.globalms(polymers[chid].sequence, fasta,
            #                        1, -1, -2, -.5)[0]
            # identical to: aln = pw.align.globalms("EPQYEEIPIYL", "EPQ?EEIPIYL", match=1, mismatch=-1, open=-2, extend=-.5)[0]
            seqres = polymers[chid].sequence
            aln = pw.align.globalds(_alignable(seqres), _alignable(fasta),   # FIX 6
                                    blosum62, -10, -0.5)[0]
            aln = _with_letters(aln, seqres, fasta)
            # FIX 5: in a run of identical residues that reaches a chain end (a His
            # tag), the aligner may put a missing terminal residue inside the run;
            # the physical-gap test below then finds two touching residues and the
            # missing residue is dropped. Move such gaps to the end of the run.
            aln = _terminal_runs(aln)
            # FIX 1: ADFR tests format_alignment(*aln).split('\n').count('-'), which counts
            # LINES equal to '-' and so is always 0: no gap was ever reported. A residue
            # missing from the structure is a '-' in the aligned structure sequence.
            if aln[1].count('-'): # there are missing residues

                i = 0
                ri = i # index in pdb sequence
                while aln[1][i]=='-': i+=1
                missingNterm = i
                j = 0
                while aln[1][-j-1]=='-': j+=1
                missingCterm = j
                #ri = i # index in pdb sequence
                gap = []
                for i in range(missingNterm, len(aln[1])-missingCterm):
                    if aln[1][i]=='-':
                        if len(gap)==0:
                            gap = [reslist[ri-1]]
                        gap.append(aln[0][i])
                    else:
                        if len(gap):
                            gap.append(reslist[ri])
                            # check that there is a physical gap (i.e. 1jff GLY34-HIS61 covalent bond)
                            c = gap[0].select('name C')
                            n = gap[-1].select('name N')
                            cc = None
                            if n and c:
                                for b in ag[int(c.getIndices()[0])].iterBonded():
                                    if b.getIndex()==n.getIndices()[0]: # C and N are bonded, i.e. no gap in structure
                                        cc = c.getCoords()[0]
                                        break
                                cc = c.getCoords()[0]
                                cn = n.getCoords()[0]
                                cutoff = 12.25 # 3.5**2. i.e. vdw(C) + vdw(N)
                            else: ## either C or N was not found, so we look at CA-CA distance
                                ca1 = gap[0].select('name CA')
                                ca2 = gap[-1].select('name CA')
                                if ca1 and ca2:
                                    cc = ca1.getCoords()[0]
                                    cn = ca2.getCoords()[0]
                                    cutoff = 16 # 4**2
                                else: # jsut pick first atom in each of the 2 residues
                                    ca1 = gap[0].all[0]
                                    ca2 = gap[-1].all[0]
                                    cc = ca1.getCoords()[0]

                            if cc is not None:
                                d2 = numpy.sum((cn-cc)*(cn-cc))
                                if d2>cutoff:
                                    # get the gap residues name, num and icode
                                    gapData = [[gap[0].getResname(), gap[0].getResnum(), gap[0].getIcode()]]
                                    gapData.extend(gap[1:-1])
                                    gapData.append([gap[-1].getResname(), gap[-1].getResnum(), gap[-1].getIcode()])
                                    gaps.append(gapData)
                            gap = []
                        ri += 1
                mrn = aln[0][:missingNterm] # sequence of missing Nterminus
                # FIX 2: aln[0][-0:] is the WHOLE sequence, so with no missing C-terminus
                # every chain reported its full sequence as a missing C-terminus
                mrc = aln[0][len(aln[0]) - missingCterm:] # sequence of missing Cterminus
                if len(mrn):
                    rr = reslist[0]
                    mrn = [[rr.getResname(), rr.getResnum(), rr.getIcode()], mrn]
                if len(mrc):
                    rr = reslist[-1]
                    mrc = [[rr.getResname(), rr.getResnum(), rr.getIcode()], mrc]
            else:
                mrn = mrc = []
                gaps = []
            if len(mrn) or len(gaps) or len(mrc):
                self._gapsPerChain[chid] = [mrn, gaps, mrc]


# ---- what the 3D view draws for them -------------------------------------
BOND = 1.9     # A: neighbours of a terminal atom, found by distance (no bond table needed)


def _short(seq):
    """gridGUI's label rule: more than 10 residues -> first 3 ... last 3."""
    return seq if len(seq) <= 10 else seq[:3] + '...' + seq[-3:]


def _res(ag, chid, num, icode):
    sel = 'segment A and chid %s resnum `%d`' % (chid if chid != ' ' else '_', num)
    if icode:
        sel += ' icode %s' % icode
    return ag.select(sel)


def geometry(gaps_per_chain, mol):
    """Port of GridGUI.addGeomsForGaps (gridGUI.py): for each gap a dotted line
    from the C of the residue before it to the N of the residue after it, with
    a sphere at its middle labelled with the missing residues; for a missing
    terminus a sphere 3 A beyond the terminal N (or C) atom, opposite its bonds,
    joined to it by a dotted line. Returns [{'line': [p, q], 'sphere': c,
    'label': str, 'kind': 'gap'|'N'|'C'}], coordinates as lists of floats."""
    import numpy as np
    ag = mol._ag
    out = []
    for chid, chgaps in gaps_per_chain.items():
        for gap in chgaps[1]:
            r1 = _res(ag, chid, *gap[0][1:])
            r2 = _res(ag, chid, *gap[-1][1:])
            if not r1 or not r2:
                continue
            c, n = r1.select('name C'), r2.select('name N')
            if c is None or n is None:
                continue
            p, q = c.getCoords()[0], n.getCoords()[0]
            out.append({'kind': 'gap', 'line': [p.tolist(), q.tolist()],
                        'sphere': ((p + q) / 2.).tolist(),
                        'label': _short(''.join(gap[1:-1]))})
        for term, atname in ((chgaps[0], 'N'), (chgaps[-1], 'C')):
            if not len(term):
                continue
            res = _res(ag, chid, *term[0][1:])
            if not res:
                continue
            atoms = res.select('name %s' % atname) or res.select('name CA')
            if atoms is None:
                continue
            at = atoms.getCoords()[0]
            near = ag.select('segment A and chid %s and within %.2f of pt' % (
                chid if chid != ' ' else '_', BOND), pt=at)
            vecs = [x - at for x in (near.getCoords() if near is not None else [])
                    if 0.1 < np.linalg.norm(x - at)]
            if not vecs:
                continue
            v = np.add.reduce(vecs)
            if np.linalg.norm(v) < 1e-6:
                continue
            sp = at - 3 * v / np.linalg.norm(v)
            out.append({'kind': atname, 'line': [sp.tolist(), at.tolist()],
                        'sphere': sp.tolist(), 'label': _short(term[1])})
    return out


def install():
    """Use the fixed _gaps for every PDBProfiler created from now on."""
    recFromPDB.PDBProfiler._gaps = _gaps


def uninstall():
    recFromPDB.PDBProfiler._gaps = _ORIGINAL
