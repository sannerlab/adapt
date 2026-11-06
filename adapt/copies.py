# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Symmetry copies follow their original chain.

An assembly built by symmetry (REMARK 350) holds copies of chains of the
asymmetric unit, each copy one segment (B, C, ..., the letters repeating after
26). The tree lists the original chains only, so every choice made in it --
receptor checkbox, alternate location, ion charge -- is made on the original's
atoms. mirror() gives each copy the same choices, atom by atom; a fragment
checked as ligand is left out of the copies (one ligand file, the one picked).
Copies of a residue sitting on a symmetry axis (an ion or a water on a 3-fold
axis) fall on the atoms of an earlier copy: they are left out.

The original of a chain is its atoms in the first segment that holds it;
every later segment holding that chain is a copy, mapped by position (ProDy
copies a chain's atoms in their order).
"""
import numpy

TOL = 0.25      # A: a copy's atom this close to the same atom of an earlier copy

ORIGINAL = 'adaptorig'      # flag of the original atoms, for the tree's selections


class CopyMap(object):

    def __init__(self, ag):
        n = len(ag)
        segs, chids = ag.getSegnames(), ag.getChids()
        names, resnums = ag.getNames(), ag.getResnums()
        starts = [0] + [i for i in range(1, n) if segs[i] != segs[i - 1]]
        runs = [numpy.arange(s, e) for s, e in zip(starts, starts[1:] + [n])]
        self.original = numpy.zeros(n, bool)
        self.copies = []                # [(copy atom indices, original atom indices)]
        self.problems = []              # chains of a copy that do not match their original
        instances = {}                  # chid -> [original indices, copy indices, ...]
        for run in runs:
            seen = []
            for c in chids[run]:
                if c not in seen:
                    seen.append(c)
            for c in seen:
                idx = run[chids[run] == c]
                if c not in instances:
                    instances[c] = [idx]
                    self.original[idx] = True
                    continue
                first = instances[c][0]
                if len(idx) == len(first) and (names[idx] == names[first]).all() \
                        and (resnums[idx] == resnums[first]).all():
                    instances[c].append(idx)
                    self.copies.append((idx, first))
                else:
                    self.problems.append((segs[idx[0]], c))
        self.duplicate = self._on_symmetry_axes(ag, instances)

    @staticmethod
    def _on_symmetry_axes(ag, instances):
        """Copy atoms of residues whose every atom lies on the same atom of an
        earlier copy of the chain."""
        dup = numpy.zeros(len(ag), bool)
        xyz, resind = ag.getCoords(), ag.getResindices()
        for inst in instances.values():
            for k in range(1, len(inst)):
                for j in range(k):
                    d = numpy.sqrt(((xyz[inst[k]] - xyz[inst[j]]) ** 2).sum(1))
                    close = d < TOL
                    if not close.any():
                        continue
                    res = resind[inst[k]]
                    for r in numpy.unique(res[close]):
                        mine = res == r
                        if close[mine].all():
                            dup[inst[k][mine]] = True
        return dup


def install(ag):
    """The CopyMap of a molecule, and its ORIGINAL flag set."""
    cmap = CopyMap(ag)
    ag.setFlags(ORIGINAL, cmap.original)
    return cmap


def mirror(ag, cmap):
    """Give every copy its original's receptor flags and charges; leave the
    copies' version of a ligand, and atoms on symmetry axes, out."""
    if not cmap.copies:
        return
    rec = ag.getFlags('receptor').copy()
    lig = ag.getFlags('ligand').copy()
    charges = ag.getCharges()
    q = charges.copy() if charges is not None else None
    for dst, src in cmap.copies:
        rec[dst] = rec[src]
        lig[dst] = False
        if q is not None:
            q[dst] = q[src]
    rec[cmap.duplicate] = False
    lig[cmap.duplicate] = False
    ag.setFlags('receptor', rec)
    ag.setFlags('ligand', lig)
    if q is not None:
        ag.setCharges(q)
