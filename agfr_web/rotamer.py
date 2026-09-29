"""Rotamer preview for side chains AGFRgui will complete.

Port of ReceptorItemsTree.addRotamerMol / _setRotamer (pdbAnalyser.py,
adcpsuite copy, including its guard for residues missing N, CA, or both CB
and C).

FAITHFULNESS NOTE: in AGFRgui the rotamer combo box only changes what is
DRAWN. runAGFR.mutateResidues always builds the best-scoring rotamer when the
receptor is written, whatever the combo box shows. This preview does the same
and the page says so; making the choice stick would be a new feature.
"""
import io

import prody
from MolKit2.AARotamer import RotamerMol


class Preview(object):
    __slots__ = ('rotmol', 'res_ca', 'best', 'scores', 'favorable', 'clashes',
                 'current')

    @classmethod
    def make(cls, tree, res, chid, resname, resnum, icode):
        if resname == 'ALA':
            return None     # e.g. 2p4n A:ALA265 with CB missing
        if resname not in RotamerMol._residueNames:
            return None     # RotamerMol asserts on these
        if not res.select('bb and name N'):
            print('atom N missing in %s, cannot use rotamers to mutate' % str(res))
            return None
        if not res.select('bb and name CA'):
            print('atom CA missing in %s, cannot use rotamers to mutate' % str(res))
            return None
        if not (res.select('name CB') or res.select('name C')):
            print('atom CB and C are missing in %s, cannot use rotamers to mutate' % str(res))
            return None
        self = cls()
        rotmol = self.rotmol = RotamerMol(resname)
        rotmol._mol.name = 'rotamer for %s:%s:%s' % (chid, resname, resnum)
        rotmol.rotamer.alignRotToResBB(res.select('bb or name CB'))
        segment = res.getSegnames()[0]
        colliders = tree._ag.select(
            'not deleted and not water and (not (segment %s chid %s resnum `%s` icode %s) '
            'or (segment %s chid %s and resnum `%s` icode %s name N O))' % (
                segment, chid, resnum, icode, segment, chid, resnum, icode))
        best, scores, favorable, clashes = rotmol.rotamer.scoreRotamers(colliders)
        self.best = int(best)
        self.scores = [float(s) for s in scores]
        self.favorable = [len(f) for f in favorable]
        self.clashes = [len(c) for c in clashes]
        self.res_ca = res.select('ca').getCoords()
        self.current = -1    # "auto"
        return self

    def n(self):
        return len(self.rotmol.rotamer.angleList)

    def describe(self):
        return {'auto': self.best, 'n': self.n(), 'current': self.current}

    def index(self, choice):
        """combo-box entry -> rotamer index; entry 0 is "auto <best+1>"."""
        return self.best if choice < 0 else int(choice)

    def message(self, index):
        return 'Score: %.1f, contacts: %d favorable %d clashes' % (
            self.scores[index], self.favorable[index], self.clashes[index])

    def pdb(self, choice):
        """PDB text of the side chain (plus CA) placed at the chosen rotamer."""
        index = self.index(choice)
        ag = self.rotmol._mol._ag
        sc = ag.select('sc').getIndices()
        coords = ag.getCoords()
        coords[sc] = self.rotmol.rotamer.getCoordsForRotamer(index)[sc]
        ag.setCoords(coords)
        ag.select('ca').setCoords(self.res_ca)
        out = io.StringIO()
        prody.writePDBStream(out, ag.select('not name N C O'))
        return out.getvalue(), index
