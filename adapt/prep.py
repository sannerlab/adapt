"""Writing the receptor and ligand PDBQT files: AGFRgui's two save buttons.

Same calls, in the same order, as PDBAnalyserWidget._processReceptor +
makeAndSaveRecFile and makeAndSaveLigFile (pdbAnalyser.py, adcpsuite ADFR):

  receptor: runAGFR.processBM() (MSE->MET, mutations and side-chain
            completion with the best rotamer, altlocs cleared), then
            processReceptor(mol.select('receptor'), charges) (reduce, OpenBabel
            Gasteiger charges, AD types, non-polar H merged), then
            saveReceptorAsPDBQT (REMARK 850 AGFR provenance + atoms).
  ligand:   processBM(), then processLigands(mol.select('ligand')), one
            ADLigand per connected fragment, written with getPDBQTlines.

File names are the ones AGFRgui proposes: <pdbid>_rec.pdbqt, <pdbid>_lig.pdbqt,
or <pdbid>_lig_NNN.pdbqt when the ligand atoms form several fragments. Each
save goes into its own numbered folder of the session so no earlier output is
ever overwritten.
"""
import os

import prody

from . import flags


class Nothing(Exception):
    """No atoms are flagged for what was asked."""


def _processed_bm(agfr):
    mol = agfr._BMProcessed
    if mol is None:
        mol = agfr.processBM()
    return mol


def _write_lines(lines, filename):
    with open(filename, 'w') as f:
        for l in lines:
            f.write('%s\n' % l)


def receptor(session, folder):
    agfr = session.agfr
    mol = _processed_bm(agfr)
    atoms = mol.select('receptor')
    if atoms is None:
        raise Nothing('no atoms selected for the receptor')
    charges = flags.charges(session.tree())
    session.profile._customCharges = charges        # PDBAnalyserWidget.getCharges
    processed = agfr.processReceptor(atoms, charges=charges)
    path = os.path.join(folder, session.profile._pdbid + '_rec.pdbqt')
    agfr.saveReceptorAsPDBQT(processed, path)
    agfr._processedReceptor = processed
    return [path]


def ligands(session, folder):
    agfr = session.agfr
    mol = _processed_bm(agfr)
    ligAtoms = mol.select('ligand')
    if ligAtoms is None:
        raise Nothing('no atoms selected for the ligand')
    nfrag = len(prody.findFragments(ligAtoms))
    if nfrag > 1:
        print('INFO: ligand atoms form %d connected fragments. '
              'We will create a file for each fragment' % nfrag)
    agfr._processedLigands = agfr.processLigands(ligAtoms)
    filepath = os.path.join(folder, session.profile._pdbid + '_lig.pdbqt')
    paths = []
    if len(agfr._adLigs) > 1:
        basename = os.path.splitext(filepath)[0]
        for num, adlig in enumerate(agfr._adLigs):
            filename = '%s_%03d.pdbqt' % (basename, num + 1)
            _write_lines(adlig.getPDBQTlines(), filename)
            agfr._processedLigands[num].filename = filename
            paths.append(filename)
    else:
        _write_lines(agfr._adLigs[0].getPDBQTlines(), filepath)
        agfr._processedLigands[0].filename = filepath
        paths.append(filepath)
    return paths
