# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Writing the receptor and ligand PDBQT files: AGFRgui's two save buttons.

Same calls, in the same order, as PDBAnalyserWidget._processReceptor +
makeAndSaveRecFile and makeAndSaveLigFile (pdbAnalyser.py, adcpsuite ADFR):

  receptor: runAGFR.processBM() (MSE->MET, mutations and side-chain
            completion with the best rotamer, altlocs cleared), then
            processReceptor(mol.select('receptor'), charges) (reduce, OpenBabel
            Gasteiger charges, AD types, non-polar H merged), then the file:
            ADFR's REMARK 850 AGFR provenance lines and atom lines, as
            runAGFR.saveReceptorAsPDBQT writes them, with ONE line added by
            ADapt: which biomolecule the receptor was built from (ADFR's
            header writer has a commented-out placeholder for it). Next to
            it, <pdbid>_rec_H.pdb: the same header, then what reduce wrote
            (its USER MOD lines and every atom with all of its hydrogens).
  ligand:   processBM(), then processLigands(mol.select('ligand')), one
            ADLigand per connected fragment, written with getPDBQTlines,
            with the same biomolecule line added after its "source" line and
            then the residues the fragment came from (ADFR's ligand file
            labels its chain A whatever the chain was, and names no source).

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
    from . import reduce_log
    reduce_log.LAST.clear()
    processed = agfr.processReceptor(atoms, charges=charges)
    path = os.path.join(folder, session.profile._pdbid + '_rec.pdbqt')
    # runAGFR.saveReceptorAsPDBQT, plus the biomolecule line after "source"
    from ADFR import PDBQTlinesForAtoms
    lines = _with_biomolecule(session.profile.getPDBQTHeaderLines(), session.tree())
    lines.extend(PDBQTlinesForAtoms(processed._ag.all))
    _write_lines(lines, path)
    processed.filename = path
    agfr._processedReceptor = processed
    paths = [path]
    # the receptor as reduce protonated it: all hydrogens, non-polar ones too,
    # and reduce's choices (USER MOD lines), under the same provenance header
    out = reduce_log.LAST.get('stdout')
    if out:
        pdb = os.path.join(folder, session.profile._pdbid + '_rec_H.pdb')
        header = _with_biomolecule(session.profile.getPDBQTHeaderLines(), session.tree())
        with open(pdb, 'w') as fh:
            for l in header:
                fh.write('%s\n' % l)
            fh.write(out if out.endswith('\n') else out + '\n')
        paths.append(pdb)
    return paths


def _with_biomolecule(lines, tree, more=()):
    """ADFR's header lines with the biomolecule line, then `more` lines,
    inserted after the "source" line: 'REMARK 850 AGFR source' in a receptor,
    'REMARK 850 source' in a ligand (ADLigand's own style)."""
    lines = list(lines)
    at = next((i + 1 for i, l in enumerate(lines)
               if l.startswith(('REMARK 850 AGFR source', 'REMARK 850 source'))), None)
    if at is None:      # no source line: after the first REMARK 850, or first
        at = next((i + 1 for i, l in enumerate(lines) if l.startswith('REMARK 850')), 0)
    lines[at:at] = [biomolecule_remark(tree)] + list(more)
    return lines


ORIGIN = 'REMARK 850 AGFR ligand from '
PER_LINE = 4        # residue keys per line, so a long peptide stays readable


def origin_remarks(frag):
    """REMARK 850 AGFR ligand from <chain:name:number:icode> ...: the residues
    of the PDB file a ligand fragment was made of, in file order, keyed as the
    other REMARK 850 AGFR lines are ('_' for no insertion code)."""
    from ADFR.recFromPDB import getResidueKey
    keys = [getResidueKey(r) for r in frag.getHierView().iterResidues()]
    return [ORIGIN + ' '.join(keys[i:i + PER_LINE]) for i in range(0, len(keys), PER_LINE)]


def biomolecule_remark(tree):
    """REMARK 850 AGFR biomolecule <tab>: chains <ids>[, <n> copies]

    The tab the receptor was written from: ASU, Biomol N or "ASU and Biomol N",
    the chains of the asymmetric unit it uses, and, for an assembly built by
    symmetry (REMARK 350), how many copies of those chains it holds."""
    chains = ' '.join(sorted(c if c != ' ' else '_' for c in tree._chaindInBM))
    copies = len(set(tree._mol._ag.getSegnames()))
    line = 'REMARK 850 AGFR biomolecule %s: chains %s' % (tree.name, chains)
    if copies > 1:
        line += ', %d copies' % copies
    return line


def ligands(session, folder):
    agfr = session.agfr
    mol = _processed_bm(agfr)
    ligAtoms = mol.select('ligand')
    if ligAtoms is None:
        raise Nothing('no atoms selected for the ligand')
    # processLigands makes one ADLigand per fragment of this same call, in this
    # order, so frags[n] is where ligand file n came from
    frags = list(prody.findFragments(ligAtoms))
    if len(frags) > 1:
        print('INFO: ligand atoms form %d connected fragments. '
              'We will create a file for each fragment' % len(frags))
    agfr._processedLigands = agfr.processLigands(ligAtoms)
    if len(agfr._adLigs) != len(frags):
        raise RuntimeError('%d ligands made from %d fragments' % (len(agfr._adLigs), len(frags)))
    filepath = os.path.join(folder, session.profile._pdbid + '_lig.pdbqt')
    paths = []
    if len(agfr._adLigs) > 1:
        basename = os.path.splitext(filepath)[0]
        for num, adlig in enumerate(agfr._adLigs):
            filename = '%s_%03d.pdbqt' % (basename, num + 1)
            _write_lines(_with_biomolecule(adlig.getPDBQTlines(), session.tree(),
                                           origin_remarks(frags[num])), filename)
            agfr._processedLigands[num].filename = filename
            paths.append(filename)
    else:
        _write_lines(_with_biomolecule(agfr._adLigs[0].getPDBQTlines(), session.tree(),
                                       origin_remarks(frags[0])), filepath)
        agfr._processedLigands[0].filename = filepath
        paths.append(filepath)
    return paths
